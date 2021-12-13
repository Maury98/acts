// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"

using namespace Acts;
using namespace Acts::VectorHelpers;

/////////////////////////////////////////////////////
// Some useful global objects, typedefs and structs
/////////////////////////////////////////////////////
const MagneticFieldContext magCtx;
const GeometryContext geoCtx;

using MultiStepperLoop =
    MultiEigenStepperLoop<StepperExtensionList<DefaultExtension>>;
using SingleStepper = EigenStepper<StepperExtensionList<DefaultExtension>>;

struct Options {
  double tolerance = 1e-4;
  double stepSizeCutOff = 0.0;
  std::size_t maxRungeKuttaStepTrials = 10;
  double mass = 1.0;
};

struct Navigation {};

template <typename stepper_state_t>
struct DummyPropState {
  stepper_state_t &stepping;
  Options options;
  Navigation navigation;
  GeometryContext geoContext;

  DummyPropState(stepper_state_t &ss)
      : stepping(ss),
        options(Options{}),
        navigation(Navigation{}),
        geoContext(geoCtx) {}
};

template <typename T>
using components_t = typename T::components;

//////////////////////////////////////////////////////
/// Test the construction of the MultiStepper::State
//////////////////////////////////////////////////////
template <typename multi_stepper_t, typename charge_t, bool Cov>
void test_multi_stepper_state() {
  static_assert(std::is_same_v<charge_t, SinglyCharged> ||
                std::is_same_v<charge_t, Neutral>);

  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  using cov_t = std::optional<BoundSymMatrix>;
  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, BoundVector::Ones(),
               Cov ? cov_t{BoundSymMatrix::Identity()} : cov_t{std::nullopt}});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  // Test charged parameters without covariance matrix
  MultiComponentBoundTrackParameters<charge_t> multi_pars(surface, cmps);

  MultiState state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                   tolerance);

  MultiStepper ms((bField));

  BOOST_CHECK_EQUAL(cmps.size(), ms.numberComponents(state));

  // Test the result & compare with the input/test for reasonable members
  auto const_iterable = ms.constComponentIterable(state);
  for (const auto &cmp : const_iterable) {
    BOOST_CHECK_EQUAL(cmp.jacTransport(), FreeMatrix::Identity());
    BOOST_CHECK_EQUAL(cmp.derivative(), FreeVector::Zero());
    if constexpr (not Cov) {
      BOOST_CHECK_EQUAL(cmp.jacToGlobal(), BoundToFreeMatrix::Zero());
      BOOST_CHECK_EQUAL(cmp.cov(), BoundSymMatrix::Zero());
    }
  }

  const auto expected_charge = []() {
    if constexpr (std::is_same_v<charge_t, Neutral>) {
      return 0.0;
    } else {
      return 1.0;
    }
  }();

  BOOST_CHECK_EQUAL(ms.charge(state), expected_charge);
  for (const auto &cmp : const_iterable) {
    BOOST_CHECK_EQUAL(cmp.charge(), expected_charge);
  }

  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  for (const auto &cmp : const_iterable) {
    BOOST_CHECK_EQUAL(cmp.pathAccumulated(), 0.);
  }

  // navDir and the covTransport in the MultiEigenStepperLoop are redundant and
  // thus not part of the interface. However, we want to check them for
  // consistency.
  if constexpr (Acts::Concepts::exists<components_t, MultiState>) {
    BOOST_CHECK(not state.covTransport);
    for (const auto &cmp : state.components) {
      BOOST_CHECK(cmp.state.covTransport == Cov);
    }

    BOOST_CHECK_EQUAL(state.navDir, ndir);
    for (const auto &cmp : state.components) {
      BOOST_CHECK_EQUAL(cmp.state.navDir, ndir);
    }
  }
}

BOOST_AUTO_TEST_CASE(multi_stepper_state_charged_no_cov) {
  test_multi_stepper_state<MultiStepperLoop, SinglyCharged, false>();
}

BOOST_AUTO_TEST_CASE(multi_stepper_state_neutral_no_cov) {
  test_multi_stepper_state<MultiStepperLoop, Neutral, false>();
}

BOOST_AUTO_TEST_CASE(multi_stepper_state_charged_cov) {
  test_multi_stepper_state<MultiStepperLoop, SinglyCharged, true>();
}

template <typename multi_stepper_t>
void test_multi_stepper_state_invalid() {
  using MultiState = typename multi_stepper_t::State;

  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  // Empty component vector
  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  BOOST_CHECK_THROW(
      MultiState(geoCtx, magCtx, bField, multi_pars, ndir, stepSize, tolerance),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_invalid) {
  test_multi_stepper_state_invalid<MultiStepperLoop>();
}

////////////////////////////////////////////////////////////////////////
// Compare the Multi-Stepper against the Eigen-Stepper for consistency
////////////////////////////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_multi_stepper_vs_eigen_stepper() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;
  const BoundVector pars = BoundVector::Ones();
  const BoundSymMatrix cov = BoundSymMatrix::Identity();

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  SingleBoundTrackParameters<SinglyCharged> single_pars(surface, pars, cov);

  MultiState multi_state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                         tolerance);
  SingleStepper::State single_state(geoCtx, bField->makeCache(magCtx),
                                    single_pars, ndir, stepSize, tolerance);

  MultiStepper multi_stepper(bField);
  SingleStepper single_stepper(bField);

  // Do some steps and check that the results match
  for (int i = 0; i < 10; ++i) {
    // Single stepper
    auto single_prop_state = DummyPropState(single_state);
    auto single_result = single_stepper.step(single_prop_state);
    single_stepper.transportCovarianceToCurvilinear(single_state);

    // Multi stepper;
    auto multi_prop_state = DummyPropState(multi_state);
    auto multi_result = multi_stepper.step(multi_prop_state);
    multi_stepper.transportCovarianceToCurvilinear(multi_state);

    // Check equality
    BOOST_REQUIRE(multi_result.ok() == true);
    BOOST_REQUIRE(multi_result.ok() == single_result.ok());

    BOOST_CHECK_EQUAL(*single_result, *multi_result);

    for (const auto &cmp : multi_stepper.constComponentIterable(multi_state)) {
      BOOST_CHECK_EQUAL(cmp.pars(), single_state.pars);
      BOOST_CHECK_EQUAL(cmp.charge(), single_state.q);
      BOOST_CHECK_EQUAL(cmp.cov(), single_state.cov);
      BOOST_CHECK_EQUAL(cmp.jacTransport(), single_state.jacTransport);
      BOOST_CHECK_EQUAL(cmp.jacToGlobal(), single_state.jacToGlobal);
      BOOST_CHECK_EQUAL(cmp.derivative(), single_state.derivative);
      BOOST_CHECK_EQUAL(cmp.pathAccumulated(), single_state.pathAccumulated);
    }
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_vs_single_eigen) {
  test_multi_stepper_vs_eigen_stepper<MultiStepperLoop>();
}

////////////////////////////////////////////////////
// Test the modifying accessors to the components
////////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_components_modifying_accessors() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps;
  for (int i = 0; i < 4; ++i) {
    cmps.push_back({0.25, BoundVector::Random(), BoundSymMatrix::Random()});
  }

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  MultiState mutable_multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                 stepSize, tolerance);
  const MultiState const_multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                     stepSize, tolerance);

  MultiStepper multi_stepper(bField);

  auto modify = [&](const auto &projector) {
    // Here test the mutable overloads of the mutable iterable
    for (auto cmp : multi_stepper.componentIterable(mutable_multi_state)) {
      using type = std::decay_t<decltype(projector(cmp))>;
      if constexpr (std::is_enum_v<type>) {
        projector(cmp) =
            static_cast<type>(static_cast<int>(projector(cmp)) + 1);
      } else {
        projector(cmp) *= 2.0;
      }
    }
  };

  auto check = [&](const auto &projector) {
    // Here test the const-member functions of the mutable iterable
    auto mutable_state_iterable =
        multi_stepper.componentIterable(mutable_multi_state);
    // Here test the const iterable
    auto const_state_iterable =
        multi_stepper.constComponentIterable(const_multi_state);

    auto mstate_it = mutable_state_iterable.begin();
    auto cstate_it = const_state_iterable.begin();
    for (; cstate_it != const_state_iterable.end(); ++mstate_it, ++cstate_it) {
      const auto mstate_cmp = *mstate_it;
      auto cstate_cmp = *cstate_it;

      using type = std::decay_t<decltype(projector(mstate_cmp))>;

      if constexpr (std::is_arithmetic_v<type>) {
        BOOST_CHECK_CLOSE(projector(mstate_cmp), 2.0 * projector(cstate_cmp),
                          1.e-8);
      } else if constexpr (std::is_enum_v<type>) {
        BOOST_CHECK_EQUAL(static_cast<int>(projector(mstate_cmp)),
                          1 + static_cast<int>(projector(cstate_cmp)));
      } else {
        BOOST_CHECK(
            projector(mstate_cmp).isApprox(2.0 * projector(cstate_cmp), 1.e-8));
      }
    }
  };

  const auto projectors = std::make_tuple(
      [](auto &cmp) -> decltype(auto) { return cmp.status(); },
      [](auto &cmp) -> decltype(auto) { return cmp.pathAccumulated(); },
      [](auto &cmp) -> decltype(auto) { return cmp.charge(); },
      [](auto &cmp) -> decltype(auto) { return cmp.weight(); },
      [](auto &cmp) -> decltype(auto) { return cmp.pars(); },
      [](auto &cmp) -> decltype(auto) { return cmp.cov(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacTransport(); },
      [](auto &cmp) -> decltype(auto) { return cmp.derivative(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacobian(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacToGlobal(); });

  std::apply(
      [&](const auto &...projs) {
        (
            [&]() {
              modify(projs);
              check(projs);
            }(),
            ...);
      },
      projectors);
}

BOOST_AUTO_TEST_CASE(multi_eigen_component_iterable_with_modification) {
  test_components_modifying_accessors<MultiStepperLoop>();
}

/////////////////////////////////////////////
// Test if the surface status update works
/////////////////////////////////////////////
template <typename multi_stepper_t>
void test_multi_stepper_surface_status_update() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto start_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});

  auto right_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0});

  auto bField = std::make_shared<NullBField>();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(2, {0.5, BoundVector::Zero(), std::nullopt});
  std::get<BoundVector>(cmps[0])[eBoundTheta] = M_PI_2;
  std::get<BoundVector>(cmps[1])[eBoundTheta] = -M_PI_2;
  std::get<BoundVector>(cmps[0])[eBoundQOverP] = 1.0;
  std::get<BoundVector>(cmps[1])[eBoundQOverP] = 1.0;

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(start_surface,
                                                               cmps);

  BOOST_REQUIRE(std::get<1>(multi_pars[0])
                    .unitDirection()
                    .isApprox(Vector3{1.0, 0.0, 0.0}, 1.e-10));
  BOOST_REQUIRE(std::get<1>(multi_pars[1])
                    .unitDirection()
                    .isApprox(Vector3{-1.0, 0.0, 0.0}, 1.e-10));

  MultiState multi_state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                         tolerance);
  SingleStepper::State single_state(geoCtx, bField->makeCache(magCtx),
                                    std::get<1>(multi_pars[0]), ndir, stepSize,
                                    tolerance);

  MultiStepper multi_stepper(bField);
  SingleStepper single_stepper(bField);

  // Update surface status and check
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *right_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::reachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::reachable);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::missed);
  }

  // Step forward now
  {
    auto multi_prop_state = DummyPropState(multi_state);
    multi_stepper.step(multi_prop_state);

    // Single stepper
    auto single_prop_state = DummyPropState(single_state);
    single_stepper.step(single_prop_state);
  }

  // Update surface status and check again
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *right_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::onSurface);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::onSurface);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::missed);
  }

  // Start surface should be unreachable
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *start_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::unreachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::unreachable);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::unreachable);
  }
}

BOOST_AUTO_TEST_CASE(test_surface_status_and_cmpwise_bound_state) {
  test_multi_stepper_surface_status_update<MultiStepperLoop>();
}

//////////////////////////////////
// Test Bound state computations
//////////////////////////////////
template <typename multi_stepper_t>
void test_component_bound_state() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto start_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});

  auto right_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0});

  auto bField = std::make_shared<NullBField>();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(2, {0.5, BoundVector::Zero(), std::nullopt});
  std::get<BoundVector>(cmps[0])[eBoundTheta] = M_PI_2;
  std::get<BoundVector>(cmps[1])[eBoundTheta] = -M_PI_2;
  std::get<BoundVector>(cmps[0])[eBoundQOverP] = 1.0;
  std::get<BoundVector>(cmps[1])[eBoundQOverP] = 1.0;

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(start_surface,
                                                               cmps);

  BOOST_REQUIRE(std::get<1>(multi_pars[0])
                    .unitDirection()
                    .isApprox(Vector3{1.0, 0.0, 0.0}, 1.e-10));
  BOOST_REQUIRE(std::get<1>(multi_pars[1])
                    .unitDirection()
                    .isApprox(Vector3{-1.0, 0.0, 0.0}, 1.e-10));

  MultiState multi_state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                         tolerance);
  SingleStepper::State single_state(geoCtx, bField->makeCache(magCtx),
                                    std::get<1>(multi_pars[0]), ndir, stepSize,
                                    tolerance);

  MultiStepper multi_stepper(bField);
  SingleStepper single_stepper(bField);

  // Step forward now
  {
    multi_stepper.updateSurfaceStatus(multi_state, *right_surface, false);
    auto multi_prop_state = DummyPropState(multi_state);
    multi_stepper.step(multi_prop_state);

    // Single stepper
    single_stepper.updateSurfaceStatus(single_state, *right_surface, false);
    auto single_prop_state = DummyPropState(single_state);
    single_stepper.step(single_prop_state);
  }

  // Check component-wise bound-state
  {
    auto single_bound_state =
        single_stepper.boundState(single_state, *right_surface, true);
    BOOST_REQUIRE(single_bound_state.ok());

    auto cmp_iterable = multi_stepper.componentIterable(multi_state);

    auto ok_bound_state =
        (*cmp_iterable.begin()).boundState(*right_surface, true);
    BOOST_REQUIRE(ok_bound_state.ok());
    BOOST_CHECK(*single_bound_state == *ok_bound_state);

    auto failed_bound_state =
        (*(++cmp_iterable.begin())).boundState(*right_surface, true);
    BOOST_CHECK(not failed_bound_state.ok());
  }
}

BOOST_AUTO_TEST_CASE(test_component_wise_bound_state) {
  test_component_bound_state<MultiStepperLoop>();
}

template <typename multi_stepper_t>
void test_combined_bound_state_function() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});
  auto bField = std::make_shared<NullBField>();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSymMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  MultiState multi_state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                         tolerance);
  MultiStepper multi_stepper(bField);

  auto res = multi_stepper.boundState(multi_state, *surface, true);

  BOOST_REQUIRE(res.ok());

  const auto [bound_pars, jacobian, pathLength] = *res;

  BOOST_CHECK(jacobian == decltype(jacobian)::Zero());
  BOOST_CHECK(pathLength == 0.0);
  BOOST_CHECK(bound_pars.parameters().isApprox(pars, 1.e-8));
  BOOST_CHECK(bound_pars.covariance()->isApprox(cov, 1.e-8));
}

BOOST_AUTO_TEST_CASE(test_combined_bound_state) {
  test_combined_bound_state_function<MultiStepperLoop>();
}

//////////////////////////////////////////////////
// Test the combined curvilinear state function
//////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_combined_curvilinear_state_function() {
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});
  auto bField = std::make_shared<NullBField>();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSymMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});
  SingleBoundTrackParameters<SinglyCharged> check_pars(surface, pars, cov);

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  MultiState multi_state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                         tolerance);
  MultiStepper multi_stepper(bField);

  const auto [curv_pars, jac, pathLength] =
      multi_stepper.curvilinearState(multi_state);

  BOOST_CHECK(
      curv_pars.fourPosition(multi_state.geoContext)
          .isApprox(check_pars.fourPosition(multi_state.geoContext), 1.e-8));
  BOOST_CHECK(
      curv_pars.unitDirection().isApprox(check_pars.unitDirection(), 1.e-8));
  BOOST_CHECK_CLOSE(curv_pars.absoluteMomentum(), check_pars.absoluteMomentum(),
                    1.e-8);
  BOOST_CHECK_CLOSE(curv_pars.charge(), check_pars.charge(), 1.e-8);
}

BOOST_AUTO_TEST_CASE(test_curvilinear_state) {
  test_combined_curvilinear_state_function<MultiStepperLoop>();
}

//////////////////////////////////////////////////
// Instatiate a Propagator with the MultiStepper
//////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(propagator_instatiation_test) {
  auto bField = std::make_shared<NullBField>();
  MultiEigenStepperLoop<> multi_stepper(bField);
  [[maybe_unused]] Propagator<MultiEigenStepperLoop<>> propagator(
      multi_stepper);
}
