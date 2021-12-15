// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/TrackFitting/detail/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfActor.hpp"
#include "Acts/Propagator/MultiStepperAborters.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"

#include <fstream>


#define RETURN_ERROR_OR_ABORT_FIT(error) \
  if (options.abortOnError) {            \
    std::abort();                        \
  } else {                               \
    return error;                        \
  }

namespace Acts {

struct GsfOptions {
  std::reference_wrapper<const GeometryContext> geoContext;
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  std::reference_wrapper<const CalibrationContext> calibrationContext;
  
  KalmanFitterExtensions extensions;
  
  LoggerWrapper logger;

  PropagatorPlainOptions propagatorPlainOptions;

  const Surface* referenceSurface = nullptr;

  bool abortOnError = true;

  std::size_t maxComponents = 4;
  
  bool applyMaterialEffects = true;
};

template <typename propagator_t>
struct GaussianSumFitter {
  GaussianSumFitter(propagator_t&& propagator,
                    const std::string& high_x0_bethe_heitler_path,
                    const std::string& low_x0_bethe_heitler_path)
      : m_propagator(std::move(propagator)),
        m_bethe_heitler_approx(
            detail::load_bethe_heitler_data<6, 5>(low_x0_bethe_heitler_path),
            detail::load_bethe_heitler_data<6, 5>(high_x0_bethe_heitler_path)) {
  }

  GaussianSumFitter(propagator_t&& propagator)
      : m_propagator(std::move(propagator)),
        m_bethe_heitler_approx(detail::bh_cdf_cmps6_order5_data) {}

  /// The navigator type
  using GsfNavigator = typename propagator_t::Navigator;

  /// The propagator instance used by the fit function
  propagator_t m_propagator;

  /// The fitter holds the instance of the bethe heitler approx
  detail::BHApprox m_bethe_heitler_approx;

  /// @brief The fit function for the Direct navigator
  template <typename source_link_it_t, typename start_parameters_t>
  Acts::Result<Acts::KalmanFitterResult> fit(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters,
      const GsfOptions& options,
      const std::vector<const Surface*>& sSequence) const {
    // Check if we have the correct navigator
    static_assert(
        std::is_same_v<DirectNavigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [&sSequence, this](const auto& opts,
                                                 const auto& logger) {
      using Actors = ActionList<GsfActor, DirectNavigator::Initializer>;
      using Aborters = AbortList<>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);
      
      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = sSequence;
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [&sSequence, this](const auto& opts,
                                                 const auto& logger) {
      using Actors = ActionList<GsfActor, DirectNavigator::Initializer>;
      using Aborters = AbortList<>;
      
      std::vector<const Surface*> backwardSequence(
          std::next(sSequence.rbegin()), sSequence.rend());
      backwardSequence.push_back(opts.referenceSurface);

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = std::move(backwardSequence);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer);
  }

  /// @brief The fit function for the standard navigator
  template <typename source_link_it_t, typename start_parameters_t>
  Acts::Result<Acts::KalmanFitterResult> fit(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters,
      const GsfOptions& options) const {
    // Check if we have the correct navigator
    static_assert(std::is_same_v<Navigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [this](const auto& opts, const auto& logger) {
      using Actors = ActionList<GsfActor>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);
      propOptions.setPlainOptions(opts.propagatorPlainOptions);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [this](const auto& opts, const auto& logger) {
      using Actors = ActionList<GsfActor>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);
      
      propOptions.setPlainOptions(opts.propagatorPlainOptions);
      
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer);
  }

  template <typename source_link_it_t, typename start_parameters_t,
            typename fwd_prop_initializer_t, typename bwd_prop_initializer_t>
  Acts::Result<Acts::KalmanFitterResult> fit_impl(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters,
      const GsfOptions& options,
      const fwd_prop_initializer_t& fwdPropInitializer,
      const bwd_prop_initializer_t& bwdPropInitializer) const {
    // The logger
    const auto& logger = options.logger;

    // Print some infos about the start parameters
    ACTS_VERBOSE("Run Gsf with start parameters: \n" << sParameters);

    auto intersectionStatusStartSurface =
        sParameters.referenceSurface()
            .intersect(GeometryContext{},
                       sParameters.position(GeometryContext{}),
                       sParameters.unitDirection(), true)
            .intersection.status;

    if (intersectionStatusStartSurface != Intersection3D::Status::onSurface) {
      ACTS_ERROR(
          "Surface intersection of start parameters with bound-check failed");
      return GsfError::StartParametersNotOnStartSurface;
    }

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(begin, end)
                              << " input measurements");
    std::map<GeometryIdentifier, std::reference_wrapper<const SourceLink>>
        inputMeasurements;
    for (auto it = begin; it != end; ++it) {
      const SourceLink& sl = *it;
      inputMeasurements.emplace(sl.geometryId(), *it);
    }

    ACTS_VERBOSE(
        "Gsf: Final measuerement map size: " << inputMeasurements.size());
    throw_assert(sParameters.covariance() != std::nullopt,
                 "we need a covariance here...");

    /////////////////
    // Forward pass
    /////////////////
    ACTS_VERBOSE("+-----------------------------+");
    ACTS_VERBOSE("| Gsf: Do forward propagation |");
    ACTS_VERBOSE("+-----------------------------+");

    auto fwdResult = [&]() {
      MultiComponentBoundTrackParameters<SinglyCharged> params(
          sParameters.referenceSurface().getSharedPtr(),
          sParameters.parameters(), sParameters.covariance());

      auto fwdPropOptions = fwdPropInitializer(options, logger);

      // Catch the actor and set the measurements
      auto& actor = fwdPropOptions.actionList.template get<GsfActor>();
      actor.m_cfg.inputMeasurements = inputMeasurements;
      actor.m_cfg.maxComponents = options.maxComponents;
      actor.m_cfg.extensions = options.extensions;
      actor.m_cfg.abortOnError = options.abortOnError;
      actor.m_cfg.applyMaterialEffects = options.applyMaterialEffects;

      fwdPropOptions.direction = Acts::forward;

      return m_propagator.propagate(params, fwdPropOptions);
    }();

    if (!fwdResult.ok()) {
      RETURN_ERROR_OR_ABORT_FIT(fwdResult.error());
    }

    auto& fwdGsfResult = (*fwdResult).template get<GsfResult>();

    if (!fwdGsfResult.result.ok()) {
      RETURN_ERROR_OR_ABORT_FIT(fwdGsfResult.result.error());
    }

    if (fwdGsfResult.processedStates == 0) {
      RETURN_ERROR_OR_ABORT_FIT(GsfError::NoStatesCreated);
    }

    ACTS_VERBOSE("Finished forward propagation");
    ACTS_VERBOSE("- visited surfaces: " << fwdGsfResult.visitedSurfaces.size());
    ACTS_VERBOSE("- processed states: " << fwdGsfResult.processedStates);
    ACTS_VERBOSE("- measuerement states: " << fwdGsfResult.measurementStates);

    //////////////////
    // Backward pass
    //////////////////
    ACTS_VERBOSE("+------------------------------+");
    ACTS_VERBOSE("| Gsf: Do backward propagation |");
    ACTS_VERBOSE("+------------------------------+");

    auto bwdResult = [&]() {
      // Use last forward state as start parameters for backward propagation
      const auto params = detail::extractMultiComponentState(
          fwdGsfResult.fittedStates, fwdGsfResult.currentTips,
          fwdGsfResult.weightsOfStates, detail::StatesType::eFiltered);

      auto bwdPropOptions = bwdPropInitializer(options, logger);

      auto& actor = bwdPropOptions.actionList.template get<GsfActor>();
      actor.m_cfg.inputMeasurements = inputMeasurements;
      actor.m_cfg.maxComponents = options.maxComponents;
      actor.m_cfg.abortOnError = options.abortOnError;
      actor.m_cfg.applyMaterialEffects = options.applyMaterialEffects;
      actor.m_cfg.extensions = options.extensions;

      // Workaround to get the first state into the MultiTrajectory seems also
      // to be necessary for standard navigator to prevent double kalman
      // update on the last surface
      actor.m_cfg.multiTrajectoryInitializer = [&fwdGsfResult](
                                                   auto& result,
                                                   [[maybe_unused]] const auto& loggerr) {
        result.currentTips.clear();
        
        throw std::runtime_error("fix this -Wshadow with logger and gcc somehow");

        // ACTS_VERBOSE(
        //     "Initialize the MultiTrajectory with information provided to the "
        //     "Actor");

        for (const auto idx : fwdGsfResult.currentTips) {
          result.currentTips.push_back(
              result.fittedStates.addTrackState(TrackStatePropMask::All));
          auto proxy =
              result.fittedStates.getTrackState(result.currentTips.back());
          proxy.copyFrom(fwdGsfResult.fittedStates.getTrackState(idx));
          result.weightsOfStates[result.currentTips.back()] =
              fwdGsfResult.weightsOfStates.at(idx);

          // Because we are backwards, we use forward filtered as predicted
          proxy.predicted() = proxy.filtered();
          proxy.predictedCovariance() = proxy.filteredCovariance();

          // Mark surface as visited
          result.visitedSurfaces.insert(proxy.referenceSurface().geometryId());
        }
      };

      bwdPropOptions.direction = Acts::backward;

      // const auto targetSurface = [&]() {
      //   const Surface* ts = nullptr;
      //   fwdGsfResult.fittedStates.visitBackwards(
      //       fwdGsfResult.currentTips.front(),
      //       [&ts](const auto& state) { ts = &state.referenceSurface();
      //       });
      //   throw_assert(ts, "target surface must not be nullptr");
      //   return ts;
      // }();

      // TODO somehow this proagation fails if we target the first
      // measuerement surface, go instead back to beamline for now
      return m_propagator
          .template propagate<decltype(params), decltype(bwdPropOptions),
                              MultiStepperSurfaceReached>(
              params, *options.referenceSurface, bwdPropOptions);
    }();

    if (!bwdResult.ok()) {
      RETURN_ERROR_OR_ABORT_FIT(bwdResult.error());
    }

    auto& bwdGsfResult = (*bwdResult).template get<GsfResult>();

    if (!bwdGsfResult.result.ok()) {
      RETURN_ERROR_OR_ABORT_FIT(bwdGsfResult.result.error());
    }

    if (bwdGsfResult.processedStates == 0) {
      RETURN_ERROR_OR_ABORT_FIT(GsfError::NoStatesCreated);
    }

    ////////////////////////////////////
    // Smooth and create Kalman Result
    ////////////////////////////////////
    ACTS_VERBOSE("Gsf: Do smoothing");
    ACTS_VERBOSE(
        "- Fwd measurement states: " << fwdGsfResult.measurementStates);
    ACTS_VERBOSE(
        "- Bwd measurement states: " << bwdGsfResult.measurementStates);

    const auto smoothResult = detail::smoothAndCombineTrajectories<true>(
        fwdGsfResult.fittedStates, fwdGsfResult.currentTips,
        fwdGsfResult.weightsOfStates, bwdGsfResult.fittedStates,
        bwdGsfResult.currentTips, bwdGsfResult.weightsOfStates, logger);

    // Cannot use structured binding since they cannot be captured in lambda
    const auto& combinedTraj = std::get<0>(smoothResult);
    const auto lastTip = std::get<1>(smoothResult);

    // Some test
    if (lastTip == SIZE_MAX) {
      RETURN_ERROR_OR_ABORT_FIT(GsfError::NoStatesCreated);
    }

    Acts::KalmanFitterResult kalmanResult;
    kalmanResult.lastTrackIndex = lastTip;
    kalmanResult.fittedStates = std::move(combinedTraj);
    kalmanResult.smoothed = true;
    kalmanResult.reversed = true;
    kalmanResult.finished = true;
    kalmanResult.lastMeasurementIndex = lastTip;

    ///////////////////////////////////////////////////////
    // Propagate back to origin with smoothed parameters //
    ///////////////////////////////////////////////////////
    ACTS_VERBOSE("+--------------------------------------+");
    ACTS_VERBOSE("| Gsf: Do propagation back to beamline |");
    ACTS_VERBOSE("+--------------------------------------+");
    auto lastResult = [&]() -> Result<std::unique_ptr<BoundTrackParameters>> {
      const auto& [surface, lastSmoothedState] =
          std::get<2>(smoothResult).front();

      throw_assert(
          detail::weightsAreNormalized(
              lastSmoothedState,
              [](const auto& tuple) { return std::get<double>(tuple); }),
          "");

      const MultiComponentBoundTrackParameters<SinglyCharged> params(
          surface->getSharedPtr(), lastSmoothedState);

      auto lastPropOptions = bwdPropInitializer(options, logger);

      auto& actor = lastPropOptions.actionList.template get<GsfActor>();
      actor.m_cfg.maxComponents = options.maxComponents;
      actor.m_cfg.abortOnError = options.abortOnError;
      actor.m_cfg.applyMaterialEffects = options.applyMaterialEffects;
      actor.m_cfg.surfacesToSkip.insert(surface->geometryId());

      lastPropOptions.direction = Acts::backward;

      auto result =
          m_propagator
              .template propagate<decltype(params), decltype(lastPropOptions),
                                  MultiStepperSurfaceReached>(
                  params, *options.referenceSurface, lastPropOptions);

      if (!result.ok()) {
        return result.error();
      } else {
        return std::move((*result).endParameters);
      }
    }();

    if (!lastResult.ok()) {
      RETURN_ERROR_OR_ABORT_FIT(lastResult.error());
    }

    kalmanResult.fittedParameters = **lastResult;

    return kalmanResult;
  }
};

}  // namespace Acts
