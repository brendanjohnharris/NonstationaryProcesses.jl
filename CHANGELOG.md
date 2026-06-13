# Changelog


## [v0.1.1]

Dependency and loading migration; no public API changes (the `Process` type,
system constructors, parameter profiles, and exported utilities keep the same
names and signatures).

### Changed
- **Replaced the lazy DifferentialEquations weak dependency** (`Requires` /
  `@require`) with direct dependencies on `OrdinaryDiffEq` and `StochasticDiffEq`.
  OrdinaryDiffEq 7 (the SciMLBase 3 build) no longer bundles RK4, RadauIIA3,
  FunctionMap, or Rodas5, so each solver is pulled in as its own
  `OrdinaryDiffEq*` subpackage.
- System definitions (`ChaoticFlows`, `ChaoticMaps`, `DeterministicFlows`,
  `Stochastic`, `Transforms`) moved from `ext/DifferentialEquationsExt/` into
  `src/` and are now always loaded on `using NonstationaryProcesses`, rather than
  only after `using DifferentialEquations`.
- Reworked CI (scheduled and manual triggers); compat and README updates.

### Added
- `SimulatorTests` covering the simulation pipeline.

### Compatibility notes
- `DifferentialEquations` is no longer a dependency and no longer needs to be
  loaded to access the chaotic systems; downstream code reaching DiffEq solvers
  *through* this package's namespace now sees OrdinaryDiffEq/StochasticDiffEq.

## [v0.1.0] --- 2024-11-26

The package carried a single registered version throughout its early history;
the milestones below are all folded into the `v0.1.0` tag.

### 2023-11 --- Flows and fixes

#### Added
- Dequan Li system.

#### Fixed
- In-place flow definitions.
- Phase-corrupted Lorenz updates; aesthetic restructuring.

### 2022 --- Core/library split

#### Changed
- **Extracted core types and machinery to `NonstationaryProcessesBase.jl`**;
  this package became the system library on top of it.
- Updated for DifferentialEquations v7.
- Lazily load DifferentialEquations (via `@require`) to cut precompilation time.
- Allow passing a single constant parameter to a `Process`.

#### Fixed
- `SVector` import bug; AR process; dependency/compat issues.

### 2021-06 to 2021-11 --- Surrogates, spectra, more systems

#### Added
- Fourier phase-corruption of existing `Process`es to generate nonstationary
  random-phase surrogates with parameterised power spectra; animation support.
- Windowed Fourier transforms (incl. single-window / non-STFT case).
- Greatest-Lyapunov-exponent estimation from a `Process`; standardised plots of
  systems with their AMI and Lyapunov responses.
- Systems: Lorenz, Ikeda, Chen, bimodal switching.

#### Changed
- Don't save transients; `mkpath` when saving timeseries.
- Don't `eval` git commit IDs when loading processes.

#### Fixed
- Phase-corruption simulation bugs (wrong times passed; transient handling);
  threshold corruption on power.

### 2021-05 --- Chaotic attractors and trajectory plotting

#### Added
- Systems: diffusionless Lorenz, 4D simplified Lorenz, double scroll, piecewise
  linear hyperchaotic, Thomas' cyclically symmetric attractor, custom
  high-dimensional system.
- High-quality 3D trajectory plotting; colour by 4th variable and by phase-space
  velocity; trajectories traced along axes.
- `rampOn`/`rampOff` and `rampInterval` parameter profiles; parameter annotation.
- File IO for timeseries.

#### Fixed
- Chaotic-flow EOM parameter bugs; high-accuracy simulation defaults.

### 2021-03 to 2021-04 --- DifferentialEquations rewrite

#### Changed
- **Ditched DynamicalSystems.jl** for streamlined simulations built directly on
  DifferentialEquations.jl; default solver RK4.
- Plot recipe for `Process`es with cleaner syntax.
- Field aliases for the `Process` constructor; `eval` String-valued fields.

#### Added
- `Discontinuous` type and `unitStep`; discontinuities propagate through
  `tuplef2ftuple` to solvers for multidimensional processes.
- Parameter profiles: ramp, triangle wave, Lorentzian, composite stepped.
- Systems: Henon, ARMA, simplest chaotic flow, bimodal switching, Gaussian
  bimodal, driven oscillators, scaly/dynamical sines, double pendulum, fmWave.
- Save-able simulation methods; time-series file IO.

### 2021-02 --- Initial

#### Added
- Initial package from template; DynamicalSystems.jl examples.
- Time-varying and discontinuous parameter investigations (e.g. Van der Pol).
