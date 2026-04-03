# Changelog

All notable changes to AtomTwin are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

---

## [Pre-release]

### Added
- Unit test suite covering `Parameter`, `ParametricExpression`, `Level`/manifolds, `Sequence`, and `System` construction.
- GitHub Actions CI workflow (`.github/workflows/CI.yml`) running on Julia 1.11 and nightly.
- `LICENSE` file (Apache-2.0).
- Narrative documentation pages for Parameters & Noise and Visualization APIs.
- Vulnerability reporting instructions added to `CONTRIBUTING.md`.

### Changed
- `Revise` removed from package dependencies (was incorrectly listed as a runtime dependency; it is a developer tool).
- `Test` moved from `[deps]` to `[extras]`/`[targets]` per Julia packaging convention.
- `julia = "1.11"` compat bound added to `Project.toml`.
- Internal DAG node types (`BeamNode`, `CouplingNode`, etc.) removed from public exports; only `AbstractNode` is exported.
- README placeholder URLs replaced with real repository and documentation links.

---

## [0.1.0] — 2026

Initial beta release.

### Features

- **Atom / level model**: `Level`, `HyperfineLevel`, `FineLevel`, `HyperfineManifold`, `FineManifold`, `Superposition`.
- **Atomic species**: `Ytterbium171Atom`, `Potassium39Atom`, `Rubidium87Atom`, `Strontium88Atom`.
- **Gaussian beam model**: `GaussianBeam`, `GeneralGaussianBeam`, `PlanarBeam`.
- **System building**: `System`, `add_coupling!`, `add_detuning!`, `add_decay!`, `add_dephasing!`, `add_interaction!`, `add_zeeman_detunings!`.
- **Parametric simulations**: `Parameter`, `ParametricExpression`, `MixedPolarization`.
- **Laser phase noise**: `LaserPhaseNoiseModel`, `NoisyField`, `laser_freq_psd`, `laser_phase_psd`.
- **Tweezer arrays**: `TweezerArray`.
- **Instruction sequences**: `Sequence`, `@sequence`, `Wait`, `Pulse`, `On`, `Off`, `MoveRow`, `MoveCol`, `RampRow`, `RampCol`, `AmplRow`, `AmplCol`, `FreqRow`, `FreqCol`, `Parallel`.
- **Simulation engine**: `play`, `compile`, `recompile!`.
- **Detectors**: `PopulationDetectorSpec`, `CoherenceDetectorSpec`, `MotionDetectorSpec`, `FieldDetectorSpec`.
- **Analysis**: `process_tomography`, `estimate_PER`, `estimate_PER_dB`.
- **Visualization**: `AtomTwin.Visualization.animate` (requires GLMakie extension).
- **Plotting recipes** via RecipesBase (Plots.jl extension).
- 11 worked examples covering Rabi oscillations, dissipation, motion, laser noise, Rydberg blockade, time-optimal gates, Ramsey/EIT, Yb-171 Raman gates, K-39 state preparation, process tomography, and atom sorting.

### Known limitations

- Visualization is currently limited to 2D atomic trajectory animation
- The General Registry registration is pending.
