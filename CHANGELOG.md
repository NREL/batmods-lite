# BATMODS Changelog

## [Unreleased](https://github.com/NREL/BATMODS-lite/)

### New Features
- `Exerpiment` class interface allows for multi-step and dynamic load profiles
- Simplified `run` and `run_step` methods to work with generalized `Experiment` class
- New solver now provides additional statistics, e.g., number of function and Jacobian evaluations
- Control-specific `Solution` classes are replaced with generalized `StepSolution` and `CycleSolution`

### Optimizations
- Replaced `scikits-odes` solver with `scikit-sundae` to improve installation
- Sped up SPM by flipping cathode pointers so that bandwidth is reduced to +/- 2
- Use a generic `bandwidth` function instead of having one per model subpackage

### Bug Fixes

### Breaking Changes
- Removal of `run_CC`, `run_CV`, and `run_CP` methods
- Renamed some attributes so they are no longer user-facing (`_sv0`, `_algidx`, etc.)

## [v0.0.1](https://github.com/NREL/BATMODS-lite/tree/v0.0.1)
This is the first release of BATMODS-lite. Main features/capabilities are listed below.

### Features
- Working SPM and P2D models
- Capability to run with current, voltage, or power demands

### Notes
- Models were verified mathematically by checking conservation equations
- Validations were completed using COMSOL equivalents
