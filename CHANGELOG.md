# Changelog

All notable changes to the ABACUS INPUT extension will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-02-24

- **relax_method**: Fixed validation for two-element vector values (e.g., "cg 1", "bfgs 2"). 
- **Two-element vector parameters**: Implemented **fully generic validation logic** that auto-detects two-element vectors by `mdType` pattern and reads constraints from JSON:
  - Automatically recognizes patterns: `"Boolean [Integer](optional)"`, `"Integer [Integer](optional)"`, `"Vector of string"`
  - Second element validation configured in `constraints.json` via `secondElement` field
  - `relax_method`: Second element enum `["1", "2"]` defined in JSON (no code changes needed)
  - `out_mat_l`: Second element type `integer` inferred from `mdType`
  - `out_xc_r`: Second element type `integer` inferred from `mdType`
  - **Future-proof**: New two-element parameters only need proper `mdType` and optional `secondElement` in constraints.json, **zero TypeScript changes required**
- **bug fixes**
## [0.1.0] - 2026-02-22

- Initial development version
