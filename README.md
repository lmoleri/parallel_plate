# Nanoscale Parallel-Plate Garfield++ Scan

This project studies a two-plate gas gap in Garfield++ using `ComponentAnalyticField`, `MediumMagboltz`, and `AvalancheMicroscopic`. A single seed electron is injected near the source plate, avalanche multiplication is enabled, and the analysis records only electrons that hit the arrival plate.

Runtime configs live in the dedicated `config/` folder as JSON files, following the same study-style pattern used in `SEE_in_vacuum`.

The baseline gas mixture is Ar:CO2 (93:7) with Penning transfer enabled. The default starter scans are:

- `d` scan at fixed `p = 760 Torr`, `E = 200 kV/cm`: `200, 250, 300, 350 nm`
- `p` scan at fixed `d = 300 nm`, `E = 200 kV/cm`: `200, 400, 600, 760 Torr`
- `E` scan at fixed `d = 300 nm`, `p = 760 Torr`: `50, 100, 200, 400 kV/cm`

## Dependencies

This repository does not commit Garfield++ itself. You can use either:

- a normal Garfield++ installation exposed via `GARFIELD_INSTALL` or `CMAKE_PREFIX_PATH`
- the shared workspace installation at `../../local/garfield`

If you are working in the same multi-project workspace, the shared install can be activated with:

```bash
source ../../local/garfield/share/Garfield/setupGarfield.sh
```

## Optional Shared-Workspace Bootstrap

If this repo lives inside the same workspace as a shared Garfield++ source checkout under `../../vendor/garfieldpp`, you can rebuild the shared install with:

```bash
./scripts/bootstrap_garfield.sh
```

The bootstrap script does not vendor Garfield++ into this repository; it targets the workspace-level shared install.

## Workspace Boundary

This project owns everything under `projects/parallel_plate`.

The shared Garfield++ dependency assets remain at the workspace root so future Garfield++ projects can reuse them:

- `../../vendor/garfieldpp`
- `../../vendor/garfieldpp-build`
- `../../local/garfield`

## Build

```bash
source ../../local/garfield/share/Garfield/setupGarfield.sh
cmake -S . -B build
cmake --build build -j4
```

If you are not using the shared workspace install, set `GARFIELD_INSTALL` or `CMAKE_PREFIX_PATH` before configuring.

## Run

Run all three scans with the default configuration:

```bash
./build/parallel_plate_scan \
  --config config/default_parallel_plate.json \
  --scan all \
  --out results
```

Run a single scan:

```bash
./build/parallel_plate_scan --scan d --out results
./build/parallel_plate_scan --scan p --out results
./build/parallel_plate_scan --scan field --out results
```

## Outputs

The `--out` argument is the parent directory for results. The executable always creates a deterministic parameter-based run subfolder inside it, for example:

- `results/scan-field__d-300nm__p-760Torr__E-50-100-200-400kVcm__n-500/`
- `results/scan-all__d-200-250-300-350nm__p-200-400-600-760Torr__E-50-100-200-400kVcm__n-500/`

Each run folder writes:

- `parallel_plate_scan.root`: ROOT file with per-point histograms and per-scan `TGraphErrors`
- `scan_summary.csv`: flat summary table with scan point statistics
- `run_config_used.json`: effective configuration used for the run after defaults are applied
- `run_manifest.json`: resolved scan/output metadata and the original command line
- `summary/*.png`: two-panel summary plots for mean arrival energy and mean arrival multiplicity
- `histograms/<scan>/*.png`: two-panel per-point plots of arrival-energy and multiplicity histograms

## Config Schema

The shipped configs are:

- `config/default_parallel_plate.json`
- `config/smoke_parallel_plate.json`

They use this JSON schema:

```json
{
  "baseline": {
    "distance_nm": 300.0,
    "pressure_torr": 760.0,
    "field_kv_per_cm": 200.0
  },
  "scan": {
    "distance_nm": [200.0, 250.0, 300.0, 350.0],
    "pressure_torr": [200.0, 400.0, 600.0, 760.0],
    "field_kv_per_cm": [50.0, 100.0, 200.0, 400.0]
  },
  "gas": {
    "enable_penning": true,
    "enable_thermal_motion": true
  },
  "simulation": {
    "num_primaries": 500,
    "initial_electron_energy_ev": 0.1,
    "max_electron_energy_ev": 200.0,
    "avalanche_size_limit": 100000,
    "temperature_k": 293.15
  },
  "histogram": {
    "energy_bins": 200,
    "multiplicity_bins_min": 20
  }
}
```

The executable classifies arrivals geometrically at the arrival boundary. For `ComponentAnalyticField`, Garfield can report a plate crossing either as `StatusHitPlane` or as `StatusLeftDriftMedium` with an endpoint on the plate boundary, so the analysis accepts both cases when the endpoint is at `x = d` within a small tolerance. Mean arrival energy is therefore defined over successfully arriving electrons only, while multiplicity is defined per primary and includes zero-arrival events.
