# Benchmark and Example Cases

The recommended approach for creating a new CHAPSim2 configuration is to locate the most similar existing case and duplicate it, then edit `input_chapsim.ini` according to your requirements. The `tests/` directory contains compact validation cases designed for smoke and regression testing. The `examples/` directory contains postprocessing scripts and reference data for selected canonical production-scale cases.

## Case Families

| Family | Representative test cases | Primary use |
| --- | --- | --- |
| Taylor-Green vortex | `tgv_iso`, `tgv_scp` | Periodic validation for core numerics, scalar/thermal coupling, and regression metrics |
| Channel flow | `channel_iso_periodic`, `channel_iso_inout`, `channel_scp_*` | Wall-bounded Cartesian configurations (periodic or inlet/outlet, isothermal or thermal) |
| Pipe flow | `pipe_iso_periodic`, `pipe_iso_inout`, `pipe_scp_*` | Cylindrical wall-bounded cases with radial/azimuthal constraints |
| Annular flow | `annular_iso_periodic`, `annular_iso_inout`, `annular_scp_*` | Cylindrical annular geometries with inner and outer wall treatment |

Suffixes in case names are used consistently:

| Suffix | Definition |
| --- | --- |
| `iso` | Isothermal flow configuration |
| `scp` | Scalar/thermal property case |
| `periodic` | Periodic streamwise direction (typically pressure-gradient or flow-rate controlled) |
| `inout` | Inlet/outlet boundary condition |
| `Tw` | Wall-temperature thermal boundary condition |
| `qw` | Wall-heat-flux thermal boundary condition |

## Use Cases and Starting Configurations

| Objective | Recommended starting case |
| --- | --- |
| Initial solver validation | `tests/tgv_iso` |
| Periodic channel Direct Numerical Simulation | `tests/channel_iso_periodic` |
| Channel with inlet/outlet conditions | `tests/channel_iso_inout` |
| Periodic pipe Direct Numerical Simulation | `tests/pipe_iso_periodic` |
| Pipe with inlet/outlet conditions | `tests/pipe_iso_inout` |
| Annular periodic configuration | `tests/annular_iso_periodic` |
| Thermal case with wall-temperature | Closest `*_Tw` case |
| Thermal case with wall-heat-flux | Closest `*_qw` case |

After copying a case, verify the following parameters:

1. `[domain] icase` and domain spatial extent
2. `[mesh] ncx`, `ncy`, `ncz`, `istret`, and `rstret`
3. `[flow] initfl`, `ren`, `idriven`, and boundary-condition selections
4. `[thermo]` only when thermal or scalar physics is enabled
5. `[io]` and `[statistics]` output frequencies before executing production simulations

## Example Postprocessing Assets

| Example path | Content |
| --- | --- |
| `examples/channel_iso_periodic/case/2_visu/` | Channel velocity/stress plotting and wall-unit postprocessing scripts. |
| `examples/channel_iso_periodic/MKM180_profiles/` | Reference profile data for channel comparison. |
| `examples/channel_iso_periodic/MKM395_profiles/` | Higher-Re channel reference profiles. |
| `examples/pipe_iso_periodic/case/2_visu/` | Pipe velocity/stress plotting scripts. |
| `examples/pipe_iso_periodic/TDL180/` | Pipe reference data at friction Reynolds number near 180. |
| `examples/pipe_iso_periodic/TDL550/` | Pipe reference data at higher friction Reynolds number. |
| `examples/3_monitor/` | Generic monitor plotting scripts. |
| `examples/4_check/` | Generic mesh-check plotting script. |

## How Benchmarks Fit the Workflow

Use benchmark cases in three sequential stages:

1. **Build confidence**: Run a small smoke test or Taylor-Green vortex case immediately after compilation
2. **Create new configuration**: Copy the most similar case and adapt the input file to your requirements
3. **Validate results**: Compare monitor histories, mean profiles, stresses, and available reference data

For automated checking, see [Regression and Smoke Tests](testing.md). For
profile plotting and visualisation, see [Postprocessing and Output Data](postprocessing.md).
