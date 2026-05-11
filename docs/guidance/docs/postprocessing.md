# Postprocessing and Output Data

CHAPSim2 produces several classes of output. The postprocessing approach depends on whether restart data, flow-field visualization, monitor histories, or statistical profiles are required.

## Output Folders

| Folder | Main content | Typical use |
| --- | --- | --- |
| `1_data/` | Restart/checkpoint data and domain files | Restart initialization, mesh interpolation, wall-unit post-processing, and low-level field validation |
| `2_visu/` | Visualization files (typically XDMF/HDF5 format) | ParaView or VisIt inspection, field slices, and instantaneous flow snapshots |
| `3_monitor/` | Monitor history files and plotting scripts | Bulk quantities, point probes, convergence verification, and transient diagnostics |
| `4_check/` | Mesh and setup validation outputs | Grid quality assessment, spatial distribution verification, and early validation |

Not all simulations produce every folder. Output availability depends on the switches and frequencies defined in `[io]`, `[statistics]`, and related input-file sections.

## Visualisation Data

Visualization is controlled primarily by `visu_idim`, `visu_nfre`, and `visu_nskip` in the input file. Full-domain and plane outputs can be visualized in ParaView or VisIt when XDMF format files are present.

Common example scripts:

| Script | Purpose |
| --- | --- |
| `examples/channel_iso_periodic/case/2_visu/plot_channel_velo_stress.py` | Plot channel mean velocity and Reynolds-stress profiles. |
| `examples/channel_iso_periodic/case/2_visu/postprocess_channel_wall_units.py` | Convert channel output to wall-unit profiles. |
| `examples/pipe_iso_periodic/case/2_visu/plot_pipe_velo_stress_v2.py` | Plot pipe velocity and stress statistics with robust figure saving. |
| `examples/pipe_iso_periodic/case/2_visu/plot_pipe_velo_stress.py` | Earlier pipe profile plotting script. |

For production postprocessing, copy the closest example script into the case
directory, update the run index, mesh/domain assumptions, and file paths, then
run it from the folder containing the output files.

## Monitor Data

Monitor output is useful for assessing run health before expensive statistics accumulation is performed. Reference examples include:

| Script | Purpose |
| --- | --- |
| `examples/3_monitor/plot_monitor_bulk_change_history.py` | Plot bulk-value evolution. |
| `examples/3_monitor/plot_monitor_points.py` | Plot monitor-point histories. |
| `examples/3_monitor/plot_monitor_points_seperate.py` | Plot monitor points in separate figures. |
| `examples/channel_iso_periodic/case/3_monitor/plot_monitor_bulk_change_history.py` | Case-local channel monitor plotting. |
| `examples/pipe_iso_periodic/case/3_monitor/plot_monitor_points.py` | Case-local pipe monitor plotting. |

Recommended checks:

- CFL remains within the intended range.
- Flow rate, pressure gradient, or driving force reaches the expected regime.
- Monitor points do not show sudden jumps after restart.
- Statistics are not started until the transient period has passed.

## Statistics Levels

The input variable `stat_level` controls the extent of statistical information collected. Select the minimum level that contains required quantities, as higher levels increase memory, I/O, and post-processing costs.

| `stat_level` | Content | Typical use |
| --- | --- | --- |
| `0` | Statistics disabled | Initial setup verification, mesh validation, and short test runs |
| `1` | Basic mean quantities | Mean-flow development and preliminary monitoring |
| `2` | Second-order statistics | Reynolds stresses and wall-bounded turbulence statistics |
| `3` | Extended statistics (where supported) | Detailed budgets and advanced diagnostics |

Control sampling behavior with `stat_istart`, `stat_nskip`, and related parameters to ensure statistics are collected only after the transient initial condition phase.

## Mesh and Setup Checks

Before trusting a new case, inspect the mesh:

```bash
python3 examples/4_check/plot_check_mesh.py
```

Case-local copies also exist, for example:

- `examples/channel_iso_periodic/case/4_check/plot_check_mesh.py`
- `examples/pipe_iso_periodic/case/4_check/plot_check_mesh.py`

Check near-wall spacing, stretching smoothness, and whether the physical domain
matches the intended geometry.

## Recommended Workflow

1. Execute a short test case and review log diagnostics
2. Plot monitor histories from `3_monitor/`
3. Inspect mesh and initial condition outputs from `4_check/`
4. Open XDMF visualization outputs from `2_visu/`
5. Enable statistics collection only after flow development is complete
6. Use the closest example plotting script as a template for profile analysis

**Note on mounted filesystems:** Figure saving can occasionally fail with I/O errors on remote or mounted file systems. A robust approach is to save figures to a local temporary directory first, verify completion, and copy to the case directory when the file system is stable.

For interactive tuning of the y-direction stretching before a run, see the
[Mesh Stretching Reviewer](mesh-reviewer.md).
