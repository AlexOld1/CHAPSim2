# CHAPSim Input File Guide

CHAPSim2 reads simulation parameters from an INI-style configuration file, typically named `input_chapsim.ini`. Input generation tools may produce files named `input_chapsim_auto.ini` or `input_chapsim_gui.ini`; before solver execution, ensure the selected file is placed or renamed as `input_chapsim.ini` in the case directory.

The file parser operates on a section-based structure; however, variable order within sections is significant. Maintain variable order as shown in generated templates.

## What This File Controls

The input file establishes the primary interface between user specifications and the Fortran solver. It controls:

- **Physical configuration**: Geometry specification, Reynolds number, thermal/magnetohydrodynamic options, and working-fluid properties
- **Numerical discretization**: Grid resolution, domain extent, mesh stretching, time-stepping, and spatial/temporal discretization schemes
- **Boundary conditions**: Periodic directions, inlet/outlet treatment, wall velocity specifications, wall temperature, and wall heat flux
- **Simulation control**: Iteration ranges for flow and thermal field computation
- **Output specification**: Restart checkpoint frequency, visualization output frequency, statistics accumulation parameters, and database plane I/O
- **Diagnostics**: Probe location definitions and monitoring-output frequency

For new configurations, begin with the most similar existing `tests/*/input_chapsim.ini` template or generate a file using `prepost/autoinput/autoinput_script.py` or `autoinput_gui.py`, then edit accordingly.

## Basic Rules

| Property | Format |
| --- | --- |
| Boolean values | Fortran logical syntax: `.true.` or `.false.` |
| Integer values | Whole-number identifiers or counters |
| Real values | Decimal or scientific notation (e.g., `1e-05`) |
| List values | Comma-separated (e.g., `veloinit= 0.0,0.0,0.0`) |
| Domain lengths | Nondimensional with respect to reference half-height, radius, or equivalent case length |
| Thermal parameters | SI units before solver nondimensionalization |
| Reynolds numbers | Based on channel half-height, pipe radius, or equivalent case reference length |
| Boundary condition rows | Format: `bc_low,bc_high,value_low,value_high` |
| Comments | Lines beginning with `#` or `;`, blank lines, or indented comment lines are ignored |

## `[process]`

Governs high-level execution mode selection.

| Variable | Fortran type | Meaning |
|---|---|---|
| `is_prerun` | `logical` | If `.true.`, run preprocessing/recommendation logic only. |
| `is_postprocess` | `logical` | If `.true.`, run postprocessing mode. |

## `[decomposition]`

Controls MPI/domain decomposition.

| Variable | Fortran type | Meaning |
|---|---|---|
| `nxdomain` | `integer` | Number of domains in `x`. Current production inputs should use `1`. |
| `p_row` | `integer` | MPI process grid rows, usually aligned with `y`. Use `0` for automatic decomposition. |
| `p_col` | `integer` | MPI process grid columns, usually aligned with `z`. Use `0` for automatic decomposition. |

## `[domain]`

Defines the physical case and domain extents.

| Variable | Fortran type | Meaning |
|---|---|---|
| `icase` | `integer` | Flow geometry/case ID. |
| `lxx` | `real` | Domain length in `x`. For channel, pipe, and annular cases this is usually streamwise length. |
| `lyt` | `real` | Upper `y`/radial boundary. Some cases reset this internally. |
| `lyb` | `real` | Lower `y`/radial boundary. Some cases reset this internally. |
| `lzz` | `real` | Domain length in `z`. For pipe/annular this is azimuthal length and is reset to `2π`. |

Case IDs:

| ID | Case |
|---:|---|
| 1 | Channel |
| 2 | Pipe |
| 3 | Annular |
| 4 | 3-D Taylor-Green vortex |
| 5 | Duct |

Notes:

- Pipe and annular cases use cylindrical coordinates internally.
- For pipe, `lyb`, `lyt`, and `lzz` are reset to `0`, `1`, and `2π`.
- For annular flow, `lyt` and `lzz` are reset to `1` and `2π`.
- For duct, `x` and `y` are wall-normal directions and `z` is streamwise.

## `[flow]`

Defines flow-field initialisation and Reynolds numbers.

| Variable | Fortran type | Meaning |
|---|---|---|
| `initfl` | `integer` | Flow-field initialisation method ID. |
| `irestartfrom` | `integer` | Restart iteration for flow when `initfl=0`; otherwise reset internally to `0`. |
| `veloinit` | `real(3)` | Constant initial velocity vector used when `initfl=4`. |
| `noiselevel` | `real` | Random perturbation amplitude added during initialisation where applicable. |
| `reni` | `real` | Initial Reynolds number used for ramping/scaling. |
| `nreni` | `integer` | Number of iterations over which the initial Reynolds setting is applied. |
| `ren` | `real` | Target Reynolds number for the run. |

Initialisation IDs:

| ID | Meaning |
|---:|---|
| 0 | Restart from saved fields |
| 2 | Random perturbation |
| 3 | Initialise from inlet data |
| 4 | Given constant values |
| 5 | Poiseuille profile |
| 6 | Analytic function, used for Taylor-Green vortex |
| 7 | Given/mixed boundary-condition initialisation |

Common choices are `initfl=5` for periodic channel, pipe, and annular cases;
`initfl=3` for inlet/outlet cases; and `initfl=6` for Taylor-Green vortex.

## `[thermo]`

Defines thermal/energy-equation settings. Include this section when solving the
energy equation. Some templates include placeholder thermal values even for
isothermal cases because the input format is shared.

| Variable | Fortran type | Meaning |
|---|---|---|
| `ithermo` | `logical` | Enables thermal/energy equation. |
| `icht` | `logical` | Enables conjugate heat transfer mode. |
| `igravity` | `integer` | Gravity direction ID. |
| `ifluid` | `integer` | Working-fluid property model ID. |
| `ref_l0` | `real` | Dimensional reference length in metres. |
| `ref_t0` | `real` | Reference temperature in Kelvin. |
| `inittm` | `integer` | Thermal-field initialisation method ID. Uses the same IDs as `initfl`. |
| `irestartfrom` | `integer` | Restart iteration for thermal field when `inittm=0`. |
| `tini` | `real` | Initial temperature in Kelvin. |
| `inout_buffer` | `real(2)` | Inlet and outlet thermal buffer lengths as `inlet,outlet`, scaled by `L0`. |
| `qw_ramp` | `logical,integer,integer` | Heat-flux ramp as `enabled,start_iter,end_iter`. |

Gravity IDs:

| ID | Direction |
|---:|---|
| 0 | No gravity |
| 1 | +x |
| -1 | -x |
| 2 | +y |
| -2 | -y |
| 3 | +z |
| -3 | -z |

Fluid IDs:

| ID | Fluid |
|---:|---|
| 1 | Supercritical water |
| 2 | Supercritical CO2 |
| 3 | Liquid sodium |
| 4 | Liquid lead |
| 5 | Liquid bismuth |
| 6 | Liquid LBE |
| 7 | Liquid water |
| 8 | Liquid lithium |
| 9 | Liquid FLiBe |
| 10 | Liquid PbLi eutectic |

## `[mhd]`

Defines magnetohydrodynamics settings. If MHD is disabled, the section may still
contain placeholder values.

| Variable | Fortran type | Meaning |
|---|---|---|
| `imhd` | `logical` | Enables MHD model. |
| `NStuart` | `logical,real` | Pair `enabled,value` for Stuart number. |
| `NHartmn` | `logical,real` | Pair `enabled,value` for Hartmann number. |
| `B_static` | `real(3)` | Static magnetic-field vector `Bx,By,Bz`. |

Exactly one of `NStuart` or `NHartmn` should be enabled for an MHD run.

## `[mesh]`

Defines grid resolution and wall-normal/radial stretching.
Use the [Mesh Stretching Reviewer](mesh-reviewer.md) to inspect the y-direction
mapping before running an expensive case.

| Variable | Fortran type | Meaning |
|---|---|---|
| `ncx` | `integer` | Number of cells in `x`. |
| `ncy` | `integer` | Number of cells in `y` or radial direction. |
| `ncz` | `integer` | Number of cells in `z` or azimuthal direction. For cylindrical cases, odd values are increased to the next even value. |
| `istret` | `integer` | Mesh stretching type ID. |
| `rstret` | `integer,real` | Stretching method and factor as `method,factor`. |

Stretching type IDs for `istret`:

| ID | Meaning |
|---:|---|
| 0 | No stretching |
| 1 | Centre clustering |
| 2 | Two-side clustering |
| 3 | Bottom-side clustering |
| 4 | Top-side clustering |

Stretching method IDs for the first value of `rstret`:

| ID | Meaning |
|---:|---|
| 1 | Five-mode spectral stretching |
| 2 | Tanh stretching method |
| 3 | Power-law stretching method |

Five-mode spectral stretching uses a smooth analytic stretching function whose
spectral representation contains only five modes, `-2`, `-1`, `0`, `1`, and `2`.
This compact support is useful when the stretching effect is handled as a
convolution in the FFT spectral domain.

Recommended defaults are two-side clustering for channel and annular flow,
top-side clustering for pipe flow, and no stretching for Taylor-Green vortex.

## `[bc]`

Defines boundary conditions and periodic-flow driving.

Each boundary-condition line has four values:

```ini
ifbcx_u= bc_xlow,bc_xhigh,value_xlow,value_xhigh
```

The prefix gives the direction (`ifbcx`, `ifbcy`, `ifbcz`) and the suffix gives
the variable (`u`, `v`, `w`, `p`, `t`). The first BC/value pair belongs to the
lower/start boundary, and the second belongs to the upper/end boundary.

| Variable | Fortran type | Meaning |
|---|---|---|
| `ifbcx_u`, `ifbcx_v`, `ifbcx_w` | `integer,integer,real,real` | Velocity BCs on x boundaries. |
| `ifbcx_p` | `integer,integer,real,real` | Pressure BC on x boundaries. |
| `ifbcx_t` | `integer,integer,real,real` | Temperature BC on x boundaries. Temperature values are Kelvin for Dirichlet; heat flux values are W/m² before nondimensionalisation for Neumann. |
| `ifbcy_u`, `ifbcy_v`, `ifbcy_w` | `integer,integer,real,real` | Velocity BCs on y/radial boundaries. |
| `ifbcy_p` | `integer,integer,real,real` | Pressure BC on y/radial boundaries. |
| `ifbcy_t` | `integer,integer,real,real` | Temperature BC on y/radial boundaries. |
| `ifbcz_u`, `ifbcz_v`, `ifbcz_w` | `integer,integer,real,real` | Velocity BCs on z/azimuthal boundaries. |
| `ifbcz_p` | `integer,integer,real,real` | Pressure BC on z/azimuthal boundaries. |
| `ifbcz_t` | `integer,integer,real,real` | Temperature BC on z/azimuthal boundaries. |
| `idriven` | `integer` | Flow-driving method ID. |
| `drivenfc` | `real` | Magnitude for wall-shear or pressure-gradient driving. Mass-flux driving normally uses `0.0`. |

Boundary condition IDs:

| ID | Meaning | Typical use |
|---:|---|---|
| 0 | Interior | Pipe axis or internal boundary |
| 1 | Periodic | Periodic directions |
| 2 | Symmetric | Symmetry plane |
| 3 | Antisymmetric | Antisymmetry plane |
| 4 | Dirichlet | Fixed value |
| 5 | Neumann | Fixed gradient or heat flux |
| 6 | Interpolation | Internal/interpolation use |
| 7 | Convective outlet | Open outlet, only supported in `x` or `z` |
| 9 | Profile inlet | 1-D/profile inlet, not supported on all faces |
| 10 | Database inlet | Inlet from stored plane data |
| 11 | Poiseuille | Nominal Poiseuille BC |
| 12 | Other/interpolation | Special cases |

Flow driving IDs for `idriven`:

| ID | Meaning |
|---:|---|
| 0 | No forcing |
| 1 | Constant streamwise mass flux in `x` |
| 2 | Constant wall shear in `x` |
| 3 | Constant pressure gradient in `x` |
| 4 | Constant streamwise mass flux in `z` |
| 5 | Constant wall shear in `z` |
| 6 | Constant pressure gradient in `z` |

Important constraints:

- Boundary rows are given for five variables in this order: `u`, `v`, `w`, `p`,
  and `T`. Even isothermal runs still include the `T` rows.
- Use driving only for periodic wall-bounded cases. In open inlet/outlet cases,
  the solver disables flow driving.
- Convective outlet in `y` is not supported.
- For pipe, the lower radial boundary is treated internally as the axis/interior.
- If any side of a variable is periodic, both sides for that variable are made periodic.
- For database inlet in `x`, the solver applies database treatment to all velocity
  components and Neumann treatment to pressure.

## `[scheme]`

Defines time integration and spatial discretisation.

| Variable | Fortran type | Meaning |
|---|---|---|
| `dt` | `real` | Time-step size. |
| `itimescheme` | `integer` | Time integration method ID. |
| `iaccuracy` | `integer` | Spatial derivative accuracy ID. |
| `iviscous` | `integer` | Viscous-term treatment ID. |
| `out_sponge_l_re` | `real(2)` | Outlet sponge layer as `length,Re_strength`. Use nonzero length for open outlet cases if needed. |

Time scheme IDs for `itimescheme`:

| ID | Meaning |
|---:|---|
| 0 | Euler |
| 1 | Adams-Bashforth 2 |
| 2 | RK3-Crank-Nicolson |
| 3 | RK3 |

Spatial accuracy IDs for `iaccuracy`:

| ID | Meaning |
|---:|---|
| 1 | 2nd-order central difference |
| 2 | 4th-order central difference |
| 3 | 4th-order compact |
| 4 | 6th-order compact |

Viscous treatment IDs for `iviscous`:

| ID | Meaning |
|---:|---|
| 1 | Explicit viscous treatment |
| 2 | Semi-implicit viscous treatment |

For cylindrical coordinates, the solver forces `iaccuracy` to 2nd-order central
difference. For channel cases, compact schemes are reduced to 4th-order central
difference.

## `[simcontrol]`

Defines iteration ranges.

| Variable | Fortran type | Meaning |
|---|---|---|
| `niterflowfirst` | `integer` | First flow iteration to run. |
| `niterflowlast` | `integer` | Last flow iteration to run. |
| `niterthermofirst` | `integer` | First thermal iteration; use `0` when thermal is disabled. |
| `niterthermolast` | `integer` | Last thermal iteration; use `0` when thermal is disabled. |

The environment variable `CHAPSIM_NITER` can override the final iteration for
short smoke tests.

## `[io]`

Controls monitor, restart, visualisation, statistics, and plane-database I/O.

| Variable | Fortran type | Meaning |
|---|---|---|
| `cpu_nfre` | `integer` | Frequency for CPU/progress output. |
| `ckpt_nfre` | `integer` | Checkpoint/restart write frequency. |
| `visu_idim` | `integer` | Visualisation mode ID. |
| `visu_nfre` | `integer` | Visualisation output frequency. |
| `visu_nskip` | `integer(3)` | Cell skip for visualisation output in `x,y,z`. |
| `stat_istart` | `integer` | Iteration at which statistics begin. |
| `stat_level` | `integer` | Statistics level. |
| `stat_nskip` | `integer(3)` | Cell skip for statistics in `x,y,z`. |
| `is_wrt_read_bc` | `logical,logical` | Pair `write_outlet,read_inlet` for plane database files. |
| `wrt_read_nfre` | `integer(3)` | Plane database frequency and range as `frequency,start,end`. |
| `io_mode` | `integer` | Existing-file handling mode ID. |

Visualisation mode IDs for `visu_idim`:

| ID | Meaning |
|---:|---|
| 0 | 3-D only |
| 1 | 2-D planes only |
| 2 | Both 3-D and 2-D outputs |

Statistics levels for `stat_level`:

| ID | Meaning |
|---:|---|
| 0 | No statistics |
| 1 | Mean/first moments |
| 2 | Second moments |
| 3 | Extended/turbulent budget statistics where supported |

I/O mode IDs for `io_mode`:

| ID | Meaning |
|---:|---|
| 0 | Overwrite existing files |
| 1 | Skip write if file exists |
| 2 | Rename existing file before writing |

For database inlet cases, set `is_wrt_read_bc= .false.,.true.` and make sure
the requested inlet database files exist.

## `[probe]`

Defines point probes for monitor output.

| Variable | Fortran type | Meaning |
|---|---|---|
| `npp` | `integer` | Number of probe points. |
| `pt1`, `pt2`, ... | `real(3)` | Probe coordinates as `x,y,z`. |

Probe coordinates use the same nondimensional coordinate system as the domain.

## Practical Advice

### Recommended Setup Workflow

1. Select the closest existing case under `tests/`.
2. Generate or copy an input file.
3. Check `[domain]`, `[mesh]`, `[flow]`, and `[bc]` first; these define the core
   physics and numerics.
4. Enable `[thermo]` or `[mhd]` only when the case needs those physics.
5. Keep the first run short by reducing `niterflowlast`, or by using the
   `CHAPSIM_NITER` environment override.
6. Increase output frequencies only after the run is stable.

### Case Setup Checks

- Keep `nxdomain= 1` unless the solver is extended to support multiple x-domains.
- Use `p_row= 0` and `p_col= 0` unless a specific MPI decomposition is required.
- For periodic channel, pipe, and annular cases, use periodic streamwise BCs and
  a driving method such as `idriven=1`.
- For open inlet/outlet cases, use `initfl=3`, `idriven=0`, and database reading
  when `ifbcx_u` or `ifbcz_u` uses BC `10`.
- For thermal wall-heating cases, confirm the sign convention of Neumann heat
  flux before launching a long production run.
- For pipe and annular cases, remember that cylindrical-coordinate constraints
  can override some mesh and accuracy choices.

### Runtime Monitoring

During early runs, watch:

- Input-reading messages for unexpected case, mesh, BC, or scheme overrides.
- CFL and time-step diagnostics.
- Mass-conservation checks.
- Wall quantities and driven-flow response.
- Restart, visualisation, statistics, and plane-database output frequencies.

Treat warnings during input reading as setup feedback. They often indicate that
the solver has corrected an unsupported or inconsistent option.
