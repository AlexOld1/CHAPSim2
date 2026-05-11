# Changelog

## Unreleased - 2026-04-27

### Added

- Added a unified pipe-axis halo and centre reconstruction routine,
  `axis_mirror_fbcy`, replacing the older even/odd-only mirroring routines.
- Added axis reconstruction modes for regular cylindrical centreline behaviour:
  `AXIS_RECON_NONE`, `AXIS_RECON_ZERO`, `AXIS_RECON_M0`, `AXIS_RECON_M1`, and
  `AXIS_RECON_M0_M2`.
- Added pipe-centre reconstruction support in momentum, energy, statistics,
  pressure-gradient, and MHD-related halo updates so scalar-like, vector-like,
  and quadratic/tensor-like terms use appropriate centreline regularity.
- Added a random-initialisation envelope that damps perturbations near channel,
  pipe, and annular boundaries instead of applying a flat perturbation level.
- Added an interactive `build/Makefile` default target that can optionally clean
  first and choose between default, GNU, Intel, Cray, and NVHPC CPU build modes.

### Changed

- Reworked pipe-centre treatment for cylindrical flows to enforce single-valued
  scalar quantities and regular first/second azimuthal-mode behaviour at the
  axis.
- Updated thermal-property interpolation and energy RHS assembly to pass
  boundary-condition halos into midpoint and derivative operations, then
  reconstruct pipe-axis values where required.
- Updated cylindrical `q/r` and radial derivative handling to reconstruct
  centreline values by azimuthal projection rather than by the previous
  lower-order estimates.
- Updated conservative/primary velocity refresh order in the momentum solver so
  thermal cases convert conservative variables back to velocities before
  enforcing velocity boundary conditions.
- Changed cylindrical visualisation binary output to little-endian and updated
  the generated XDMF metadata accordingly.
- Tidied module `use` lists across the Fortran sources: unused modules were
  removed and remaining imports were alphabetically ordered for readability and
  easier review.

### Fixed

- Fixed thermal restart handling so restart cases do not require the inlet
  thermal boundary condition checks that apply only to fresh initialisation.
- Fixed statistics restart/read handling by allowing statistics arrays to be
  populated in `STATS_READ` mode, not only accumulated in `STATS_TAVG` mode.
- Fixed the statistics averaging count after restart by using
  `iter - dm%stat_istart` instead of adding one extra sample.
- Fixed `is_IO_off` behaviour so initial visualisation, mesh/check-file output,
  monitor/probe history files, outlet-record output, folder creation, and
  initial thermo/flow visualisation are skipped when I/O is disabled.
- Fixed the Makefile object list by removing stale merge-conflict markers around
  `eq_continuity.o`.
