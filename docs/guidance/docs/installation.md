# Installation and First Run

This page gives the minimal path from source tree to a short CHAPSim2 run.
Platform-specific compiler and MPI module names vary between local machines and
HPC systems, so use the commands here as a checklist rather than a universal
script.

## Prerequisites

CHAPSim2 is a Fortran Direct Numerical Simulation (DNS) solver with MPI parallelism and FFT-based components.
A working build environment normally needs:

- A Fortran compiler (gfortran version 10 or newer, or compatible Fortran 90 compiler)
- An MPI implementation with Fortran support (OpenMPI, MPICH, or compatible)
- Standard build tools: `make` and optionally `cmake`
- The bundled 2decomp-fft library dependency

### Platform-Specific Installation

**Debian/Ubuntu workstation:**

```bash
sudo apt-get update
sudo apt-get install gfortran make cmake openmpi-bin libopenmpi-dev
```

**macOS (using Homebrew):**

```bash
brew install gcc cmake open-mpi
```

**HPC systems:**

Prefer the site-provided compiler and MPI modules. Load them via your module system:

```bash
module load compiler/gcc
module load mpi/openmpi
```

Consult your HPC documentation for available compiler and MPI versions.


## Download

```bash
git clone https://github.com/CCP-NTH/CHAPSim2.git
cd CHAPSim2
```

## Build

From the repository root:

1. Make the build script executable:

```bash
chmod +x build_chapsim.sh
```

2. Run the build script:

```bash
./build_chapsim.sh
```

The script will prompt you for configuration options:
- Whether to refresh and rebuild the 2decomp-fft library
- Whether to run `make clean` before compilation
- Whether to build in debug mode

Upon successful completion, the compiled executable `CHAPSim` will be available in the `bin/` directory.

The compiled solver from the build tree is used when running test cases under `tests/`.

## Run Tests

### Regression and Smoke Tests

To execute comprehensive validation of all test cases:

```bash
./run_regression.sh
```

This runs automated regression testing with metrics validation. Comparison tolerances are defined in `tests/tools/tolerances.json`.

### Manual Test Case

To run a single test case manually:

```bash
cd tests/<case_name>
mpirun -np 4 ../../../bin/CHAPSim
```

Inspect `input_chapsim.ini` in the test directory to understand the case configuration.

### Quick Smoke Test

For rapid validation, limit iterations using the environment variable:

```bash
CHAPSIM_NITER=20 mpirun -np 4 ../../../bin/CHAPSim
```

## Verifying Your Installation

During a first run, monitor the solver output for:

- Input section parsing success without errors
- Appropriate Courant–Friedrichs–Lewy (CFL) and time-step values
- Global mass conservation diagnostics
- Boundary-condition configuration status
- Timely creation of restart, visualization, and statistics outputs

Review the timestamped log file for comprehensive diagnostics, warnings, and errors.

For input-file details and configuration options, see [CHAPSim Input File Guide](input-file.md).

## Next Steps

Once installation is complete:

1. Review the [Project Structure](../index.md#project-structure) to understand the repository layout
2. Examine test cases in `tests/` to learn configuration patterns
3. Use `prepost/autoinput/` for Python scripts to generate input files
4. Use `prepost/useful_scripts/run_local.sh` to set up and run your own simulations
