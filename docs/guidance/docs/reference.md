# Reference

Reference documentation is organized for quick lookup rather than sequential reading.

## User-Facing Reference

| Page | Use |
| --- | --- |
| [CHAPSim Input File Guide](input-file.md) | Variable meanings, Fortran types, IDs, and common setup rules for `input_chapsim.ini`. |
| [Benchmark and Example Cases](benchmark-cases.md) | Case naming conventions and recommended starting cases. |
| [Postprocessing and Output Data](postprocessing.md) | Output folders, statistics levels, monitor scripts, and visualisation scripts. |

## Repository File Map

| Path | Description |
| --- | --- |
| `src/` | Main Fortran solver source. |
| `tests/` | Smoke/regression cases and metric-checking tools. |
| `examples/` | Example postprocessing scripts, reference data, and case-local plotting utilities. |
| `prepost/input_generator/` | Python and shell tools for generating or modifying input files. |
| `prepost/mesh_reviewer/` | Interactive mesh-stretching inspection tool. |
| `docs/guidance/docs/` | Markdown source for user guidance. |
| `docs/guidance/html/` | Dependency-free static HTML preview generated from the Markdown source. |
| `docs/code_structure/` | Generated FORD code-structure documentation. |

## Code Structure Reference

The generated FORD documentation is separate from the user guide:

```text
docs/code_structure/index.html
```

From the documentation home page, choose **Code Structure** to browse modules,
procedures, derived types, and source files.
