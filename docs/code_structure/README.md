# CHAPSim2 FORD Reference Documentation

This directory contains the FORD-generated Fortran code-structure reference and
the configuration used to regenerate it.

## Configuration

Use the single maintained FORD configuration file:

```text
docs/code_structure/ford.yaml
```

## Generate the Reference

From the repository root, generate into a temporary folder first:

```bash
repo="$(pwd)"
ford -d "$repo/src" -o /tmp/chapsim2_ford_out "$repo/docs/code_structure/ford.yaml"
```

Then copy the completed generated files into `docs/code_structure/`.

## FORD Comment Style

Use FORD documentation comments directly before the entity they describe:

```fortran
!> Short one-line summary.
!>
!> Longer explanation of purpose, assumptions, and workflow.
!> - dm (in): Domain descriptor.
!> - fl (inout): Flow state.
subroutine example_routine(fl, dm)
```

Recommended practice:

- The maintained marker configuration is `predocmark: ">"` and
  `docmark: "!"`, which means FORD recognises leading `!>` blocks and
  continuation `!!` comments.
- The alternate markers are intentionally empty because CHAPSim2 does not use a
  second documentation-comment style.
- Put a short `!> ...` summary block before each public module, derived type,
  subroutine, and function.
- Use simple Markdown bullets such as `!> - dm (in): Domain descriptor.` for
  important arguments.
- Keep ordinary implementation notes as normal `!` comments inside routines.
- Avoid long banner comments as the only documentation; FORD renders structured
  `!>` blocks much more clearly.

## View the Reference

Open the generated HTML file in a browser, for example:

```bash
google-chrome docs/code_structure/index.html
```
