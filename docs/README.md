# CHAPSim2 Documentation

CHAPSim2 documentation is organized into two components:

- **User Guidance**: Comprehensive documentation covering installation, benchmark cases, practical workflows, input reference, numerical methodology, and troubleshooting
- **Code Structure Reference**: Automatically generated FORD documentation providing Fortran API and source-code structure details

## Accessing the Documentation

Browser access to documentation:

```bash
google-chrome docs/guidance/html/index.html
```

## Documentation Regeneration

**User Guidance (MkDocs format):**

Regerate the user guidance documentation from source Markdown files (requires MkDocs):

```bash
cd docs/guidance/
mkdocs build
```

**Code Structure Reference (FORD format):**

Regenerate the Fortran API reference from source annotations:

```bash
cd docs/code_structure/
ford ford.yaml
```

**Static HTML Preview (no dependencies required):**

If MkDocs is unavailable, regenerate static HTML previews:

```bash
python3 docs/guidance/build_static.py
```
