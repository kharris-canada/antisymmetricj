# Contributing

Use a local editable install with development dependencies:

```bash
python -m pip install -e ".[dev]"
```

Before committing, run:

```bash
ruff check .
black .
isort .
mypy src
pytest
python -m build
```

Please keep scientific convention changes explicit in code comments, tests, and
`docs/theory.md`.
