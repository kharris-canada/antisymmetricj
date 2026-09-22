# Examples

These commands assume you are running from the repository root:

```bash
cd /Users/krisharris/Documents/AntisymmetricJ/Code/antisymmetricj
```

If the package has not been installed in editable mode, prefix the commands with
`PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli`. If it has been
installed, the shorter `antisymmetricj` command is equivalent.

## Calculate an f_zy Distribution

This uses the default in-memory ZCW orientation set with 28,656 orientations and
writes only the distribution table:

```bash
PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli calc_fzy \
  examples/output/fzy_distribution.txt
```

To write a plot at the same time, use `--plot_fzy`:

```bash
PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli calc_fzy \
  examples/output/fzy_distribution.txt \
  --plot_fzy examples/output/fzy_distribution.jpg
```

For an f_zy squared distribution, add `--squared`:

```bash
PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli calc_fzy \
  examples/output/fzy_squared_distribution.txt \
  --plot_fzy examples/output/fzy_squared_distribution.jpg \
  --squared
```

## Generate and Reuse a ZCW Orientation Set

Write an orientation file:

```bash
PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli make_zcw_set \
  examples/output/orientations_outfile.txt \
  --min-orientations 28656
```

Then use that file as the input to `calc_fzy`:

```bash
PYTHONPATH=src .venv/bin/python -m antisymmetricj.cli calc_fzy \
  examples/output/orientations_outfile.txt \
  examples/output/fzy_distribution_from_file.txt \
  --plot_fzy examples/output/fzy_distribution_from_file.jpg
```

## Reproduce Figure 2

The Figure 2 reproduction is intentionally an explicit Python example rather
than a CLI mode:

```bash
PYTHONPATH=src .venv/bin/python examples/reproduce_figure2.py
```

It writes the Figure 2 text tables and JPEGs into `examples/output/`.
