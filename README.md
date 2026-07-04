# Aerofoil Generator

Python scripts to generate aerofoil definitions.

## Scripts

- `bin/naca4series.py` creates NACA 4-series aerofoils.
- `bin/naca5series.py` creates NACA 5-series aerofoils, including standard
  mean lines such as `23012` and reflex mean lines such as `23112`.

Both scripts write CSV coordinates to stdout by default:

```sh
bin/naca5series.py 23012 > naca23012.csv
```

Use `-h` or `--help` on either script to see the available resolution,
chord-length, point-spacing, plane, z-coordinate, and mean-camber-line options.
