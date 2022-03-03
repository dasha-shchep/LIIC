# Linear interpolation of internal coordinates in Julia
Execute by running `julia liic.jl --flags path_to_geom1 path_to_geom2`

Flags:

`--cartesian`/`-c` LIIC in cartesian coordinates: default

`--internal`/`-i` LIIC in internal coordinates

`--distance`/`-d` LIIC using distance matrix

`--steps`/`-s` Number of steps in LIIC

Internal coordinate scan requires a working installation of OpenBabel, for now.

A Kabsch optimisation routine can preceed the LIIC and is recommended in the case of cartesian coordinates, and also gives the RMSD between the two structures.
It is invoked by the `--kabsch`/`-k` flag.


