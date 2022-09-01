# pyaoa - Automatic Angle-of-Attack Analysis

Automatic pre- and postprocessing of external aerodynamic simulations of arbitrary numbers of objects to analyze at different angles-of-attack $\alpha$ in your preferred CFD solver.

## Installation

```sh
pip install pyaoa
```

<!-- # Requirements
pip install numpy pandas matplotlib rich pyyaml
# pyaoa -->

## Examples

Clone the repository which includes the examples:

```sh
git clone https://github.com/inflowencer/pyaoa.git ~/
```

Four examples are provided in the [examples](examples/) directory.

1. [Subsonic NACA0012](examples/NACA0012/) in SU2
2. [Transonic RAE2822](examples/RAE2822/) in OpenFOAM
3. [Subsonic Selig S1210](examples/S1210) in Fluent
4. [High-lift airfoil variations](examples/NACA0056mod) in SU2

## Userguide

See the [Quickstart](doc/README.md).