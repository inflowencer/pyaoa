# pyaoa - Quickstart

### Boundary conditions:

| Symbol | Name (unit) | Description |
| --: | ------ | --------------------- |
| Re | Reynolds number (-) | Free stream Reynolds number $\text{Re} = \frac{u L \rho}{\mu}$ |
| d | 2D Object depth (m) | For calculating the forces and area in 2D simulations |
| T | Free stream temperature (K) | Static temperature in the free stream |
| p | Free stream pressure (Pa) | Static pressure in the free stream |
| I | Turbulent intensity (%) | [Free stream turbulent intensity](https://www.cfd-online.com/Wiki/Turbulence_intensity) |
| inlet | Inlet name (-) | Specify the name and type of the farfield inlet BC as defined in the mesh for correct adjustment of the velocity components |

### Numerics:

| Key | Value (options) | Description |
| --: | ----| ----------- |
| solver | openfoam (su2, fluent) | Set the CFD solver |
| dim | 2d (3d) | Dimension of the simulation |
| type | incompressible (compressible, transonic, supersonic) | Type of simulation; required for setting farfield boundary conditions |
| iter | 200 | Number of iterations |
| np | 4 | Number of processors used for simulation; required for the run-script. |

### Parameters

### I/O

### Objects

### Plot