# Dissertation Code: Stokes Flow Simulations of Fractal Aggregates in MATLAB

This repository contains MATLAB scripts used in my dissertation project, which investigates drag on self-similar fractal aggregates of spheres in Stokes flow, in the continuum regime. The project uses two numerical methods: the **Method of Fundamental Solutions (MFS)** and the **Method of Reflections (MOR)**, supported by various utility scripts for efficient matrix construction and sphere positioning.

---

## Acknowledgment and Citation

The work presented here was conducted in collaboration with, and under the guidance of, my supervisor, **Duncan A. Lockerby**. I have his permission to share this code. The methodologies implemented build upon insights from the following paper, which you are encouraged to cite if extending this work:

**Jordan, J.J.P., & Lockerby, D.A. (2024).** *The method of fundamental solutions for multi-particle Stokes flows: Application to a ring-like array of spheres.* Journal of Computational Physics. [https://doi.org/10.1016/j.jcp.2024.111235](https://doi.org/10.1016/j.jcp.2024.111235)

---

## Methods Overview

### Method of Fundamental Solutions (MFS)
The **MFS** is a boundary method that can be applied to Stokes flow due to the property of linearity of the equations that
arises from neglecting inertial forces in the fluid. To determine the unique solution for different geometries, boundary conditions are satisfied using a superposition of fundamental solutions (Stokeslets). It is a mesh-free method, requiring discretisation only on the surface of the geometry. 

- **Advantages**: Significant reduction in computational cost and dimensionality.
- **Challenges**: Ill-conditioning of the MFS matrix if site positions are not carefully chosen.
- **Usage**: Ideal for single-particle drag calculations or small systems.

### Method of Reflections (MOR)
The **MOR** extends the MFS to simulate multiple particles by iteratively calculating hydrodynamic interactions. It uses smaller, more manageable matrices for each particle, iteratively improving accuracy until convergence.

- **Advantages**: Computational efficiency for multi-particle systems.
- **Extensions**: Implements relaxation factors and cutoffs to balance accuracy and speed.
- **Usage**: Recommended for large systems with many interacting particles.

---

## Key Scripts and Their Purposes

### Shared Utility Scripts
- **`getSpherePositions.m`**: Generates positions of spheres for aggregates. This script is critical for defining the spatial layout of aggregates, such as chains, squares, or fractals. The generated sphere positions must be updated for each aggregate you wish to simulate.
- **`pointonsphere.m`**: Computes evenly distributed points on a sphere for collocation or singularity sites.
- **`constructRHS.m`**: Constructs the right-hand side of the linear system for Stokes flow simulations.
- **`matrixConstructFast.m`**: Efficiently assembles the system matrix.

### Workflow Scripts
#### MFS
- **`MFSv1.m`**: Main script for calculating drag on single particles using the Method of Fundamental Solutions.

#### MOR
- **`MORv4.m`**: The most up-to-date script for simulating drag on multiple interacting particles using the Method of Reflections.
- **Other versions**: Explore earlier iterations of the MOR algorithm for comparison.

---

### Important Note on Aggregate Simulations
To run simulations for different aggregates using the MOR workflow:
1. Update the sphere positions in **`getSpherePositions.m`** to define the layout of your desired aggregate (e.g., chain, square, or fractal geometries).
2. The updated sphere positions are then used in MOR scripts (e.g., `MORv4.m`) to calculate drag and other properties for the specified aggregate.
