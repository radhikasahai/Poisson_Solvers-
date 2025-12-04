# Comparing Poisson Solvers

## Overview

This project implements and compares three different numerical methods for solving Poisson's equation in 2D and 3D space. The comparison evaluates the performance, accuracy, and computational efficiency of **Jacobi Iteration**, **Successive Over-Relaxation (SOR)**, and **Finite Difference** methods.

## Motivation

Poisson solvers are essential tools in quantum chemistry calculations, providing a way to determine the electrostatic environment experienced by electrons. This is crucial for solving the Schrödinger equation and understanding the electronic structure of molecules [1].

## Background: Poisson Equation for Electrostatic Potential

### Discretization and Grid

The simulation space is discretized into a grid of points, where the electrostatic potential (φ) is represented by values assigned to each grid point.

### Grid Size Impact

- **Fine Grid**: More points with smaller spacing lead to more accurate solutions, capturing detailed potential variations
- **Computational Cost**: Finer grids increase computational expense
- **Grid Convergence**: Solutions converge when further refinement no longer produces significant changes, indicating grid-independent accuracy

### Fundamental Principle

Poisson's equation links charge distribution to electrostatic potential, enabling calculation of electrostatic fields from known charge distributions.

## Numerical Methods

### 1. Jacobi Iteration

The simplest method for solving Poisson's equation.

**Algorithm:**
- Iteratively updates the potential at each grid point
- Based on the average potential of four nearest neighbors (up, down, left, right)
- Simple update rule: `new_phi = 0.25 * (neighbor_sum + charge_term)`

**Advantages:**
- Easy to implement
- Straightforward to parallelize

**Disadvantages:**
- Slowest convergence among the three methods
- Requires many iterations to reach desired accuracy

**Implementation:** 2D and 3D working

### 2. Successive Over-Relaxation (SOR)

An acceleration technique built upon Jacobi iteration.

**Algorithm:**
- Introduces a relaxation parameter **ω** (omega) between 1 and 2
- Blends current potential value with weighted average from neighbors
- Update rule: `new_phi = (1 - ω) * old_phi + ω * 0.25 * (neighboring_potentials + charge_term)`

**Advantages:**
- Faster than Jacobi iteration
- Fewer iterations required for convergence
- Controlled acceleration through ω parameter

**Disadvantages:**
- Optimal ω is problem-dependent
- Difficult to automatically determine optimal relaxation parameter
- Requires manual tuning for different problems

**Implementation:** 2D working; 3D in progress

### 3. Finite Difference

A broader, more general method applicable to various partial differential equations.

**Algorithm:**
- Approximates derivatives using differences between potential values at neighboring grid points
- Iterates over interior grid points (excluding boundaries)
- Updates potential based on average of four nearest neighbors
- Spreads charge influence to surrounding points

**Advantages:**
- Fastest method among the three
- Most widely applicable (general PDE solver)
- Minimal iterations required for convergence

**Disadvantages:**
- More complex implementation
- Requires careful boundary condition handling

**Implementation:** 2D working; 3D in progress

## Test Case: Water Molecule (H₂O)

The comparison uses a water molecule geometry as test input:

| Atom | Element | X (Å) | Y (Å) | Z (Å) |
|------|---------|-------|-------|-------|
| 1    | O       | 0     | 0     | -0.9584 |
| 2    | H       | 0.8085| 0     | -0.4792 |
| 3    | H       | -0.8085| 0    | -0.4792 |

The water molecule is positioned with the oxygen atom near the origin and hydrogen atoms placed symmetrically along the x-axis at approximately 0.8 Å from the oxygen, with all atoms slightly below the xy-plane.

## Performance Comparison

### Convergence Time Results (ms)

| Method | Average Time |
|--------|--------------|
| Finite Difference (2D) | 152.463 ms |
| Jacobi Iteration (2D) | 642.882 ms |
| SOR (2D) | 286.951 ms |

### Performance Summary

| Criterion | Jacobi | SOR | Finite Difference |
|-----------|--------|-----|-------------------|
| **Speed** | Slower | Medium | **Fastest** |
| **Ease of Implementation** | **Easiest** | Moderate | Most Complex |
| **Parameter Tuning** | None | Required (ω) | None |
| **Iteration Count** | Highest | Lower | **Lowest** |

**Key Finding**: Finite Difference method is ~4.2× faster than Jacobi and ~1.9× faster than SOR for 2D problems.


## Future Work

1. **3D Implementation**: Complete 3D implementations of SOR and Finite Difference methods
2. **Generalization**: Extend code to handle more complex molecular inputs beyond water
3. **Parallelization**: Parallelize solvers for improved performance on multi-core systems
4. **Validation**: Implement convergence analysis and numerical validation against analytical solutions

## References

[1] Herbert, J. M. (2021). Dielectric continuum methods for quantum chemistry. *WIREs Computational Molecular Science*, 11(4). https://doi.org/10.1002/wcms.1519

[2] Poisson's Equation. Retrieved from https://acme.byu.edu/00000179-d3f1-d7a6-a5fb-ffff6a1c0001/poissonequation-1-pdf

## Installation & Usage

To use this project, ensure you have the necessary dependencies installed and run the solvers with your desired molecular geometry as input.

```bash
# Example usage (implementation-specific)
python poisson_solver.py --method finite_difference --dimension 2d --input molecule.xyz
```



## Author

Radhika Sahai
radhikasahai786@gmail.com