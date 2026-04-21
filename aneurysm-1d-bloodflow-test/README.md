# 1D Blood Flow Model (Aneurysm Project)

## Overview 
This repository contains a simplified 1D numerical model for wave propagation in blood flow, developed as part of an aneurysm-related research project. The goal is to study how waves evolve over time in a simplified arterial segment using a finite-difference scheme.<br>

The model tracks two coupled variables:

- Q(z,t): primary wave variable (e.g. flow/pressure-like quantity)
- A(z,t): auxiliary coupled variable

From these, characteristic variables are constructed:

- W1 = Q + A (right-moving waves)
- W2 = Q - A (left-moving waves)

---

## Numerical Method
The system is solved using the Lax–Wendroff scheme, a second-order finite difference method for hyperbolic PDEs.<br>

Key features:

- Spatial domain: z ∈ [0, 1]
- Uniform grid discretisation
- Explicit time stepping
- CFL-stable timestep choice (dt proportional to dz)

A damping term is included in the Q-equation to model energy loss:

- This leads to gradual energy decay over time.

---

## Initial Conditions
The simulation starts from two Gaussian perturbations:

- One in Q centered at z = 0.4
- One in A centered at z = 0.7

These evolve into travelling waves due to the hyperbolic nature of the system.

---

## Boundary Conditions
Coupled derivative-based boundary conditions are used:

- Left boundary: Q_z + A_z = 0
- Right boundary: Q_z - A_z = 0

These allow outgoing waves while minimizing artificial reflections.

---

## Outputs
The simulation generates snapshots at selected times:

### Primary variables
Q(z,t)
A(z,t)

### Characteristic variables
W1 = Q + A (right-moving component)
W2 = Q - A (left-moving component)

### Energy
The total energy is computed as: E(t) = ∫ (Q² + A²) dz

This is used to observe numerical damping and stability behaviour.

---

## Key Observations
- Waves propagate and split into left- and right-moving components
- W1 and W2 clearly separate directional wave behaviour
- Energy shows a gradual decrease due to damping
- Boundary conditions reduce strong reflections
- Small numerical dispersion is visible in wave shape over time

---

## Files
- code/bloodflow_lax_wendroff_final.py → main improved model
- code/bloodflow_lax_wendroff_boundary_conditions.py → boundary condition test
- figures/ → generated plots (Q, A, W1, W2, energy)

---

## Purpose of This Version
This version represents a cleaned and improved numerical experiment focusing on:

- stable wave propagation
- characteristic decomposition
- energy behaviour
- physically consistent boundary treatment

It serves as a foundation for further 1D vascular modelling extensions.<br>

---

## Notes
This is a simplified educational / research prototype, not a full physiological model of aneurysm blood flow. It is intended for numerical experimentation and qualitative insight into wave behaviour in 1D systems.
