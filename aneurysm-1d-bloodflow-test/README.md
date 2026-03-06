# Numerical Simulation of a 1D Blood Flow Test Problem

This repository contains a numerical implementation of a simplified **1D blood flow model** described in the manuscript:

*Mathematical Modelling of 1D Blood Flow in Arteries with Aneurysm Using the Navier–Stokes Equations.*

The goal is to solve the simplified 1D system

Q_t + A_z + (1/5)Q = 0
A_t + Q_z = 0

on the spatial domain

z ∈ (0,1)

using finite difference schemes.

---

# Numerical Methods

Three schemes were implemented:

• Lax–Friedrichs
First‑order, highly diffusive, used as a baseline.
• Lax–Wendroff
Second‑order accurate, less diffusive, produces sharper wave profiles
• Improved Lax–Wendroff with correct boundary conditions
Final version used for the main results.

The Lax–Wendroff scheme is the primary method recommended in the manuscript and is used for the final plots.

---

# Initial Conditions

Two localized pulses are used:

Flow pulse near z = 0.4
Area pulse near z = 0.7

Both pulses are implemented as Gaussian functions with width ε = 0.02.

These pulses generate left‑ and right‑moving waves that interact through the coupled system.

---

# Boundary Conditions

The boundary conditions are

Q_z(0,t) + A_z(0,t) = 0
Q_z(1,t) − A_z(1,t) = 0

They are implemented using finite difference approximations consistent with the numerical scheme.

---

# Output Times

Solutions are evaluated at

t = 0.25
t = 0.5
t = 1.0
t = 1.7
t = 2.6

Both Q(z, t) and A(z, t) are stored and plotted.

---

# Running the Code

Install dependencies:

pip install -r requirements.txt

Run the improved solver:

python code/lax_wendroff_improved_bc.py

---

# Output

The simulation generates plots of

Q(z,t)
A(z,t)

which are saved in the `figures/` directory.

---

# Report

A short explanation of the model, numerical method, and results is included in the `report/` directory.
