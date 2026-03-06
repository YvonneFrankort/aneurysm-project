March 2026

Implemented test problem from manuscript.

System:

Qt + Az + (1/5)Q = 0
At + Qz = 0

Methods tried:
- Lax-Friedrichs
- Lax-Wendroff

Numerical schemes tested:
    - Lax–Friedrichs (baseline, diffusive)
    - Lax–Wendroff (second-order, sharper)
    - Lax–Wendroff with corrected boundary conditions

Observations:
    - Lax–Wendroff produces significantly sharper wave profiles.
    - LF is much more diffusive, as expected.
    - Boundary conditions implemented using finite differences.
    - Wave propagation behaves as expected (left/right-moving waves, damping in Q).
