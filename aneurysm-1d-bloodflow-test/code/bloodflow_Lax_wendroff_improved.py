import numpy as np
import matplotlib.pyplot as plt

# Spatial grid
N = 600
z = np.linspace(0.0, 1.0, N+1)
dz = z[1] - z[0]

# Time settings
dt = 0.0005
T_end = 2.6

snapshot_times = [0.25, 0.5, 1.0, 1.7, 2.6]
snap_Q = {t: None for t in snapshot_times}
snap_A = {t: None for t in snapshot_times}

# Initial bumps
eps = 0.02
def bump(z, center, eps):
    return np.exp(-((z - center)/eps)**2)

Q = bump(z, 0.4, eps)
A = bump(z, 0.7, eps)

t = 0.0

while t < T_end + 1e-12:

    Q_new = Q.copy()
    A_new = A.copy()

    # ---- Lax–Wendroff update ----
    Q_new[1:-1] = (
        Q[1:-1]
        - dt/(2*dz) * (A[2:] - A[:-2])
        + (dt**2)/(2*dz**2) * (Q[2:] - 2*Q[1:-1] + Q[:-2])
        - dt * Q[1:-1] / 5.0
    )

    A_new[1:-1] = (
        A[1:-1]
        - dt/(2*dz) * (Q[2:] - Q[:-2])
        + (dt**2)/(2*dz**2) * (A[2:] - 2*A[1:-1] + A[:-2])
    )

    # ---- Boundary conditions ----

    # Left boundary: Qz + Az = 0  ->  Q0 + A0 = Q1 + A1
    Q_new[0] = Q_new[1]
    A_new[0] = Q_new[1] + A_new[1] - Q_new[0]

    # Right boundary: Qz - Az = 0 ->  QN - AN = Q(N-1) - A(N-1)
    Q_new[-1] = Q_new[-2]
    A_new[-1] = Q_new[-1] - (Q_new[-2] - A_new[-2])

    Q, A = Q_new, A_new
    t += dt

    for s in snapshot_times:
        if snap_Q[s] is None and t >= s:
            snap_Q[s] = Q.copy()
            snap_A[s] = A.copy()

print("Snapshots stored:")
for s in snapshot_times:
    print(s, "OK" if snap_Q[s] is not None else "MISSING")

# ---- Plot Q ----
plt.figure(figsize=(10,6))
for s in snapshot_times:
    plt.plot(z, snap_Q[s], label=f"Q at t={s}")
plt.legend()
plt.grid()
plt.xlabel("z")
plt.ylabel("Q")
plt.title("Lax–Wendroff Solution for Q(z,t)")
plt.show()

# ---- Plot A ----
plt.figure(figsize=(10,6))
for s in snapshot_times:
    plt.plot(z, snap_A[s], label=f"A at t={s}")
plt.legend()
plt.grid()
plt.xlabel("z")
plt.ylabel("A")
plt.title("Lax–Wendroff Solution for A(z,t)")

plt.savefig("figures/Q_snapshots.png", dpi=300)
plt.savefig("figures/A_snapshots.png", dpi=300)

plt.show()