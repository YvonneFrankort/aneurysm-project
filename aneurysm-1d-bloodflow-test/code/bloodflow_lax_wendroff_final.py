import numpy as np
import matplotlib.pyplot as plt
import os
os.makedirs("figures", exist_ok=True)

# Spatial grid
N = 600
z = np.linspace(0.0, 1.0, N+1)
dz = z[1] - z[0]

dt = 0.45 * dz
T_end = 2.6

snapshot_times = [0.25, 0.5, 1.0, 1.7, 2.6]
snap_Q = {t: None for t in snapshot_times}
snap_A = {t: None for t in snapshot_times}

# add energy storage for the energy-decay validation plot
energy = []   # list of (time, E)


# Initial bumps
eps = 0.02
def bump(z, center, eps):
    return np.exp(-((z - center)/eps)**2)

Q = bump(z, 0.4, eps)
A = bump(z, 0.7, eps)

energy.append((0.0, np.trapezoid(Q**2 + A**2, z)))

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
    # Left boundary: Q_z + A_z = 0
    Q_new[0] = Q_new[1] + A_new[1] - A_new[0]   # your line, correct
    A_new[0] = A_new[1]                           # your line, correct

    # Right boundary: Q_z - A_z = 0
    Q_new[-1] = Q_new[-2] + A_new[-1] - A_new[-2]  # your line, correct
    A_new[-1] = A_new[-2]                            # your line, correct

    Q, A = Q_new, A_new
    t += dt

    # record energy at each step
    energy.append((t, np.trapezoid(Q**2 + A**2, z)))

    for s in snapshot_times:
        if snap_Q[s] is None and t >= s:
            snap_Q[s] = Q.copy()
            snap_A[s] = A.copy()

print("Snapshots stored:")
for s in snapshot_times:
    print(s, "OK" if snap_Q[s] is not None else "MISSING")

for s in snapshot_times:
    print(s, np.max(snap_Q[s]), np.max(snap_A[s]))


# ---- Compute characteristic variables ----
snap_W1 = {}
snap_W2 = {}

for s in snapshot_times:
    snap_W1[s] = snap_Q[s] + snap_A[s]
    snap_W2[s] = snap_Q[s] - snap_A[s]


# ---- Plot Q ---- 
plt.figure(figsize=(10,6))
for s in snapshot_times:
    plt.plot(z, snap_Q[s], label=f"Q at t={s}")
plt.legend()
plt.grid()
plt.xlabel("z")
plt.ylabel("Q")
plt.title("Lax–Wendroff Solution for Q(z,t)")
plt.savefig("figures/Q_snapshots.png", dpi=150)
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
plt.savefig("figures/A_snapshots.png", dpi=150)
plt.show()

# ---- Plot W1 ---- 
colors = ["blue", "green", "orange", "red", "purple"]

plt.figure(figsize=(10,6))
for i, s in enumerate(snapshot_times):
    plt.plot(z, snap_W1[s], color=colors[i], label=f"W1 at t={s}")
plt.legend()
plt.grid()
plt.xlabel("z")
plt.ylabel("W1")
plt.title("Characteristic Variable W1 = Q + A")
plt.savefig("figures/W1_snapshots.png", dpi=150)
plt.ylim(-1, 1)
plt.show()

# ---- Plot W2 ---- 
plt.figure(figsize=(10,6))
for i, s in enumerate(snapshot_times):
    plt.plot(z, snap_W2[s], color=colors[i], label=f"W2 at t={s}")
plt.legend()
plt.grid()
plt.xlabel("z")
plt.ylabel("W2")
plt.title("Characteristic Variable W2 = Q - A")
plt.savefig("figures/W2_snapshots.png", dpi=150)
plt.ylim(-1, 1)
plt.show()

# Energy decay
t_arr, E_arr = zip(*energy)
plt.figure(figsize=(8,4))
plt.plot(t_arr, E_arr, color="blue", linewidth=1.2)
plt.xlabel("t")
plt.ylabel("E(t) = ∫(Q² + A²) dz")
plt.title("Energy E(t) — should decrease due to Q-damping")
plt.grid()
plt.savefig("figures/energy_decay.png", dpi=150)
plt.show()