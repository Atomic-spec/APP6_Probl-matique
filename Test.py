import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import cmath

# ================= CONFIGURATION =================
# System Parameters
V_ll = 230e3
V_phase = V_ll / np.sqrt(3)
SIL_total = 204.6e6 / 3  # SIL per phase

# Load Parameters
PF = 1.0  # Power Factor (1.0 = Resistive, <1.0 = Lagging)
phi_load = np.arccos(PF)  # Angle in radians

# Thermal Constraints
Itherm = 1250.0
num_conductors = 4  # Bundle size
I_lim = 0.3 * Itherm * num_conductors  # Total limit for the phase

# Range
P_start = 0.8 * SIL_total
P_end = 350e6  # Extended slightly to see voltage collapse
P_step = 2e6
P_values = np.arange(P_start, P_end + 1, P_step)
comp_rates = np.arange(0, 0.71, 0.10)

# ================= MATRICES =================
# Line Parameters
Y_ligne = 1j * 11.76e-4
Z_ligne = 3 + 1j * 78.63
A_lin = 1 + Y_ligne * Z_ligne / 2
B_lin = Z_ligne
C_lin = Y_ligne + (Y_ligne ** 2) / 4 * Z_ligne
D_lin = A_lin
mat_line = np.array([[A_lin, B_lin], [C_lin, D_lin]])

# Transformer Parameters
R_fe = 25190
X_mu = 25190
Z_shunt = (R_fe * 1j * X_mu) / (R_fe + 1j * X_mu)
Y_transfo = 1 / Z_shunt

Z1_t = 0.1 + 1j * 4.03
Z2_t = 0.2166 + 1j * 8.71

A_t = 1 + Y_transfo * Z1_t
B_t = Z1_t + Z2_t + Y_transfo * Z1_t * Z2_t
C_t = Y_transfo
D_t = 1 + Y_transfo * Z2_t
mat_transfo = np.array([[A_t, B_t], [C_t, D_t]])


# ================= SOLVER =================
def solve_system(P_target, A, B, Vs_target, prev_Vr_guess=None):
    """
    Solves for Vr (Voltage Receiving) given P_load and Vs (Voltage Sending).
    Includes Power Factor logic.
    """

    def equations(vars):
        Vr_real, Vr_imag = vars
        Vr = complex(Vr_real, Vr_imag)

        # Avoid division by zero
        if abs(Vr) < 1.0: return [1e9, 1e9]

        # Calculate Current based on Power and PF
        # S = P + jQ -> I = conj(S/V)
        # But simpler: |I| = P / (|V| * PF)
        # Angle of I = Angle of V - phi_load

        Ir_mag = P_target / (abs(Vr) * PF)
        angle_Vr = cmath.phase(Vr)
        angle_Ir = angle_Vr - phi_load  # Lagging current

        Ir = cmath.rect(Ir_mag, angle_Ir)

        # Calculate Source Voltage using A and B parameters
        Vs_calc = A * Vr + B * Ir

        # Residuals: Real and Imaginary parts of difference from known Source Mag
        # We assume Source Angle is reference, but here we just match Magnitude
        return [abs(Vs_calc) - Vs_target, 0]  # Simplified constraint logic for fsolve

    # Better approach for residual:
    # We want |Vs_calc| = Vs_target.
    # However, fsolve needs N equations for N vars.
    # Let's solve for Vr complex directly assuming Vs is reference angle 0?
    # No, Vr angle changes.

    # Revised Equations for stability
    def precise_equations(vars):
        v_r, v_i = vars
        Vr = complex(v_r, v_i)

        Ir_mag = P_target / (abs(Vr) * PF)
        Ir_angle = cmath.phase(Vr) - phi_load
        Ir = cmath.rect(Ir_mag, Ir_angle)

        Vs_calc = A * Vr + B * Ir

        # Error 1: Magnitude match
        err_mag = abs(Vs_calc) - Vs_target
        # Error 2: We can enforce Vs to be real-aligned (angle 0) for convenience,
        # or just let the solver find any valid state.
        # Ideally: Vs_calc.imag = 0 (Reference bus)
        err_angle = Vs_calc.imag

        return [err_mag, err_angle]

    guess = [Vs_target, 0] if prev_Vr_guess is None else [prev_Vr_guess.real, prev_Vr_guess.imag]

    root, info, ier, msg = fsolve(precise_equations, guess, full_output=True)

    if ier != 1: return None, None

    Vr_sol = complex(root[0], root[1])
    Ir_mag = P_target / (abs(Vr_sol) * PF)
    Ir_sol = cmath.rect(Ir_mag, cmath.phase(Vr_sol) - phi_load)

    return Vr_sol, Ir_sol


# ================= SIMULATION LOOP =================
results = []

print(f"Starting Simulation with PF = {PF}")

for comp in comp_rates:
    # Compensation Matrix (Series Capacitor)
    # Z_comp = -j * Xc
    z_comp_val = -comp * (1j * 78.63)  # 100% of line reactance
    mat_comp = np.array([[1, z_comp_val], [0, 1]])

    # Total Matrix: Source -> Transfo -> Line -> Comp -> Load
    mat_total = mat_transfo @ mat_line @ mat_comp
    A_tot = mat_total[0, 0]
    B_tot = mat_total[0, 1]

    prev_Vr = complex(V_phase, 0)

    for P in P_values:
        Vr, Ir = solve_system(P, A_tot, B_tot, V_phase, prev_Vr)

        if Vr is None:
            break  # Stop if voltage collapse (nose curve reached)

        prev_Vr = Vr  # Update guess for next step (continuation method)

        # Calculate metrics
        Vr_mag = abs(Vr)
        VB3_pu = Vr_mag / V_phase
        var_pct = (V_phase - Vr_mag) / V_phase * 100  # Voltage Regulation
        Il_mag = abs(Ir)

        # Calculate Delta (Angle difference)
        Vs_calc = A_tot * Vr + B_tot * Ir
        delta = np.degrees(np.angle(Vs_calc) - np.angle(Vr))

        results.append([comp * 100, P / 1e6, VB3_pu, var_pct, Il_mag, delta])

# ================= RESULTS PROCESSING =================
res = np.array(results)
if len(res) == 0: raise ValueError("No solution found.")

comp_col = res[:, 0]
P_MW_col = res[:, 1]
VB3_col = res[:, 2]
var_col = res[:, 3]
Il_col = res[:, 4]
delta_col = res[:, 5]

# ================= PLOTTING (Example: VB3 pu) =================
plt.figure(figsize=(10, 6))
for comp in comp_rates * 100:
    idx = np.isclose(comp_col, comp)
    if np.any(idx):
        plt.plot(P_MW_col[idx], VB3_col[idx], linewidth=2, label=f"K={int(comp)}%")

plt.title(f"Voltage Profile (VB3) vs Power Transfer (PF={PF})", fontweight='bold')
plt.xlabel("Power P (MW)")
plt.ylabel("Voltage (p.u.)")
plt.axhline(0.95, color='r', linestyle='--', label='Min Limit (0.95)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ================= CONSTRAINTS CHECK =================
print("\n===== MAXIMUM PERMISSIBLE POWER (HYDRO-QUEBEC CRITERIA) =====")
for comp in comp_rates * 100:
    idx = np.isclose(comp_col, comp)

    # Filter data for this compensation level
    P_ = P_MW_col[idx]
    V_ = VB3_col[idx]
    Var_ = var_col[idx]
    I_ = Il_col[idx]
    D_ = delta_col[idx]

    # Boolean mask for constraints
    valid_mask = (
            (V_ >= 0.95) & (V_ <= 1.05) &
            (Var_ <= 5.0) &
            (D_ <= 35.0) &
            (I_ <= I_lim)
    )

    if np.any(valid_mask):
        max_P = np.max(P_[valid_mask])
        print(f"Comp {int(comp)}%: Max Power = {max_P:.1f} MW")
    else:
        print(f"Comp {int(comp)}%: No valid operating point found.")