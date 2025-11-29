import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Paramètres de la ligne
Y_ligne = 1j * 11.76e-4
Z_ligne = 3 + 1j * 78.63

A_ligne = 1 + Y_ligne * Z_ligne / 2
B_ligne = Z_ligne
C_ligne = Y_ligne + (Y_ligne ** 2) / 4 * Z_ligne
D_ligne = 1 + Y_ligne * Z_ligne / 2

ligne = np.array([[A_ligne, B_ligne],
                  [C_ligne, D_ligne]])
print("A_ligne est égal à:", A_ligne)

#Paramètre du transformateur
Z = 25190 * 25190 * 1j / (25190 + 25190 * 1j)
Y_transfo = 1 / Z

Z1_transfo = 0.1 + 1j * 4.03
Z2_transfo = 0.2166 + 1j * 8.71

A_transfo = 1 + Y_transfo * Z1_transfo
B_transfo = Z1_transfo + Z2_transfo + Y_transfo * Z1_transfo * Z2_transfo
C_transfo = Y_transfo
D_transfo = 1 + Y_transfo * Z2_transfo

transfo = np.array([[A_transfo, B_transfo],
                    [C_transfo, D_transfo]])

# Paramètre du système
V_ll = 230e3
V_phase = V_ll / np.sqrt(3)

#SIL par phase
SIL_total = 204.6e6 / 3

# Plage de puissance
P_start = 0.8 * SIL_total
P_end = 300e6
P_step = 2e6

P_values = np.arange(P_start, P_end + 1, P_step)

# Compensation de 0% à 70%
comp_rates = np.arange(0, 0.71, 0.10)



#Calcul de Vr
def solve_Vr(P_total, A, B, Vs_mag, prev_solution=None):
    P_phase = P_total

    def residual(vars):
        Vr_real, Vr_imag = vars
        Vr = complex(Vr_real, Vr_imag)
        Vr_mag = abs(Vr)

        if Vr_mag < 1e-9:
            Vr_mag = 1e-9

        Ir_mag = P_phase / Vr_mag
        Ir = Ir_mag * (Vr / Vr_mag)

        Vs_calc = A * Vr + B * Ir
        return [Vs_calc.real - Vs_mag, Vs_calc.imag]

    # Hypothèse
    if prev_solution is not None:
        guess = [prev_solution.real, prev_solution.imag]
    else:
        guess = [Vs_mag * 0.95, 0]

    sol, infodict, ier, msg = fsolve(residual, guess, full_output=True)

    #si échec
    if ier != 1:
        guess2 = [Vs_mag * 0.9, -Vs_mag * 0.1]
        sol, infodict, ier, msg = fsolve(residual, guess2, full_output=True)


        if ier != 1:

            return None, None

    Vr = complex(sol[0], sol[1])
    Ir = (P_phase / abs(Vr)) * (Vr / abs(Vr))

    return Vr, Ir



results = []

for comp in comp_rates:

    zcomp = -comp * (1j * 78.63)
    comp_mat = np.array([[1, zcomp],
                         [0, 1]])

    ABCD = transfo @ ligne @ comp_mat
    A_tot = ABCD[0, 0]
    B_tot = ABCD[0, 1]

    prev_Vr = None

    for P in P_values:
        Vr, Ir = solve_Vr(P, A_tot, B_tot, V_phase, prev_Vr)


        if Vr is None or Ir is None:
            continue

        prev_Vr = Vr

        Vr_mag = abs(Vr)
        VB3_pu = Vr_mag / V_phase
        var_pct = abs(V_phase - Vr_mag) / V_phase * 100
        Il = abs(Ir)

        # Calcul du déphasage entre Vs et Vr
        Vs_calc = A_tot * Vr + B_tot * Ir
        phase_Vs = np.angle(Vs_calc, deg=True)
        phase_Vr = np.angle(Vr, deg=True)
        dphi = phase_Vs - phase_Vr  # en degrés


        results.append([comp * 100, P / 1e6, VB3_pu, var_pct, Il, dphi])


if len(results) == 0:
    raise RuntimeError("Aucune solution trouvée : vérifier les paramètres (SIL, V_phase, etc.).")

#Conversion en tableau
res = np.array(results)
comp_col = res[:, 0]
P_MW = res[:, 1]
VB3_pu = res[:, 2]
var_pct = res[:, 3]
Il = res[:, 4]
dphi = res[:, 5]


#Graphique

# ===== Figure 1 - Variation de tension =====
plt.figure()
for comp in comp_rates:
    idx = comp_col == (comp * 100)
    if np.any(idx):
        plt.plot(P_MW[idx], var_pct[idx], linewidth=2, label=f"{int(comp * 100)}%")

plt.title("Variation de tension |(VB2 - VB3) / VB2| (%)", fontsize=12, fontweight='bold')
plt.xlabel("Puissance P (MW)")
plt.ylabel("Variation (%)")
plt.grid(True, alpha=0.3)
plt.axhline(5, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()

# ===== Figure 2 - Courant de ligne =====
plt.figure()
for comp in comp_rates:
    idx = comp_col == (comp * 100)
    if np.any(idx):
        plt.plot(P_MW[idx], Il[idx], linewidth=2, label=f"{int(comp * 100)}%")

plt.title("Courant de ligne en fonction de P", fontsize=12, fontweight='bold')
plt.xlabel("Puissance P (MW)")
plt.ylabel("I_ligne (A)")
plt.grid(True, alpha=0.3)
plt.axhline(1500, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()

# ===== Figure 3 - VB3 (pu) =====
plt.figure()
for comp in comp_rates:
    idx = comp_col == (comp * 100)
    if np.any(idx):
        plt.plot(P_MW[idx], VB3_pu[idx], linewidth=2, label=f"{int(comp * 100)}%")

plt.title("VB3 (pu) en fonction de la puissance transportée", fontsize=12, fontweight='bold')
plt.xlabel("Puissance P (MW)")
plt.ylabel("VB3 (pu)")
plt.grid(True, alpha=0.3)
plt.axhline(1.05, linestyle='--')
plt.axhline(0.95, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()

# ===== Figure 4 - Déphasage Vs - Vr =====
plt.figure()
for comp in comp_rates:
    idx = comp_col == (comp * 100)
    if np.any(idx):
        plt.plot(P_MW[idx], dphi[idx], linewidth=2, label=f"{int(comp * 100)}%")

plt.title("Déphasage entre Vs et Vr (°)", fontsize=12, fontweight='bold')
plt.xlabel("Puissance P (MW)")
plt.ylabel("Déphasage (°)")
plt.grid(True, alpha=0.3)
plt.axhline(35, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()


# ------------ Analyse numérique : P_max pour chaque taux de compensation ------------

Itherm = 1250.0  # A par sous-conducteur (tableau A.4)
I_lim = 0.3 * Itherm  # 30 % de la limite thermique

print("\n===== RÉSUMÉ PAR TAUX DE COMPENSATION =====")
for comp in comp_rates:
    idx = comp_col == (comp * 100)

    P_comp = P_MW[idx]  # Puissances (MW)
    VB3_comp = VB3_pu[idx]  # VB3 (pu)
    var_comp = var_pct[idx]  # variation (%)
    Il_comp = Il[idx]  # courant de ligne (A)
    delta_comp = dphi[idx]  # angle δ (°)

    # Contraintes Hydro-Québec
    ok = (
            (VB3_comp >= 0.95) & (VB3_comp <= 1.05) &  # VB3 dans ±5 %
            (var_comp <= 5.0) &  # variation ≤ 5 %
            (delta_comp <= 35.0) &  # δ ≤ 35°
            (Il_comp / 4.0 <= I_lim)  # courant par sous-conducteur ≤ 30 % Itherm
    )

    if not np.any(ok):
        print(f"Taux {int(comp * 100)} % : aucune puissance ne respecte toutes les contraintes.")
        continue

    # Puissance maximale qui respecte toutes les contraintes
    P_max = P_comp[ok].max()
    i_max = np.where((P_comp == P_max) & ok)[0][0]

    print(f"\nTaux de compensation : {int(comp * 100)} %")
    print(f"  P_max admissible  = {P_max:.1f} MW")
    print(f"  VB3 (pu)          = {VB3_comp[i_max]:.3f}")
    print(f"  Variation         = {var_comp[i_max]:.2f} %")
    print(f"  I_ligne           = {Il_comp[i_max]:.1f} A  (par cond ≈ {Il_comp[i_max] / 4:.1f} A)")
    print(f"  δ                 = {delta_comp[i_max]:.2f} °")
