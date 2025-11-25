import numpy as np
import cmath
from sympy import symbols, solve

# =============================================================================
# DONNÉES GLOBALES ET CONSTANTES (Sources: Annexe 1 & 4)
# =============================================================================
CONF = {
    'V_nom_Ligne': 230000.0,  # Tension nominale ligne (230 kV) [Source: 29]
    'L_km': 240.0,  # Longueur ligne [Source: 32]
    'f': 60.0,  # Fréquence
    'epsilon0': 8.854e-12,
    'mu0': 4 * np.pi * 1e-7,

    # Géométrie Pylône [Source: Annexe 1]
    'D_ab': 12.0,
    'D_bc': 12.0,
    'D_ac': 24.0,

    # Conducteur ACSR 54/19 Martin 1351 kcmil (Table A.4 Power System)
    'd_faisceau': 0.45,  # Espacement faisceau 45cm [Source: 152]
    'nb_cond': 4,  # 4 conducteurs par phase
    'r_cond': 0.01809,  # Rayon conducteur (m) (Diam 1.424 in)
    'GMR_cond': 0.0142,  # GMR un seul brin (m) (0.0466 ft)
    'R_ac_ohm_km': 0.044,  # Résistance AC à 50°C (approx)

    # Générateur et Transfo (Estimations pour l'impédance totale)
    # Ces valeurs doivent être ajustées selon vos données exactes de laboratoire/simulink
    'S_base': 100e6,  # Base 100 MVA
    'Z_gen_pu': complex(0, 0.2),  # Réactance synchrone approx
    'Z_tx_pu': complex(0, 0.1),  # Réactance fuite transfo approx
    'Eg_module': 235000.0  # Tension interne Générateur (ramenée au secondaire)
}


# =============================================================================
# FONCTIONS UTILITAIRES
# =============================================================================
def get_GMD(dab, dbc, dac):
    """Calcule la Distance Géométrique Moyenne (GMD) [Source: 478]"""
    return (dab * dbc * dac) ** (1 / 3)


def get_GMR_faisceau(gmr_brin, d, n):
    """Calcule le GMR équivalent pour un faisceau de 4 conducteurs [Source: 474]"""
    if n == 4:
        return 1.091 * (gmr_brin * (d ** 3)) ** (1 / 4)
    return gmr_brin  # Fallback


def get_r_eq_capa(r_brin, d, n):
    """Calcule le rayon équivalent pour la capacité (faisceau 4) [Source: 603]"""
    if n == 4:
        return 1.091 * (r_brin * (d ** 3)) ** (1 / 4)
    return r_brin


# =============================================================================
# MANDAT 1 : CARACTÉRISTIQUES DE LA LIGNE
# =============================================================================
def mandat_1_parametres_ligne():
    print("\n" + "=" * 60)
    print("MANDAT 1 : CALCUL DES PARAMÈTRES DE LA LIGNE")
    print("=" * 60)

    w = 2 * np.pi * CONF['f']

    # 1. Calculs Géométriques
    GMD = get_GMD(CONF['D_ab'], CONF['D_bc'], CONF['D_ac'])

    # Pour Inductance (utilise GMR du brin)
    GMR_L = get_GMR_faisceau(CONF['GMR_cond'], CONF['d_faisceau'], CONF['nb_cond'])

    # Pour Capacité (utilise rayon réel du brin)
    r_eq_C = get_r_eq_capa(CONF['r_cond'], CONF['d_faisceau'], CONF['nb_cond'])

    # 2. Paramètres linéiques (Par phase) [Source: 480, 552]
    L_per_m = 2e-7 * np.log(GMD / GMR_L)
    C_per_m = (2 * np.pi * CONF['epsilon0']) / np.log(GMD / r_eq_C)
    R_per_m = CONF['R_ac_ohm_km'] / 1000.0 / CONF['nb_cond']  # Résistance divise par 4 brins

    # 3. Paramètres Totaux (Modèle Pi Nominal)
    Z_total = complex(R_per_m * CONF['L_km'] * 1000, w * L_per_m * CONF['L_km'] * 1000)
    Y_total = complex(0, w * C_per_m * CONF['L_km'] * 1000)  # G négligé [Source: 41]

    print(f"GMD = {GMD:.4f} m")
    print(f"GMR (Inductance) = {GMR_L:.4f} m")
    print(f"L linéique = {L_per_m * 1e3:.6f} mH/km")
    print(f"C linéique = {C_per_m * 1e9:.6f} nF/km")
    print(f"R linéique (faisceau) = {R_per_m * 1e3:.6f} Ohm/km")
    print("-" * 20)
    print(f"Z série total (R + jX) = {Z_total:.4f} Ohm")
    print(f"Y shunt total (G + jB) = {Y_total:.4e} Siemens")

    return Z_total, Y_total, L_per_m, C_per_m


# =============================================================================
# MANDAT 2 : SIL ET VALIDATION
# =============================================================================
def mandat_2_calcul_SIL(Z_line, Y_line, L_pm, C_pm):
    print("\n" + "=" * 60)
    print("MANDAT 2 : CALCUL ET VALIDATION DU SIL")
    print("=" * 60)

    # 1. Calcul du SIL [Source: 61]
    # Zc (Impédance caractéristique sans perte pour SIL approx)
    Zc = np.sqrt(L_pm / C_pm)
    SIL_W = (CONF['V_nom_Ligne'] ** 2) / Zc

    print(f"Impédance Caractéristique Zc = {Zc:.2f} Ohms")
    print(f"SIL (Puissance Naturelle)    = {SIL_W / 1e6:.2f} MW")

    # 2. Validation (Calcul de Vb3 si Vb2 = 230kV et Charge = SIL)
    # Matrice ABCD de la ligne seule (Pi Nominal) [Source: 200]
    A_L = 1 + (Z_line * Y_line) / 2
    B_L = Z_line
    # On assume VB2 = 230kV (Ligne-Ligne) -> VB2_pn = 230/sqrt(3)
    V_B2_ref = CONF['V_nom_Ligne'] / np.sqrt(3)

    # Courant SIL (Charge résistive, FP=1)
    # Approximation: Si Vb3 est proche de Vnom (ce qu'on veut prouver)
    # Pour une validation rigoureuse, on inverse la matrice:
    # [Vs] = [ABCD] [Vr]
    # Pour le SIL, théoriquement Vs = Vr (si sans perte).
    # Calculons Vb3 en supposant Vb2 fixé.

    print("\n--- Validation : Tension Vb3 lorsque P = SIL ---")
    # On utilise une méthode simplifiée : Au SIL, la ligne est "plate".
    # Calculons la chute approximative.

    # Equation: V_B2 = A*V_B3 + B*(SIL / (3*V_B3))
    # A*V_B3^2 - V_B2*V_B3 + B*P_phase = 0
    P_phase_SIL = SIL_W / 3

    # Coefficients pour ax^2 + bx + c = 0
    coeff_a = A_L
    coeff_b = -V_B2_ref
    coeff_c = B_L * (P_phase_SIL / V_B2_ref)  # Approx courant constant pour l'estimation

    # Racines
    delta = coeff_b ** 2 - 4 * coeff_a * coeff_c
    root1 = (-coeff_b + cmath.sqrt(delta)) / (2 * coeff_a)
    root2 = (-coeff_b - cmath.sqrt(delta)) / (2 * coeff_a)

    V_B3_est = abs(root2) * np.sqrt(3)  # On prend la racine logique

    print(f"Si on injecte le SIL ({SIL_W / 1e6:.1f} MW) :")
    print(f"Tension VB3 calculée ≈ {V_B3_est / 1000:.2f} kV")
    print("Note: Devrait être très proche de 230 kV (profil plat).")

    return SIL_W


# =============================================================================
# MANDAT 3 : OPTIMISATION (COMPENSATION SÉRIE)
# =============================================================================
def get_matrices_ABCD_globales(Z_line, Y_line, k_comp):
    """
    Construit la matrice ABCD totale (Gen -> Charge).
    Intègre Gen (Thévenin), Transfo, et Ligne compensée.
    k_comp: Taux de compensation série (0.0 à 0.7)
    """
    # 1. Impédance de compensation série [Source: 83]
    # Xc = k * X_ligne. On soustrait jXc de l'impédance série de la ligne.
    X_L_total = Z_line.imag
    Z_comp = complex(0, -1 * k_comp * X_L_total)
    Z_line_comp = Z_line + Z_comp  # Z_serie réduite

    # 2. Matrice Ligne (Pi Nominal) [Source: 200]
    A_L = 1 + (Z_line_comp * Y_line) / 2
    B_L = Z_line_comp
    C_L = Y_line * (1 + (Z_line_comp * Y_line) / 4)
    D_L = A_L
    Mat_Ligne = np.array([[A_L, B_L], [C_L, D_L]])

    # 3. Matrice Générateur + Transfo (Impédances séries)
    # Note: Le guide demande d'inclure tout le réseau.
    # Z_amont = Z_gen + Z_transfo (ramenés en Ohm au niveau 230kV)
    Z_base_230 = (CONF['V_nom_Ligne'] ** 2) / CONF['S_base']
    Z_amont_ohm = (CONF['Z_gen_pu'] + CONF['Z_tx_pu']) * Z_base_230

    # Matrice Série Amont
    Mat_Amont = np.array([[1, Z_amont_ohm], [0, 1]])

    # 4. Matrice Totale
    Mat_Total = np.dot(Mat_Amont, Mat_Ligne)

    return Mat_Total


def mandat_3_optimisation(Z_line, Y_line, SIL_W):
    print("\n" + "=" * 60)
    print("MANDAT 3 : OPTIMISATION (BOUCLES DE PUISSANCE ET COMPENSATION)")
    print("=" * 60)
    print("Résolution des équations non-linéaires avec SymPy [Source: 1138]")

    # Plages de calcul
    taux_comp_list = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]  # 0% à 70%

    # Variable symbolique pour Sympy (VB3 phase-neutre)
    x = symbols('x', positive=True, real=True)

    # En-tête tableau
    print(
        f"{'Comp(%)':<8} | {'P_trans(MW)':<12} | {'V_B3(kV)':<10} | {'Delta(deg)':<10} | {'Var V(%)':<10} | {'I_B3(A)':<10}")
    print("-" * 75)

    Eg = CONF['Eg_module'] / np.sqrt(3)  # Phase-Neutre référence

    for k_comp in taux_comp_list:
        # Matrice ABCD pour ce taux de compensation
        Mat_T = get_matrices_ABCD_globales(Z_line, Y_line, k_comp)
        A_tot = Mat_T[0, 0]
        B_tot = Mat_T[0, 1]

        # Définition de la plage de puissance (0.8 SIL à Pmax approx)
        # On va jusqu'à 1.5 SIL pour tester
        puissances_MW = np.arange(0.8 * SIL_W / 1e6, 2.0 * SIL_W / 1e6, 50)

        for P_MW in puissances_MW:
            P_W = P_MW * 1e6
            # Courant I_B3 = P / (3 * V_B3) car FP=1
            # Expression vectorielle : Eg = A*Vb3 + B*Ib3
            # Ici Vb3 est la référence (angle 0). Ib3 est en phase (FP=1).
            # Eg^2 = |A*x + B*(P/3x)|^2

            # Préparation des constantes complexes pour l'équation Sympy
            # Forme: Eg^2 = | (ar + j*ai)*x + (br + j*bi)*(P/3x) |^2
            # Eg^2 = (Real_part)^2 + (Imag_part)^2

            P_3 = P_W / 3.0

            # Equation selon "Guide pour mandat 3" [Source: 1148]
            # (aa - ((bb*x) - cc/x)**2 - ((dd*x) + ee/x)**2) = 0
            # Note: L'équation du guide semble être Eg^2 - Real^2 - Imag^2 = 0

            # Partie Réelle de (A*x + B*P/3x)
            # = Ar*x + Br*P/3x
            Term1_R = float(np.real(A_tot))
            Term2_R = float(np.real(B_tot) * P_3)

            # Partie Imaginaire
            Term1_I = float(np.imag(A_tot))
            Term2_I = float(np.imag(B_tot) * P_3)

            eq = (Eg ** 2) - (Term1_R * x + Term2_R / x) ** 2 - (Term1_I * x + Term2_I / x) ** 2

            # Résolution numérique (plus rapide que solve formel dans une boucle)
            # On transforme en polynôme X^4 pour éviter la lenteur de Sympy dans une boucle
            # Multiplier tout par x^2 pour enlever division, puis encore pour carré...
            # Mais suivons la demande d'utiliser solve/nsolve si possible, ou polyroots pour vitesse.
            # Pour respecter la consigne "Code Python pour calculer... avec solve", on utilise sympy
            # mais attention, c'est lent. Je vais utiliser nsolve ou conversion poly pour efficacité.

            # --- Approche Optimisée (Racines Polynôme) ---
            # L'équation est de la forme: K - (C1*x + C2/x)^2 - (C3*x + C4/x)^2 = 0
            # Multiplions par x^2: K*x^2 - (C1*x^2 + C2)^2 - ...
            # C'est un polynôme en x^4 (ou Y=x^2).

            try:
                # On utilise solve de sympy comme demandé explicitement dans le guide
                solutions = solve(eq, x)

                # Filtrage solutions réelles positives [Source: 1113]
                sol_reelles = [float(s) for s in solutions if s.is_real and s > 0]

                if not sol_reelles:
                    continue  # Pas de solution (instabilité)

                # Sélection de la bonne solution (proche de 230kV/sqrt(3) = 132kV)
                # Le guide dit de vérifier avec 0.8 SIL.
                v_target = 230000 / np.sqrt(3)
                best_v = min(sol_reelles, key=lambda v: abs(v - v_target))

                # Calculs des paramètres résultants
                V_B3_final = best_v
                V_B3_LL = V_B3_final * np.sqrt(3)

                I_B3_final = P_3 / V_B3_final

                # Calcul de VB2 (Entrée ligne) pour Delta et Var
                # VB2 = A_L * VB3 + B_L * IB3 (Matrice Ligne SEULEMENT pour VB2)
                # Attention: A_tot incluait le générateur. Ici on veut VB2 de la ligne.
                # Recalcul ABCD ligne seule
                X_L = Z_line.imag
                Z_line_k = Z_line - complex(0, k_comp * X_L)
                A_L = 1 + (Z_line_k * Y_line) / 2
                B_L = Z_line_k

                V_B2_cplx = A_L * V_B3_final + B_L * I_B3_final
                V_B2_LL = abs(V_B2_cplx) * np.sqrt(3)

                # Angle Delta (Vb2 vs Vb3)
                delta_rad = cmath.phase(V_B2_cplx) - cmath.phase(complex(V_B3_final, 0))
                delta_deg = np.degrees(delta_rad)

                # Variation de tension
                var_pct = abs((V_B2_LL - V_B3_LL) / V_B2_LL) * 100

                # Affichage si conditions OK (Exemple: V dans +/- 5% et Delta < 35)
                if 0.95 * 230000 <= V_B3_LL <= 1.05 * 230000:
                    print(
                        f"{k_comp * 100:<8.0f} | {P_MW:<12.1f} | {V_B3_LL / 1000:<10.2f} | {delta_deg:<10.2f} | {var_pct:<10.2f} | {I_B3_final:<10.1f}")

            except Exception as e:
                pass  # Erreur de convergence ou math


# =============================================================================
# EXÉCUTION PRINCIPALE
# =============================================================================
if __name__ == "__main__":
    # 1. Mandat 1
    Z_line, Y_line, L_pm, C_pm = mandat_1_parametres_ligne()

    # 2. Mandat 2
    SIL_W = mandat_2_calcul_SIL(Z_line, Y_line, L_pm, C_pm)

    # 3. Mandat 3
    mandat_3_optimisation(Z_line, Y_line, SIL_W)