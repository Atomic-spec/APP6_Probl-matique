import numpy as np
import cmath


def calcul_reseau_electrique():
    print("--- CALCULS DES PARAMÈTRES DE LA LIGNE ET DU SIL ---")

    # ==========================================
    # 1. DONNÉES DU PROBLÈME (Guide de l'étudiant & Annexe 1)
    # ==========================================
    # Constants
    epsilon_0 = 8.854e-12
    mu_0 = 4 * np.pi * 1e-7
    f = 60
    w = 2 * np.pi * f

    # Paramètres de la ligne (Page 5 et Annexe 1) [cite: 40, 46, 150]
    L_ligne_km = 240.0  # Longueur en km
    V_nom = 230000.0  # Tension nominale (Ligne-Ligne)

    # Géométrie du pylône (Annexe 1) [cite: 150, 478]
    # Distance entre phases (Dab = 12m, Dbc = 12m, Dac = 24m)
    D_ab = 12.0
    D_bc = 12.0
    D_ac = 24.0
    GMD = (D_ab * D_bc * D_ac) ** (1 / 3)  # Distance Géométrique Moyenne

    # Caractéristiques Conducteur (Faisceau de 4 - ACSR Martin 1351) [cite: 46, 151]
    # NOTE: Ces valeurs viennent typiquement de la table A.4 du livre "Power System".
    # J'utilise ici les valeurs standards pour un ACSR 54/19 Martin 1351 kcmil.
    d_faisceau = 0.45  # Espacement faisceau (Annexe 1) [cite: 152]
    diametre_cond_in = 1.424  # pouces (standard Martin)
    r_cond = (diametre_cond_in * 0.0254) / 2  # rayon en mètres
    GMR_brin = 0.0466 * 0.3048  # GMR d'un brin (pieds vers mètres, approx standard)
    R_ac_50C_per_km = 0.044  # Ohm/km (Valeur approx à 50C pour Martin)

    # ==========================================
    # 2. CALCUL DES PARAMÈTRES R, L, C (Annexe 1, p.22-26)
    # ==========================================

    # Calcul du GMR et r_eq pour faisceau de 4 conducteurs [cite: 474, 603]
    # Pour l'inductance (L), on utilise le GMR du brin (GMR')
    GMR_b_L = 1.091 * (GMR_brin * (d_faisceau ** 3)) ** (0.25)

    # Pour la capacité (C), on utilise le rayon réel (r) dans la formule équivalente
    GMR_b_C = 1.091 * (r_cond * (d_faisceau ** 3)) ** (0.25)

    # Inductance (L) et Capacité (C) linéiques [cite: 480, 552]
    L_per_m = 2e-7 * np.log(GMD / GMR_b_L)
    C_per_m = (2 * np.pi * epsilon_0) / np.log(GMD / GMR_b_C)

    # Valeurs totales pour la ligne de 240 km
    L_total = L_per_m * L_ligne_km * 1000
    C_total = C_per_m * L_ligne_km * 1000
    R_total = R_ac_50C_per_km * L_ligne_km

    # Impédance série (Z) et Admittance shunt (Y) totales [cite: 220, 223]
    Z_ligne = complex(R_total, w * L_total)
    Y_ligne = complex(0, w * C_total)  # On néglige la conductance G

    print(f"Paramètres calculés:")
    print(f"  Z ligne total = {Z_ligne:.4f} Ohm")
    print(f"  Y ligne total = {Y_ligne:.4e} Siemens")

    # ==========================================
    # 3. CALCUL DU SIL (Surge Impedance Loading) [cite: 58, 61]
    # ==========================================
    # Zc (Impédance caractéristique, approx sans perte pour SIL)
    Zc = cmath.sqrt(L_per_m / C_per_m)

    # SIL = V_nom^2 / |Zc| (Puissance triphasée)
    SIL_W = (V_nom ** 2) / abs(Zc)
    print(f"\n>>> SIL (Puissance Naturelle) = {SIL_W / 1e6:.2f} MW")
    print(f"    Impédance caractéristique Zc = {abs(Zc):.2f} Ohm")

    # ==========================================
    # 4. CRÉATION DES 4 MATRICES ABCD [cite: 159, 203]
    # ==========================================

    # Matrice 1: Générateur (Source Thévenin)
    # Note: L'énoncé dit de simuler le groupe par son équivalent Thévenin.
    # Ici, on met une impédance arbitraire faible car on s'intéresse surtout à la ligne B2-B3
    Z_gen = complex(0, 2.0)  # Valeur hypothétique faible ou à récupérer des données Labo
    Mat_Gen = np.array([[1, Z_gen], [0, 1]])

    # Matrice 2: Transformateur (B1 vers B2)
    # 1050 MVA, Zt standard approx 10-12% pu? On prend Zt série rapporté au secondaire.
    # Z_base_230 = 230^2 / 1050 = 50.3 Ohms. Si X=0.1pu -> X_ohm = 5 Ohms.
    Z_transfo = complex(0, 5.0)  # Valeur estimée, à ajuster selon données exactes
    Mat_Transfo = np.array([[1, Z_transfo], [0, 1]])

    # Matrice 3: Ligne de Transport (Modèle en Pi Nominal)
    # A = D = 1 + YZ/2
    # B = Z
    # C = Y(1 + YZ/4)
    A_L = 1 + (Y_ligne * Z_ligne) / 2
    B_L = Z_ligne
    C_L = Y_ligne * (1 + (Y_ligne * Z_ligne) / 4)
    D_L = A_L
    Mat_Ligne = np.array([[A_L, B_L], [C_L, D_L]])

    # Matrice 4: Matrice Équivalente Totale (Produit matriciel)
    # Note: Mat_Total = Mat_Gen * Mat_Transfo * Mat_Ligne
    Mat_Total = Mat_Gen @ Mat_Transfo @ Mat_Ligne

    print("\n--- MATRICES ABCD ---")
    print("Matrice Ligne (B2 -> B3):")
    print(f"  A = {Mat_Ligne[0, 0]:.4f}")
    print(f"  B = {Mat_Ligne[0, 1]:.4f}")
    print(f"  C = {Mat_Ligne[1, 0]:.4e}")
    print("-" * 30)

    # ==========================================
    # 5. VALIDATION DE VB3 LORSQUE P = SIL
    # ==========================================
    # On suppose que la charge consomme exactement le SIL avec Fp = 1 (résistif)
    # On fixe V_B2 = 230 kV (Tension entrée ligne) et on cherche V_B3.
    # Relation: V_B2 = A*V_B3 + B*I_B3

    print("\n--- VERIFICATION VB3 (Charge = SIL) ---")

    # Pour trouver Vb3 exact, on peut résoudre numériquement ou itérer car I_B3 dépend de Vb3.
    # Approximation: Si la ligne est sans perte, V_B3 = V_B2 = 230 kV.
    # Avec pertes (R > 0), V_B3 sera légèrement inférieur.

    # Testons si V_B3 = 230 kV est une solution proche
    V_B3_phase_ref = V_nom / cmath.sqrt(3)  # 230kV L-N reference 0 deg

    # Courant correspondant au SIL à cette tension (Charge résistive, I en phase avec V)
    I_SIL = SIL_W / (np.sqrt(3) * V_nom)

    # Calcul de V_B2 nécessaire pour avoir 230kV en sortie
    V_B2_calc_phase = Mat_Ligne[0, 0] * V_B3_phase_ref + Mat_Ligne[0, 1] * I_SIL
    V_B2_calc_LL = abs(V_B2_calc_phase) * np.sqrt(3)

    print(f"Si V_B3 = 230 kV :")
    print(f"  Tension requise à B2 (calculée) = {V_B2_calc_LL / 1000:.2f} kV")

    # Inversement: Si on fixe V_B2 = 230 kV, combien vaut V_B3?
    # V_B2 = A*V_B3 + B*(P/sqrt(3)V_B3) -> Equation quadratique complexe
    # V_B2 * V_B3 = A * V_B3^2 + B * P_phase
    # A * V_B3^2 - V_B2 * V_B3 + B * P_phase = 0

    # Résolution de l'équation quadratique pour V_B3 (Phase-Neutre)
    # Coefficients pour ax^2 + bx + c = 0
    P_phase = SIL_W / 3
    V_B2_ref = V_nom / np.sqrt(3)  # On suppose V_B2 fixé à 230kV

    # Note: C'est une approximation car V_B3 et V_B2 ont un déphasage.
    # Méthode itérative simple pour trouver V_B3 exact
    V_B3_iter = V_B2_ref
    for i in range(10):
        I_iter = P_phase / V_B3_iter  # Courant (Fp=1)
        V_B2_found = A_L * V_B3_iter + B_L * I_iter
        ratio = V_B2_ref / abs(V_B2_found)
        V_B3_iter = V_B3_iter * ratio  # Ajustement magnitude

    V_B3_final_LL = abs(V_B3_iter) * np.sqrt(3)

    print(f"\nRésultat Final :")
    print(f"Si on impose V_B2 = 230 kV (Entrée ligne) et Charge = SIL :")
    print(f">>> V_B3 (Sortie ligne) = {V_B3_final_LL / 1000:.2f} kV")
    print(f"    (Théoriquement proche de 230 kV car courbe de tension plate au SIL)")

    print("\nValidation:")
    ecart = abs(230 - V_B3_final_LL / 1000)
    if ecart < 5:
        print("SUCCESS: La tension est maintenue (~230 kV) comme attendu.")
    else:
        print("WARNING: Écart significatif, vérifier les paramètres R (pertes).")


if __name__ == "__main__":
    calcul_reseau_electrique()