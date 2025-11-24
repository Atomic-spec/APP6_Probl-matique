from APP6_Mandat3 import*
import math
import matplotlib.pyplot as plt


def solve_Vr_single(Vs, A, B, P):
    """
    Résout l'équation :
        Vs = A*Vr + B*(P/(Vr*sqrt(3)))

    => A*Vr^2 - Vs*Vr + (B*P)/sqrt(3) = 0

    Retourne UNE solution physique Vr>0, sinon None.
    """
    a = A
    b = -Vs
    c = (B * P) / math.sqrt(3)

    D = b * b - 4 * a * c  # discriminant
    if D < 0:
        return None

    sqrtD = math.sqrt(D)

    Vr1 = (-b + sqrtD) / (2 * a)
    Vr2 = (-b - sqrtD) / (2 * a)

    candidates = [Vr for Vr in (Vr1, Vr2) if Vr > 0]
    return max(candidates) if candidates else None


# --------------------------------------------------------
# PROGRAMME PRINCIPAL
# --------------------------------------------------------

Vs = 230000  # (si c'est ligne-ligne phase-> Vs = 230000/math.sqrt(3))

listA = [0.939283759772047] * 8
listB = [
    90.8550206, 83.4753272, 76.0967820, 68.7197545,
    61.3447924, 53.9727421, 46.6049858, 39.2439417
]

listP = [
    165e6, 175e6, 185e6, 195e6, 205e6, 215e6,
    225e6, 235e6, 245e6, 255e6, 265e6, 275e6, 285e6
]

# Vr_matrix[i][j] = Vr pour la paire (A[i],B[i]) à la puissance P[j]
Vr_matrix = [[None for _ in listP] for _ in listA]
Ir_matrix = [[None for _ in listP] for _ in listA]

for i, (A, B) in enumerate(zip(listA, listB)):
    for j, P in enumerate(listP):
        Vr = solve_Vr_single(Vs, A, B, P)
        Vr_matrix[i][j] = Vr
        Ir_matrix[i][j] = None if Vr is None else P / (Vr * math.sqrt(3))

print("=== Résultats Vr & Ir ===\n")
for i in range(len(listA)):
    print(f"--- Ligne A={listA[i]}, B={listB[i]} ---")
    for j in range(len(listP)):
        P = listP[j]
        Vr = Vr_matrix[i][j]
        Ir = Ir_matrix[i][j]
        if Vr is None:
            print(f"P={P:.2e} W → pas de solution réelle")
        else:
            print(f"P={P:.2e} W → Vr={Vr:.2f} V, Ir={Ir:.2f} A")
    print()

# ------------------ Graphe Vr(P) ------------------
plt.figure()
for i in range(len(listA)):
    plt.plot(listP, Vr_matrix[i], marker='o',
             label=f"A={listA[i]:.4f}, B={listB[i]:.2f}")
plt.xlabel("Puissance P (W)")
plt.ylabel("Vr (V)")
plt.title("Vr(P) pour chaque paire A,B")
plt.grid(True)
plt.legend()
plt.show()

# ------------------ Graphe Ir(P) ------------------
plt.figure()
for i in range(len(listA)):
    plt.plot(listP, Ir_matrix[i], marker='o',
             label=f"A={listA[i]:.4f}, B={listB[i]:.2f}")
plt.xlabel("Puissance P (W)")
plt.ylabel("Ir (A)")
plt.title("Ir(P) pour chaque paire A,B")
plt.grid(True)
plt.legend()
plt.show()
