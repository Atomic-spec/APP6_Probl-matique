from APP_6_Mandat_1 import Xl, R1
from APP6_Mandat3 import Rt1, Rt2, X1, X2
import numpy as np


print(f"\n-------------------------------MANDAT 5----------------------------------")
Zbase = (230000)**2/(1050e6)
Ibase = 1050e6 / (np.sqrt(3) * 230e3)
Xd = 0.15

Xsys_pu = (X1 + X2 + Xd + Xl) / Zbase
Rsys_pu = (Rt1 + Rt2 + R1) / Zbase

# Impédance correcte
Zsys_pu = Rsys_pu + 1j * Xsys_pu
print(f"Zsys = {abs(Zsys_pu*Zbase):.6f} Ohms")
print("RTotal", Rsys_pu*Zbase)
print("Xtotal", Xsys_pu*Zbase)

# Courant symétrique
Isym_pu = 1 / abs(Zsys_pu)
print(f"Ibase = {Ibase:.6f}")
# Constante de temps (en p.u.)
tau = Xsys_pu / Rsys_pu

# Temps d'ouverture
t_op = 3 / 60

# facteur d'asymétrie
k = np.sqrt(1 + 2 * np.exp(-4*np.pi * t_op / tau))

Iasym_pu = Isym_pu * k

# Conversion dans le système réel
Isym = Isym_pu * Ibase
Iasym = Iasym_pu * Ibase

print(f"Isym = {Isym:.6f} A")
print(f"Iasym = {Iasym:.6f} A")
print("R1", R1, "RT1", Rt1, "Rt2", Rt2, "X1", X1, "X2", X2)
print("--------------------------------------------------------------------------")