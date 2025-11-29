import numpy as np
from APP_6_Mandat_1 import Xl
from APP6_Mandat_2 import A
Vrll = 230e3
Vsll = 230e3
f = 60.0


Aabs = np.abs(A)

#Calcul de l'angle
Pr = 10e6
Rep_angle = -np.arccos((Pr*Xl)/(Vrll*Vsll)) + (np.pi/2)
print("Rep_angle (rad) =", Rep_angle)
print("Rep_angle (deg) =", np.degrees(Rep_angle))

#Calcul Q
P = (Vrll * Vsll) / Xl * np.cos(Rep_angle)
P2 = (Aabs * Vrll**2) / Xl
Q_total = P - P2
H = ""
print("Q_total (vars) = ", Q_total)

