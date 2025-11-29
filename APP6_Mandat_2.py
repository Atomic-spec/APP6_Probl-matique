import numpy as np
from APP_6_Mandat_1 import*


#Calcul de Zc
Zc = (L1_m/C1_m)**0.5
print("-------------------------Calcul de Zc------------------------------\n")
print("La valeur de Zc en ohm est:",Zc,"Ohm")
print("-------------------------------------------------------------------\n")

#Calcul de Snominal
Sil = (V_nominal**2)/Zc#204e6
print("--------------------------Calcul de SIL----------------------------\n")
print("La valeur de SIL en VA est:",Sil,"VA")
print("-------------------------------------------------------------------\n")

#Calcul de Vs

print("-----------------------Calcul de Vs--------------------------------\n")
Z = 2.994 + 79.17j
Y = (2*np.pi*f*C1)*1j
A = D_mat = 1 + (Z*Y)/2
B = Z
C = Y+(Y/2)**2 *Z
MP = np.array(([A,B],[C,D]), dtype=complex)
print("Y:",Y,"\tZ",Z)
Bp = np.abs(B)

Vr = (V_nominal/np.sqrt(3))
Ir = (Sil/(np.sqrt(3)*V_nominal))

VA = Vr*A
VB = B*Ir
Vs = VA+VB
Vs_mod = np.sqrt(Vs.real**2 + Vs.imag**2)
delta = np.atan(Vs.imag/Vs_mod)
delta_degree = delta*360/(2*np.pi)
Vsll = Vs_mod*np.sqrt(3)


print("La valeur de A en pu est:",Bp)
print("La valeur de Vs en V est:",Vs,"V")
print("La valeur de Vs est de:", Vs_mod,"avec un angle de déphasage de", delta_degree)
print("La valeur de Vs ligne-ligne est de:", Vsll,"V")
print("La valeur de Is en A est:",Ir,"A")
print("-------------------------------------------------------------------\n")

#Calcul de Is
print("-----------------------Calcul de Is--------------------------------\n")
Is = C.imag*Vr + D_mat*Ir
Is_mod = np.sqrt(Is.real**2 + Is.imag**2)
deltaIs = np.atan(Is.imag/Is_mod)
deltaIs_degree = deltaIs*360/(2*np.pi)
print("La valeur de Is en A est:",Is,"A")
print("la vrai valeur de Is en A est:", Is_mod,"A avec un déphasage de", deltaIs_degree)
print("-------------------------------------------------------------------\n")

#Calcul de variation de tension

deltaTension = ((Vs_mod-Vr)/Vs_mod)*100
print("-----------------Calcul de Variation de tension--------------------\n")
print("La tension varie de :",deltaTension,"%")
print("-------------------------------------------------------------------\n")

#Calcul de Variation de tension entre VB3 et Vnominal
VB2 = 2.29*10**5
VB3 = 2.255*10**5
DeltaVB3 = 100*(VB2-VB3)/VB2
print("--------------Calcul de Variation de tension VB3-------------------\n")
print("La tension VB3 varie de :",DeltaVB3,"%")
print("-------------------------------------------------------------------\n")





