import numpy as np


#Données de base
H = 20                                      #Mètre
D = 12                                      #Mètre
L = 25                                      #Mètre
d = 0.45                                    #Mètre
T = 50                                      #Degré Celsius
r = d/2                                     #Mètre
l = 240                                     #Kilomètre
V_nominal = 230                             #KiloVolts
E0 = 8.854*10**-12                          #Perméabilité du vide
f=60                                        #Fréquence (Hz)


#Calcul de L1
GMD = (D**2 * 2*D)**(1/3)
GMR0 = 0.0146
GMR0_prime = GMR0*0.7788
GMR_prime = (GMR0_prime*d**3)**(1/4)
L1_m = (2*10**-7)*np.log(GMD/GMR_prime)     #En Henry/mètre
L1 = L1_m * l*1000
print("\n------------------------Calcul de L1-------------------------------\n")
print("La valeur de L1 en H/m est:",L1_m,"H/m")
print("La valeur de L1 est:",L1,"H")
print("-------------------------------------------------------------------\n")

#Calcul de C1
C1_m = (2*np.pi*E0)/(np.log(GMD/GMR_prime)) #Farad/mètre
C1 = C1_m * l*1000
print("------------------------Calcul de C1-------------------------------\n")
print("La valeur de C1 en F/m est:",C1_m,"F/m")
print("La valeur de C1 en F est:",C1,"F")
print("-------------------------------------------------------------------\n")

#Calcul de R1
R1_m= 0.0803/(4*1.6093)                     #Ohm/mètre
R1 = R1_m * l
print("------------------------Calcul de R1-------------------------------\n")
print("La valeur de R1 en ohm/m est:",R1_m,"Ohm/m")
print("La valeur de R1 en F est:",R1,"Ohm")
print("-------------------------------------------------------------------\n")

#Calcul de Xl et de B
Xl = L1*2*np.pi*f
B = (C1*2*np.pi*f)*10
print("----------------------Calcul de Xl et B----------------------------\n")
print("La valeur de Xl en ohm est:",Xl,"Ohm")
print("La valeur de B en S est:",B,"S")
print("-------------------------------------------------------------------\n")

#Calcul de Emax
Emax = 2
print("------------------------Calcul de Emax-----------------------------\n")
print("La valeur de Emax en V est:",Emax,"V")
print("-------------------------------------------------------------------\n")