from http.cookies import BaseCookie

from APP_6_Mandat_1 import*
from APP6_Mandat_2 import *
from sympy import symbols
from sympy.solvers import solve

Rt1pu = 0.002                         #pu
Rt2pu = 0.0043218                     #pu
Lt1pu = 0.08                          #L1 en pu
Lt2pu = 0.17287                       #L2 en pu
Ltmpu = 500                           #LTm en pu
Rtmpu = 500.01                        #pu
S_nominal = 1050e6
Lm = 0.24055
FP = 1
IB3= 516                              #Ampère

#Calcul de la puissance disponible fournie par la centrale
print("-------------Calcul de puissance disponible----------------\n")
Pdispo = 0.9*(3*350e6 - 100e6)
Pligne = Pdispo/3
print("la puissance totale disponible est de", Pdispo/1e6,"MW donc",Pligne/1e6,"MW par lignes")
print("-----------------------------------------------------------\n")


#Calcul de Zbase secondaire
print("--------------Calcul de Zbase secondaire-------------------\n")
Zbase_s = ((V_nominal*1000)**2)/S_nominal
a = Vr/Vs_mod
print("la valeur de Vr est:", Vr)
print("La valeur de a est:", a)
print("La valeur de Zbase secondaire est de:",Zbase_s)
print("-----------------------------------------------------------\n")

#Calcul des résistance vu du secondaire
print("----------Calcul de R1, R2, X1, X2, Rm et Xm---------------\n")
Rt1 = Rt1pu * Zbase_s
Rt2 = Rt2pu * Zbase_s

Lt1 = Lt1pu * Zbase_s
Lt2 = Lt2pu * Zbase_s

X1 = (Lt1 *2*np.pi*f)*1j
X2 = (Lt2*2*np.pi*f)*1j

Rm = (Rtmpu * Zbase_s)
Ltm = Ltmpu * Zbase_s
Xm = (Ltm*2*np.pi*f)*1j
Rm_équi = ((1/Rm)+(1/Xm))**-1
print("La valeur de R1 =",Rt1,"et R2 =",Rt2)
print("La valeur de L1 =",Lt1, "et L2 =",Lt2)
print("La valeur de X1 =",X1, "et X2 =",X2)
print("La valeur de Rm =",Rm)
print("La valeur de Ltm =",Ltm, "et Xm =",Xm)
print("La valeur de Rm équivalente est de:", Rm_équi,"ohm")
print("-----------------------------------------------------------\n")

#Calcul Matriciel avec Circuit en T
print("--------------Calcul de matrice de réseau------------------\n")
Zt1 = Rt1 + X1
Zt2 = Rt2 + X2
Yt = Rm_équi
Am = (1+(Yt*Zt1))
Bm = (Zt1+Zt2+(Yt*Zt1*Zt2))
Cm = Yt
Dm = (1+Yt*Zt2)
MT = np.array(([Am,Bm],[Cm,Dm]), dtype=complex)

#Matrice alternateur
AC = 1
BC = 0
CC = 0
DC = 1
MA = np.array(([AC,BC],[CC,DC]), dtype=complex)


#Matrice de compensateur
kc = 0.1
kc_array = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])

#Index
k_index = np.where(kc_array == kc)[0][0]
Zc_k = kc_array[k_index]

MC = np.array(([1,Zc_k],[0,1]), dtype =complex )

# Calcul de Matrice Alternateur * Matrice T
print("début\n")
M1 = np.dot(MA,MT)
M2 = np.dot(M1,MP)
M3 = np.dot(M2,MC)
print("La matrice 1 est:\n",M1,"\n")
print("La matrice 2 est:\n",M2,"\n")
print("La matrice de réseau est:\n",M3,"\n")
print("fini")
print("-----------------------------------------------------------\n")


#Calcul de Puissance
print("----------------Calcul de Puissance------------------------\n")
P = 3*VB3*IB3*FP
print("La valeur de la puissance est de:",P/1000,"kW")
print("-----------------------------------------------------------\n")


#Calcul de Eg

aa = 230000**2
bb = M3[0][0]
cc = M3[0][1]
dd = M3[1][0]
ee = M3[1][1]
x = symbols('x',positive=True)          #sélectionne univquement les valeurs réelle positive
sol = solve(aa-((bb*x)-cc/(x))**2 - ((dd*x)+ee/(x))**2, x)

print("-----------------------Calcul de Eg------------------------\n")

print(sol)
print(0.8*Sil)
if not sol:
    print("No solution")
else:
    LaBonneSolution = sol[1]            #avec 0.8xSIL j'ai repéré que c'est la 2e valeur du tableau qui est la bonne
    print("Avec 0.8xSIL j'ai repéré que c'est la 2e valeur du tableau qui est la bonne et je fais les calculs avec la bonne solution\n")
    print("La bonne solution est:",LaBonneSolution)
    print("LaBonneSolution x2=", LaBonneSolution*2)

print("-----------------------------------------------------------\n")
