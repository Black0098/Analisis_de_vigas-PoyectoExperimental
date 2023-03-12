import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from funtions_1 import *
                                 
Qys = []                                        # Vector para los primeros momentos de area 
coordenadas_p = []                              # Vector para las coordenadas de cargas puntuales
v_cortantes = []                                # Vector para la los cortantes 
eitetha = []                                    #Vector para las pendientes ss
tetha = []                          
eideflexiones = []
elm = 1000
delt = 0.001

print('Ingrese la longitud de la viga:')
L = float(input(""))
LR = np.arange(0 ,L ,delt)
Lg = int(L+1)                                   #Longitud de la barra para graficar


Vpc = np.repeat(float(0),len(LR))               #Se almacenan las cargas en cada punto
momentums = np.repeat(float(0),len(LR))         #Se almacenan los momentos en cada punto 

Aa, Ab = apoyos()                               #Se ingresan los apoyos



#Menu de ususario
i = 0
while (i==0):
    print('Escoja la carga que quiere aplicar')
    print('1. Carga puntual')
    print('2. Carga distribuida')
    n = int(input())

#Cargas puntuales y distribuidas------------------------------------------------------------------------------------------------------
    if n == 1:
       
        Qy , cordenada = Cargas_puntuales_f(Vpc)
        Qys.append(Qy)
        coordenadas_p.append(cordenada)

    elif n == 2:
       
        Qy , area = Cargas_distribuidas_f(elm , LR , Vpc)
        Qys.append(Qy)

    print('Desea agregar otra carga?')
    print('1. Si  |  2. No')
    r = input()
    if (r=='1'):
        i = 0
    else:
        i = 1

momentums = Momentos_f(elm , momentums)

#reaccciones--------------------------------------------------------------------------------------------------------------------------
xel = Xel_total(Qys,Vpc)

R2 = (sum(Vpc)*(xel-Aa) - sum(momentums)) / (Ab-Aa)
R1 = sum(Vpc) - R2

print("x de elemento: {}".format(xel))
print("La reaccion en el primer apoyo es: {}".format(R1))
print("La reaccion en el segundo apoyo es: {}".format(R2))
print("La carga total es: ", sum(Vpc))

# Las reacciones entran en Vpc-------------------------------------------------------------------------------------------------------
Vpc = -Vpc
if Aa == 0:
    Vpc[0] += R1
    Vpc[(int(Ab*elm))-1] += R2
else:
    Vpc[(int(Aa*elm))-1] += R1
    Vpc[(int(Ab*elm))-1] += R2

#cortantes y flectores----------------------------------------------------------------------------------------------------------------

v_cortantes = integrar_num(v_cortantes , Vpc,tetha)


m_flectores = np.repeat(float(0),len(Vpc))

for t in range(len(Vpc)-1):
    if(momentums[t]<0):
        m_flectores[t] -= momentums[t] #incrementando el valor de M en dicho punto en una cantidad igual a la magnitud del par, si este último tiene un sentido a favor del novimiento de las manecillas del reloj
    else:        
        m_flectores[t] -= momentums[t] 
    m_flectores[t+1] = (m_flectores[t] + (v_cortantes[t+1]/elm))

#resistencia--------------------------------------------------------------------------------------------------------------------------

df = pd.read_json("vigas_data.json")
print(df)
print(df.iloc[:, 0:3])
print(df.iloc[1,2])

#pendientes---------------------------------------------------------------------------------------------------------------------------
eitetha = integrar_num(eitetha , m_flectores, tetha)

if Aa == 0:
    vector = m_flectores[0:int((Ab*elm))]
else:
    vector = m_flectores[int((Aa*elm)):int((Ab*elm))]

print(m_flectores)
print(vector)    

xvec = np.flip(np.arange(0, Ab-Aa,delt))
primermomento = integrate.trapz(np.multiply(xvec, vector*(1/(2000*2070000000000)), xvec))
pendi = -primermomento/(Ab-Aa)
##
if Aa == 0:
    c1 = np.repeat(2000*2070000000000*pendi-eitetha[0],len(eitetha))
else:
    c1 = np.repeat(2000*2070000000000*pendi-eitetha[int((Aa*elm))-1],len(eitetha))
##
print(c1[0])
tetha = (eitetha+c1)*(1/(2000*2070000000000))
#hacerlo para las tres secciones 

#deflexiones--------------------------------------------------------------------------------------------------------------------------

eideflexiones = integrar_num(eideflexiones, Vpc, tetha)

c2 = np.repeat(-eideflexiones[int((Ab*elm))-1],len(eideflexiones))  

deflexiones = (eideflexiones+c2)



print(c2[0])
print(deflexiones)

#Grafico de cargas--------------------------------------------------------------------------------------------------------------------


# Magnitud del vector
x = 0
y = -2

fig, ax = plt.subplots()
ax.plot([1]*Lg, color = "k", linewidth = 3)

#Cargas puntuales
for m in range(len(coordenadas_p)):
    if (coordenadas_p[m]>=0):
        ax.quiver(coordenadas_p[m], 3, x, y, scale_units='xy', scale=1, color = "g")
#Apoyos
ax.plot([(Aa-1), Aa, (Aa+1),(Aa-1)], [0,1,0, 0],[(Ab-1), Ab, (Ab+1),(Ab-1)], [0,1,0,0], color = "k")
ax.set_yticks(np.arange(0, 10, step=1))
ax.set_xticks(np.arange(0, L+1, step=1))

#Diagramas cortantes y flexionantes
fig2,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

ax1.plot(LR, v_cortantes, color = "aquamarine")
ax1.set_title("Momentos cortantes")

ax2.plot(LR, m_flectores, color = "violet")
ax2.set_title("Momentos flectores")

ax3.plot(LR, tetha, color = "red")
ax3.set_title("Momentos flectores")

ax4.plot(LR, deflexiones, color = "blue")
ax4.set_title("Momentos flectores")

plt.show()

