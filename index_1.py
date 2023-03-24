import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from funtions_1 import *

#----------declaración de variables iniciales ------------------------------------------------------------------------------------
df = pd.read_json("vigas_data.json")            # Base de datos de tipos de vigas
Qys = []                                        # Vector para los primeros momentos de area 
coordenadas_p = []                              # Vector para las coordenadas de cargas puntuales
v_cortantes = []                                # Vector para la los cortantes 
eitetha = []                                    # Vector para la integracion numérica para las pendientes (E*I*Tetha)
tetha = []                                      # Vector para las pendientes
eideflexiones = []                              # Vector para la integracion numérica para las deflexiones (E*I*Tetha)

elm = 1000                                      # numero de elemnetos diferenciales 
delt = 0.001                                    # paso del método numerico

print('Ingrese la longitud de la viga:')        # hace parte del menú de usuario
L = float(input(""))                            # almacena la longitud de la viga
LR = np.arange(0 ,L ,delt)                      # crea el eje x del sistema de referencia (vector de posicion)
Lg = int(L+1)                                   # Longitud de la barra para graficar

Vpc = np.repeat(float(0),len(LR))               # Se almacenan las cargas en cada punto
momentums = np.repeat(float(0),len(LR))         # Se almacenan los momentos en cada punto 

Aa, Ab = apoyos()                               # Se ingresan los apoyos

# ----------- Menu de ususario ---------------------------------------------------------------------------------------------------------
i = 0                                               # iterador
while (i==0):                                       # ciclo de ingreso de carga
    print('Escoja la carga que quiere aplicar')     # hace parte del menu de usuario
    print('1. Carga puntual')                       # "             "
    print('2. Carga distribuida')                   # "             "
    n = input()                                     # variable de switch 

#Cargas puntuales ------------------------------------------------------------------------------------------------------
    if n == '1':
        Qy , cordenada = Cargas_puntuales_f(Vpc)    # calcula el momento de area e ingresa la cordenada de la carga 
        Qys.append(Qy)                              # almacena el momento de area
        coordenadas_p.append(cordenada)             # almacena la carga

#Cargas distribuidas ---------------------------------------------------------------------------------------------------
    elif n == '2':
        Qy , area = Cargas_distribuidas_f(elm , LR , Vpc) # calcula el momento de area y calcula la carga equivalente
        Qys.append(Qy)                                    # almacena el momento de area
# ----
    else:                                               # hace parte del menu de usuario
        print('Por favor ingrese un numero valido')     # "             "

    print('Desea agregar otra carga?')                  # "             "
    print('1. Si  |  2. No')                            # "             "
    r = input()                                         # variable de switch 
    if (r=='1'):                                        # se sale del ciclo de ingreso de carga
        i = 0                                           # "             "
    else:                                               # repite el ciclo de ingreso de carga
        i = 1                                           # "             "

momentums = Momentos_f(elm , momentums)

print('Seleccione el numero del material de la viga: ')
print('1. Acero')
print('2. Aluminio')
m = input()

print('Seleccione el numero del tipo de perfil de la viga: ')
print('1. W')
p = input()

if (p == '1'):
    print(df.iloc[:, 0:3])
    print('Seleccione la designación de la viga: ')
    print(df.iloc[1:, 0:3])
    d = int(input())
    E_i = (df.iloc[d,2])*(df.iloc[d,3])
else:
    print('Por favor ingrese un numero valido')

#---------------- Analisis -----------------------------------------------------------------------------------------------------------

#reaccciones----------------------------------------------------------------------------------------
xel = Xel_total(Qys,Vpc)                                #calculo del centroide la distribución de cargas

R2 = (sum(Vpc)*(xel-Aa) - sum(momentums)) / (Ab-Aa)     # calculo de la reaccion en el segundo apoyo
R1 = sum(Vpc) - R2                                      # calculo de la reaccion en el primer apoyo

#Verificacion
print("x de elemento: {}".format(xel))
print("La reaccion en el primer apoyo es: {}".format(R1))
print("La reaccion en el segundo apoyo es: {}".format(R2))
print("La carga total es: ", sum(Vpc))
#

Vpc = -Vpc                          # Crea el eje y del sistema, se direccionan las cargas
if Aa == 0:                         # si el primer apoyo se encuentra en el origen 
    Vpc[0] += R1                    # se agrega la reaccion por medio de superposicion
    Vpc[(int(Ab*elm))-1] += R2      # se agrega la reaccion por medio de superposicion
else:                               # si el primer apoyo se encuentra a un distancia diferente de 0 del origen 
    Vpc[(int(Aa*elm))-1] += R1      # se agrega la reaccion por medio de superposicion
    Vpc[(int(Ab*elm))-1] += R2      # se agrega la reaccion por medio de superposicion

#Estatica: fuerzas cortantes y momentos flectores ----------------------------------------------------------------------------------------------------------------

v_cortantes = integrar_num(v_cortantes , Vpc,tetha)     # calculo de las fuerzas cortantes

# calculo de los momentos flectores ------------
m_flectores = np.repeat(float(0),len(Vpc))              # se crea un vector para los momentos cortantes, inicialmente iguales a 0          

for t in range(len(Vpc)-1):                             # calculo de los momentos flectores por medio de integración numérica
    if(momentums[t]<0):                                 # ingreso de los momentos puntuales por medio de superposicion
        m_flectores[t] -= momentums[t]                  # se incrementa el valor del momento en dicho punto en una cantidad igual a la magnitud del par, si este último tiene un sentido a favor del movimiento de las manecillas del reloj
    else:        
        m_flectores[t] -= momentums[t]                  # se disminuye el valor del momento en dicho punto en una cantidad igual a la magnitud del par, si este último tiene un sentido en contra del movimiento de las manecillas del reloj
    
    m_flectores[t+1] = (m_flectores[t] + (v_cortantes[t+1]/elm)) # integración numérica, de divide el valor de la fuerza cortante por mil debido a la discretizacion

#Resistencia --------------------------------------------------------------------------------------------------------------------------

#----- pendientes -------------------------------------------------------------
eitetha = integrar_num(eitetha , m_flectores, tetha)                        # Integración numérica de los momentos flectores sin resolver el problema de valor inicial

#solución del P.V.I teniendo en cuenta las condiciones de frontera y por medio de los teoremas de los momentos de area, se utiliza la región comprendida entre los dos apoyos -----

if Aa == 0:                                                                 # seleccion de los momentos flectores comprendidos entre los dos apoyos
    vector = m_flectores[0:int((Ab*elm))]
else:
    vector = m_flectores[int((Aa*elm)):int((Ab*elm))]

# se calcula la desviacion tangencial del segundo apoyo con respecto al primer apoyo, aplicando el teorema del segundo momento de área con respecto a un eje vertical que pasa por el segundo apoyo
xvec = np.flip(np.arange(0, Ab-Aa,delt))                                    # Definicion del eje de cordenadas con respecto al segundo apoyo
des_tan = integrate.trapz(np.multiply(xvec, vector*(1/(E_i)), xvec))        # calculo de la desviacion tangencial por medio del segundo momento del area
theta_A = -des_tan/(Ab-Aa)                                                  # pendiente en el primer apoyo, cociente entre la desviacion tangencial del segundo apoyo con respecto al segundo, y la destancia entre estos (se utiliza la definicion de pendiente)
##

#calculo de la constante c de integración utilizando las condiciones de frontera
if Aa == 0: 
    c1 = np.repeat(E_i*theta_A-eitetha[0],len(eitetha))
else:
    c1 = np.repeat(E_i*theta_A-eitetha[int((Aa*elm))-1],len(eitetha))
##

tetha = (eitetha+c1)*(1/(E_i))/elm                                          # calculo de las pendientes utilizando la ecuacion matemática

#deflexiones-------------------------------------------------------------

eideflexiones = integrar_num(eideflexiones, Vpc, tetha)                     # Integración numérica de las pendientes sin resolver el problema de valor inicial

c2 = np.repeat(-eideflexiones[int((Ab*elm))-1],len(eideflexiones))  

deflexiones = (eideflexiones+c2)/elm


#Grafico de cargas--------------------------------------------------------------------------------------------------------------------
print("prueba_Sync_3." , sum(Vpc))

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