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
eitetha = []                                    # V ector para las Pen_refentes ssssss                        
eideflexiones = []
a_g = []                                        #Vector que almacena los limites inferiores            
b_g = []                                        #Vector que almacena los limites superiores
xa_g = []                                       #Vector que almacena las primeras posiciones    
xb_g = []                                       #Vector que almacena las sugundas posiciones
y2_g = []                                       #Vector que almacena las funciones     
yE_g = []                                       #Vector que almacena las funciones evaluadas    
elm = 1000
delt = 0.001
E = 10000000                                    # Modulo de elasticidad 
I = 10                                          # Momento de inercia
                                                

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
        Qy , area, a, b, xa, xb, y2 = Cargas_distribuidas_f(elm , LR , Vpc) # calcula el momento de area y calcula la carga equivalente
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

v_cortantes = integrar_num(v_cortantes , Vpc)

# calculo de los momentos flectores ------------
m_flectores = np.repeat(float(0),len(Vpc))              # se crea un vector para los momentos cortantes, inicialmente iguales a 0          

for t in range(len(Vpc)-1):                             # calculo de los momentos flectores por medio de integración numérica
    if(momentums[t]<0):                                 # ingreso de los momentos puntuales por medio de superposicion
        m_flectores[t] -= momentums[t]                  # se incrementa el valor del momento en dicho punto en una cantidad igual a la magnitud del par, si este último tiene un sentido a favor del movimiento de las manecillas del reloj
    else:        
        m_flectores[t] -= momentums[t]                  # se disminuye el valor del momento en dicho punto en una cantidad igual a la magnitud del par, si este último tiene un sentido en contra del movimiento de las manecillas del reloj
    
    m_flectores[t+1] = (m_flectores[t] + (v_cortantes[t+1]/elm)) # integración numérica, de divide el valor de la fuerza cortante por mil debido a la discretizacion

#Resistencia --------------------------------------------------------------------------------------------------------------------------


#Pen_refentes---------------------------------------------------------------------------------------------------------------------------
eitetha = integrar_num(eitetha , m_flectores)

if Aa == 0:
    vector_M = m_flectores[0:int((Ab*elm))]               #almacena momentos necesarios para calcular la tangente de referencia (se encuentran entre los dos apoyos)
else:
    vector_M = m_flectores[int((Aa*elm)):int((Ab*elm))]

print(m_flectores)
print(vector_M)    

xvec = np.flip(np.arange(0, Ab-Aa, delt))
primermomento = integrate.trapz(np.multiply(xvec, vector_M*(1/(E *I)), xvec))
Pen_ref = -primermomento/(Ab-Aa)
##
if Aa == 0:
    c1 = np.repeat(E * I*Pen_ref-eitetha[0],len(eitetha))
else:
    c1 = np.repeat(E * I*Pen_ref-eitetha[int((Aa*elm))-1],len(eitetha))
##
print(c1[0])
tetha = (eitetha+c1)*(1/(E * I))/elm
#hacerlo para las tres secciones 

tetha = (eitetha+c1)*(1/(E_i))/elm                                          # calculo de las pendientes utilizando la ecuacion matemática

eideflexiones = integrar_num(eideflexiones,tetha)

c2 = np.repeat(-eideflexiones[int((Ab*elm))-1],len(eideflexiones))  

deflexiones = (eideflexiones+c2)/elm


#Grafico de cargas--------------------------------------------------------------------------------------------------------------------


# Magnitud del vector_M
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

#Cargas distribuidas
a_g.append(a)                                                                #Se agregan los limites inferiores       
b_g.append(b)                                                                #Se agregan los limites superiores                   
xa_g.append(xa)                                                              #Se agregan las posiciones iniciales   
xb_g.append(xb)                                                              #Se agregan las posiciones finales   
y2_g.append(y2)                                                              #Se agregan las funciones   

for e in y2_g:
    if (type(y2_g[e]) == int)|(type(y2_g[e]) == float):                     #Si la función en una constante
        y2_g = np.repeat(y2_g[e], elm)                                      #Se repite ese valor         
        yE_g.append(y2_g)                                                   #Se agrega al vector de funciones evaluadas
    else:
        func = lambda x: eval(y2_g[e])                                      #Evalua las funciones
        yE_g.append(func)                                                   #Se agregan al vector de funciones evaluadas

for f in range(len(yE_g)):
    x = np.linspace(a_g[f], b_g[f], elm)                                    
    ax.plot(x, yE_g[f](x), label=f'Function {f+1}')
    

#Diagramas cortantes y flexionantes
fig2,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

ax1.plot(LR, v_cortantes, color = "aquamarine")
ax1.set_title("Momentos cortantes")

ax2.plot(LR, m_flectores, color = "violet")
ax2.set_title("Momentos flectores")

ax3.plot(LR, tetha, color = "red")
ax3.set_title("Pendientes")

ax4.plot(LR, deflexiones, color = "blue",)
ax4.set_title("Deflexiones")

plt.show()