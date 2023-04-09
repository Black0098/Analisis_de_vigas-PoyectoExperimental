import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from funtions_1 import *
#from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon

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
coordenadas_m = []                              #Vector que almacena las coordenadas de los momentos
elm = 1000
delt = 0.001

print('Escoja el tipo de viga: \n 1. Simplemente apoya\n 2. Empotrada a la izaquierda\n 3. Empotrada a la derecha')
tipo = input()

print('Ingrese en que sistema va a trabajar')
print("1. SI    |   2. Ingles")
sys = input()

print('Ingrese la longitud de la viga:')        # hace parte del menú de usuario
L = float(input(""))                            # almacena la longitud de la viga
LR = np.arange(0 ,L ,delt)                      # crea el eje x del sistema de referencia (vector de posicion)
Lg = int(L+1)                                   # Longitud de la barra para graficar

Vpc = np.repeat(float(0),len(LR))               # Se almacenan las cargas en cada punto
momentums = np.repeat(float(0),len(LR))         # Se almacenan los momentos en cada punto 

if (tipo == '1'):
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
        
        Qys.append(Qy)                                                               # almacena el momento de area
        a_g.append(a)                                                                # almacena los limites inferiores       
        b_g.append(b)                                                                # almacena  los limites superiores                   
        xa_g.append(xa)                                                              # almacena  las posiciones iniciales   
        xb_g.append(xb)                                                              # almacena  las posiciones finales   
        y2_g.append(y2)                                                              # almacena  las funciones   
        
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
    print(df.iloc[:, 1:4])
    print('Seleccione la designación de la viga: ')
    d = int(input())
    E_i = (df.iloc[d,4])*(df.iloc[d,5])
else:
    print('Por favor ingrese un numero valido')
    E_i = 1000000

#---------------- Analisis -----------------------------------------------------------------------------------------------------------

#reaccciones----------------------------------------------------------------------------------------
xel = Xel_total(Qys,Vpc)                                #calculo del centroide la distribución de cargas

if (tipo == '1'):
    c_tot = sum(Vpc)                                        # Carga total antes de agregar las reacciones 
    R2 = (sum(Vpc)*(xel-Aa) - sum(momentums)) / (Ab-Aa)     # calculo de la reaccion en el segundo apoyo
    R1 = sum(Vpc) - R2                                      # calculo de la reaccion en el primer apoyo
    Vpc = -Vpc                          # Crea el eje y del sistema, se direccionan las cargas
    if Aa == 0:                         # si el primer apoyo se encuentra en el origen 
        Vpc[0] += R1                    # se agrega la reaccion por medio de superposicion
        Vpc[(int(Ab*elm))-1] += R2      # se agrega la reaccion por medio de superposicion
    else:                               # si el primer apoyo se encuentra a un distancia diferente de 0 del origen 
        Vpc[(int(Aa*elm))-1] += R1      # se agrega la reaccion por medio de superposicion
        Vpc[(int(Ab*elm))-1] += R2      # se agrega la reaccion por medio de superposicion
    
    #Verificacion
    print(" \n \n \n \n \n Resultados \n ")
    sys_result(sys, xel, R1, R2, c_tot)

elif (tipo == '2'):
    R = sum(Vpc)
    M_A = sum(Vpc)*(xel) - sum(momentums)
    Vpc = -Vpc
    Vpc[0] += R
    momentums[0] += M_A
    print("La reaccion en A es: {}".format(R))
    print("El momento en A es: {}".format(M_A))

elif (tipo == '3'):
    R = sum(Vpc)
    M_A = sum(Vpc)*(L-xel) - sum(momentums)
    Vpc = -Vpc
    Vpc[int(L*elm)-1] += R
    momentums[int(L*elm)-1] += M_A
    print("La reaccion en B es: {}".format(R))
    print("El momento en B es: {}".format(M_A))

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

if (tipo == '1'):
    if Aa == 0:
        vector_M = m_flectores[0:int((Ab*elm))]               #almacena momentos necesarios para calcular la tangente de referencia (se encuentran entre los dos apoyos)
    else:
        vector_M = m_flectores[int((Aa*elm)):int((Ab*elm))]
   

    xvec = np.flip(np.arange(0, Ab-Aa, delt))
    primermomento = integrate.trapz(np.multiply(xvec, vector_M*(1/(E_i)), xvec))
    Pen_ref = -primermomento/(Ab-Aa)
    ##
    if Aa == 0:
        c1 = np.repeat(E_i*Pen_ref-eitetha[0],len(eitetha))
    else:
        c1 = np.repeat(E_i*Pen_ref-eitetha[int((Aa*elm))-1],len(eitetha))
    ##
    tetha = (eitetha+c1)*(1/(E_i))/elm                                          # calculo de las pendientes utilizando la ecuacion matemática
    
    eideflexiones = integrar_num(eideflexiones,tetha)
    c2 = np.repeat(-eideflexiones[int((Ab*elm))-1],len(eideflexiones))
    deflexiones = (eideflexiones+c2)/elm

elif (tipo == '2'):
    vector_M = m_flectores[0:int(L*elm)]
    xvec = np.flip(np.arange(0, L, delt))
    Pen_ref = integrate.trapz(vector_M*(1/(E_i)), xvec)
    c1 = np.repeat(E_i*Pen_ref-eitetha[0],len(eitetha))
    tetha = (eitetha+c1)*(1/(E_i))/elm
    
    eideflexiones = integrar_num(eideflexiones,tetha)
    c2 = np.repeat(-eideflexiones[0],len(eideflexiones))  
    deflexiones = (eideflexiones+c2)/elm
    print("La deflexion en L es: {}".format(deflexiones[int(L*elm)-1]))

elif (tipo == '3'):
    vector_M = m_flectores[0:int(L*elm)]
    xvec = np.flip(np.arange(0, L, delt))
    Pen_ref = integrate.trapz(vector_M*(1/(E_i)), xvec)
    c1 = np.repeat(E_i*Pen_ref-eitetha[int(L*elm)-1],len(eitetha))
    tetha = (eitetha+c1)*(1/(E_i))/elm
    
    eideflexiones = integrar_num(eideflexiones,tetha)
    c2 = np.repeat(-eideflexiones[int(L*elm)-1],len(eideflexiones))  
    deflexiones = (eideflexiones+c2)/elm
    print("La deflexion en 0 es: {}".format(deflexiones[0]))

#Grafico de cargas--------------------------------------------------------------------------------------------------------------------

if (tipo == '1'):
    # Magnitud del vector_M
    x = 0
    y = -2

    viga = Rectangle(xy = (0,0), height = -0.3, width = Lg, edgecolor='lightslategray', facecolor='lightslategray')
    fig, ax = plt.subplots()
    ax.add_patch(viga)

    #Cargas puntuales
    for m in range(len(coordenadas_p)):
        if (coordenadas_p[m]>=0):
            ax.quiver(coordenadas_p[m], 2, x, y, scale_units='xy', scale=1, color = "g")

    #Apoyos
    ax.plot([(Aa-1), Aa, (Aa+1),(Aa-1)], [-1,-0.3,-1, -1],[(Ab-1), Ab, (Ab+1),(Ab-1)], [-1,-0.3,-1, -1], color='peru')
    ax.fill_between([Aa-1, Aa, Aa+1, Aa-1], [-1, -0.3, -1, -1], [Ab-1, Ab, Ab+1, Ab-1], facecolor='peru', alpha=0.3)


    ax.set_yticks(np.arange(0, 10, step=1))
    ax.set_xticks(np.arange(0, L+1, step=1))

    #Cargas distribuidas

    for i in range(len(y2_g)):
        if "x" not in y2_g[i]:                                             # Si la función en una constante
            exec(f"j{i} = np.repeat(int(y2_g[i]),elm)")                    # crea variables que van aumentando el numero en funcion de las F_Constantes, cada variable se repite elm_veces
            y2_g[i] = f"j{i}"                                              # Se añade la repeticion en la posicion por cada F_Constante
        
    for i in y2_g:
        func = lambda x: eval(y2_g[i])                                      #Evalua las funciones
        yE_g.append(func)                                                   #Se agregan al vector de funciones evaluadas

    for i in range(len(yE_g)):
        x = np.linspace(xa_g[i], xb_g[i], elm)                                    
        ax.plot(x, yE_g[i](x), label=f'Function {i+1}')
        ax.fill_between(x, yE_g[i](x), 0, where = yE_g[i](x)>0, interpolate = True, alpha=0.2) #Rellena la grafica

elif (tipo == '2'):
    # Magnitud del vector_M
    x = 0
    y = -2

    viga = Rectangle(xy = (0,0), height = -0.3, width = Lg, edgecolor='lightslategray', facecolor='lightslategray')
    empotrada = Rectangle(xy=(0,-(Lg/6)), height = Lg, width= -0.2, edgecolor='lightslategray', facecolor='lightslategray')
    fig, ax = plt.subplots()
    
    ax.add_patch(viga)
    ax.add_patch(empotrada)
    
    #Cargas puntuales
    for m in range(len(coordenadas_p)):
        if (coordenadas_p[m]>=0):
         ax.quiver(coordenadas_p[m], 2, x, y, scale_units='xy', scale=1, color = "g")

    #Cargas distribuidas
    for i in range(len(y2_g)):
        if "x" not in y2_g[i]:                                             # Si la función en una constante
            exec(f"j{i} = np.repeat(int(y2_g[i]),elm)")                    # crea variables que van aumentando el numero en funcion de las F_Constantes, cada variable se repite elm_veces
            y2_g[i] = f"j{i}"                                              # Se añade la repeticion en la posicion por cada F_Constante
        
    for i in y2_g:
        func = lambda x: eval(y2_g[i])                                      #Evalua las funciones
        yE_g.append(func)                                                   #Se agregan al vector de funciones evaluadas

    for i in range(len(yE_g)):
        x = np.linspace(xa_g[i], xb_g[i], elm)                                    
        ax.plot(x, yE_g[i](x), label=f'Function {i+1}')
        ax.fill_between(x, yE_g[i](x), 0, where = yE_g[i](x)>0, interpolate = True, alpha=0.2) #Rellena la grafica
  
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

