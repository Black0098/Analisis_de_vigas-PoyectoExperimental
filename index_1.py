import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from funtions_1 import *
from matplotlib.patches import Rectangle, FancyArrowPatch

#----------declaración de variables iniciales ------------------------------------------------------------------------------------
df = pd.read_json("vigas_data.json")            # Base de datos de tipos de vigas
Qys = []                                        # Vector para los primeros momentos de area 
coordenadas_p = []                              # Vector para las coordenadas de las cargas puntuales
coordenadas_m = []                              # Vector que almacena las coordenadas de los momentos puntuales
v_cortantes = []                                # Vector para las fuerzas cortantes 
eitetha = []                                    # Vector para la integracion numerica del producto E*I*Theta                       
eideflexiones = []                              # Vector para la integracion numerica del producto E*I*Y

maximos = []                                    # Maximos valores de las funciones ingresadas
xa_g = []                                       # Vector que almacena las primeras posiciones    
xb_g = []                                       # Vector que almacena las sugundas posiciones
y_axiss = []                                    # Vector que almacena todos los valores de las cargas distribuidas

elm = 1000                                      # Discretizacion del sistema, elementos por unidad
delt = 0.001                                    # Espacio entre posicion

print('Escoja el tipo de viga: \n1. Simplemente apoyada \n2. Empotrada a la izaquierda \n3. Empotrada a la derecha')
tipo = input()

print('Ingrese en que sistema va a trabajar')
print("1. SI    |   2. Ingles")
sys = input()

print('Ingrese la longitud de la viga:')
L = float(input(""))                            # almacena la longitud de la viga
LR = np.arange(0 ,L ,delt)                      # crea el eje x del sistema de referencia (vector de posicion)

Vpc = np.repeat(float(0),len(LR))               # Creacion del vector en el que se almacenan las cargas en cada punto
momentums = np.repeat(float(0),len(LR))         # Creacion del vector en el que se almacenan los momentos en cada punto 

if (tipo == '1'):                                   # Si la viga es simplemente apoyada
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
        Qy , area, xa, xb, y_axis = Cargas_distribuidas_f(elm , LR , Vpc) # calcula el momento de area y calcula la carga equivalente
        
        Qys.append(Qy)                                                               # almacena el momento de area                   
        xa_g.append(xa)                                                              # almacena  las posiciones iniciales   
        xb_g.append(xb)                                                              # almacena  las posiciones finales      
        y_axiss.append(y_axis)                                                       # almacena todos los valores de la carga distribuida

    else:                                               # hace parte del menu de usuario
        print('Por favor ingrese un numero valido')     # "             "

    print('Desea agregar otra carga?')                  # "             "
    print('1. Si  |  2. No')                            # "             "
    r = input()                                         # variable de switch 
    if (r=='1'):                                        # se sale del ciclo de ingreso de carga
        i = 0                                           # "             "
    else:                                               # repite el ciclo de ingreso de carga
        i = 1                                           # "             "

momentums, coordenadas_m, ms = Momentos_f(elm , momentums)
E_i = t_viga(sys,df)
    

#---------------- Analisis -----------------------------------------------------------------------------------------------------------

#reaccciones----------------------------------------------------------------------------------------
xel = Xel_total(Qys,Vpc)                                #calculo del centroide la distribución de cargas

if (tipo == '1'):
    c_tot = sum(Vpc)
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
    sys_result(tipo, sys, xel, R1, R2, R=0, M_A=0, c_tot=c_tot)
    
elif (tipo == '2'):
    R = sum(Vpc)
    M_A = sum(Vpc)*(xel) - sum(momentums) - momentums[0]
    Vpc = -Vpc
    Vpc[0] += R
    momentums[0] += M_A
    sys_result(tipo, sys, xel, R1=0, R2=0, R=R, M_A=M_A, c_tot=0)
    
elif (tipo == '3'):
    R = sum(Vpc)
    M_A = sum(Vpc)*(xel-L) - sum(momentums) - momentums[int(L*elm)-1]
    print('momento prueba')
    print(M_A)
    Vpc = -Vpc
    Vpc[int(L*elm)-1] += R
    momentums[int(L*elm)-1] += M_A
    sys_result(tipo, sys, xel, R1=0, R2=0, R=R, M_A=M_A, c_tot=0)


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
viga = Rectangle(xy = (0,0), height = -0.3, width = 10, edgecolor='lightslategray', facecolor='lightslategray')

fig, ax = plt.subplots()

# Magnitud del vector_M
x = 0
y = -2

ax.set_title("Análisis Estructural")
ax.add_patch(viga)

if (tipo == '1'):                                       #graficas viga apoyada
    
    #Cargas puntuales
    vector_puntual(coordenadas_p, ax,x,y, L)
           
    #Apoyos
    ax.plot([(Aa*10/L-1), Aa*10/L, (Aa*10/L+1),(Aa*10/L-1)], [-2,-0.3,-2, -2],[(Ab*10/L-1), Ab*10/L, (Ab*10/L+1),(Ab*10/L-1)], [-2,-0.3,-2, -2], color='peru')
    ax.fill_between([(Aa*10/L-1), Aa*10/L, (Aa*10/L+1),(Aa*10/L-1)], [-2,-0.3,-2, -2], facecolor='peru', alpha=0.3)
    ax.fill_between([(Ab*10/L-1), Ab*10/L, (Ab*10/L+1),(Ab*10/L-1)], [-2,-0.3,-2, -2], facecolor='peru', alpha=0.3)

elif (tipo == '2'):                                     #graficas viga empotrada izq

    empotrada = Rectangle(xy=(0,-(10/6)), height = 6, width= -0.3, edgecolor='lightslategray', facecolor='lightslategray')
    ax.add_patch(empotrada)
    
    #Cargas puntuales
    vector_puntual(coordenadas_p, ax,x,y,L)

elif (tipo == '3'):                                     #graficas viga apoyada der

    empotrada = Rectangle(xy=(10,-(10/6)), height = 6, width= 0.3, edgecolor='lightslategray', facecolor='lightslategray')
    ax.add_patch(empotrada)
    
    #Cargas puntuales
    vector_puntual(coordenadas_p, ax,x,y,L)

for i in range(len(y_axiss)):
        if ((type(y_axiss[i]) == np.int32)|(type(y_axiss[i]) == np.float32)):
            print(y_axiss[i])
            y_axiss[i] = np.repeat(y_axiss[i],int(abs(xb_g[i]*elm)-xa_g[i]*elm))
            print(y_axiss[i])
        maximos.append(max(y_axiss[i]))
if (len(y_axiss)!=0): 
    maximo = max(maximos)
    print('maximo')
    print(maximo)

for i in range(len(y_axiss)):
    rec = int(abs(LR[int(xb_g[i]*elm)-1]*elm-LR[int(xa_g[i]*elm)-1]*elm))
    if xa_g[i] == 0:
        rec = int(abs(LR[int(xb_g[i]*elm)-1]*elm-LR[int(xa_g[i]*elm)]*elm))+1

    xx = np.linspace(xa_g[i], xb_g[i], rec)
    ax.plot(xx*10/L, y_axiss[i]*8/maximo, label=f'Function {i+1}') 
    ax.fill_between(xx*10/L, y_axiss[i]*8/maximo, 0, where = y_axiss[i]*8/maximo>0, interpolate = True, alpha=0.2) #Rellena la grafica"""
ax.set_axis_off()

#Grafico de momentos
for k in  range(len(coordenadas_m)):
    # Definir los puntos de la flecha
    if (ms[k]>1):
        x0, y0 = coordenadas_m[k]-1.15, -1.2
        x1, y1 =  coordenadas_m[k] - 1.15, 1.2
    elif (ms[k]<1):
        x0, y0 = coordenadas_m[k]+1.15, -1.2
        x1, y1 =  coordenadas_m[k] + 1.15, 1.2
    
    # Dibujar la flecha recta
    ax.arrow(x0, y0, x1, y1, head_width=0.1, head_length=0.2, fc='none', ec='none')
    
    if (ms[k]>1):
        curvature = 1 # Ajusta el valor para cambiar la curvatura del arco

    elif(ms[k]<1):
        curvature = -1

    # Dibujar la flecha curvada
    arrow = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='->,head_width=0.2,head_length=0.4',
                            mutation_scale= 10, lw=3, color='darkorchid', alpha=0.7,
                            connectionstyle=f'arc3,rad={curvature}', zorder=10)
    ax.add_patch(arrow)
  
#Diagramas cortantes y flexionantes ---------------------------------------------------------------------------------------------------
fig2,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

ax1.plot(LR, v_cortantes, color = "aquamarine")
ax1.set_title("Fuerzas cortantes")
ax1.grid()

ax2.plot(LR, m_flectores, color = "violet")
ax2.set_title("Momentos flectores")
ax2.grid()

ax3.plot(LR, tetha, color = "red")
ax3.set_title("Pendientes")
ax3.grid()

ax4.plot(LR, deflexiones, color = "blue",)
ax4.set_title("Deflexiones")
ax4.grid()

if(sys=='1'):
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('V(x) (N)')
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('M(x) (N*m)')
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel('Theta(x) (rad)')
    ax4.set_xlabel('x (m)')
    ax4.set_ylabel('Y(x) (m)')
elif(sys=='2'):
    ax1.set_xlabel('x (in)')
    ax1.set_ylabel('V(x) (lb)')
    ax2.set_xlabel('x (in)')
    ax2.set_ylabel('M(x) (lb*in)')
    ax3.set_xlabel('x (in)')
    ax3.set_ylabel('Theta(x) (rad)')
    ax4.set_xlabel('x (in)')
    ax4.set_ylabel('Y(x) (in)')

plt.tight_layout()
plt.show()
