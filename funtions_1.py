import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


def Xel_total(Qys, Vpc):                                        ##Coordenada del centriode general   
    a = sum(Qys)
    b = sum(Vpc) 
    Xel= a/b
    return Xel

def apoyos():                                                   ##Apoyos de la barra                 
    print('Ingrese la posicion del primer apoyo:')
    Aa = float(input(""))
    print('Ingrese la posicion del segundo apoyo:')
    Ab = float(input(""))
    return Aa , Ab

def get_interval():                                             ##Intervalo de la función            
    a = float(input("Limite inferior de la función: "))
    b = float(input("Limite superior de la función: "))
    return a, b

def Fun_Position():                                             ##Posición de la función en la viga  
    xa = float(input("ingrese la posición donde inicia la carga: "))
    xb = float(input("ingrese la posición donde finaliza la carga: "))
    return xa,xb

def Cargas_puntuales_f(Vpc):                                    ##Cargas y primer momento puntual    
    cordenada = float(input("ingrese la coordenada en X: "))
    load = float(input("ingrese la carga puntual: "))
    Qy = cordenada * load
    if cordenada == 0:
        Vpc[0] = Vpc[0] + load
    else:
        Vpc[int(cordenada*1000)-1] += load  
    return Qy , cordenada

def Cargas_distribuidas_f(elm,LR,Vpc , function_num = [0]):     ##Cargas y primer momento distribuido
    
    function_num[0] +=1
    
    print(f"Function {function_num[0]}:")
    
    
    y2 = input("ingrese una funcion: ")
    a, b = get_interval()                          #Limites de la funcion                   
    xa, xb= Fun_Position()                         #Posicion de la funcion

    Rectangles = int(abs(LR[int(xb*elm)-1]*elm-LR[int(xa*elm)-1]*elm))

    if xa == 0:
        Rectangles = int(abs(LR[int(xb*elm)-1]*elm-LR[int(xa*elm)]*elm))+1

    dx = (b - a) / Rectangles
    x = np.linspace(a, b, Rectangles)
    y = []

    op = {'sin': np.sin, 'cos': np.cos, 'x': x, 'exp': np.exp, 'pi': np.pi, 'sqrt': np.sqrt}
    y.append(eval(y2, op))                              #Evaluar la funcion


    if (type(y[0]) == int)|(type(y[0]) == float):
        y = np.repeat(y,Rectangles)
        Cargas_puntuales = y * dx
    else:
            Cargas_puntuales = y[0] * dx

    j = int((xa*elm))
    k = 0

    while k<len(Cargas_puntuales)-1:
        Vpc[j+k] += Cargas_puntuales[k]
        k+=1

    area = sum(Cargas_puntuales)
    distancia = np.linspace(xa, xb, Rectangles)
    Qy = integrate.trapz(np.multiply(distancia,Cargas_puntuales), distancia)*elm
    return Qy , area, xa, xb, y[0]

def Momentos_f(elm,momentums):                                  ##Vector de momentos en la viga
                                     ##Vector de momentos puntuales       s
    coordenadas_m = []                              #Vector que almacena las coordenadas de los momentos
    ms = []  
    i = 0
    print('Desea ingresar un momento?')
    print('1. Si  |  2. No')
    j = int(input())
    while (i==0):
        
        if j == 1:

            cordenada = float(input("Ingrese la coordenada en X "))
            m = float(input("Ingrese el momento "))
    
            if(cordenada==0):
                momentums[int(cordenada*elm)] += m
            else:
                momentums[int(cordenada*elm)-1] += m
            coordenadas_m.append(cordenada)
            ms.append(m)

            print('Desea agregar otro momento?')
            print('1. Si  |  2. No')
            r = input()

            if (r=='1'):
                i = 0
            else:
                i = 1
        elif j == 2:
            i = 1
    return momentums, coordenadas_m, ms

def integrar_num(V_in, V_e):                                    ##Integración numérica
    
    V_in = np.repeat(float(0),len(V_e))
    V_in[0] = V_e[0]
    for t in range(len(V_e)-1):
         V_in[t+1] = V_in[t] + V_e[t+1]
    
    return V_in

def sys_result(type, sys, xel, R1, R2, R, M_A, c_tot):          ##Resultados según el sistema
    print("\n \n \n \n \n Resultados \n")
    if type == "1":
        if sys == "1":
            print("x de elemento: {} m".format(xel))
            print("La reaccion en el primer apoyo es: {} N".format(R1))
            print("La reaccion en el segundo apoyo es: {} N".format(R2))
            print("La carga total es: {} N ". format(c_tot))
        elif sys == "2":
            print("x de elemento: {} in".format(xel))
            print("La reaccion en el primer apoyo es: {} lb".format(R1))
            print("La reaccion en el segundo apoyo es: {} lb".format(R2))
            print("La carga total es: {} lb ". format(c_tot))

    elif type == "2":
        if sys == "1":
            print("La reaccion en A es: {} N".format(R))
            print("El momento en A es: {} N·m".format(M_A))
        elif sys =="2":
            print("La reaccion en A es: {} lb".format(R))
            print("El momento en A es: {} lb·in".format(M_A))

    elif type == "3":
        if sys == "1":
            print("La reaccion en B es: {} N".format(R))
            print("El momento en B es: {} N·m".format(M_A))
        elif sys =="2":
            print("La reaccion en B es: {} lb".format(R))
            print("El momento en B es: {} lb·in".format(M_A))

def t_viga(sys,df):                                             ##Tipo de viga

    print('Seleccione el  material de la viga y tipo de perfil: \f')
    print('Acero:')
    print('     1. W')
    print('     2. S')
    print('Aluminio:')
    print('     3. I')
    p = input()

    if (p == '1'):
        if (sys == '1'):
            print(df.iloc[53:70, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        elif (sys == '2'):
            print(df.iloc[2:20, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        else:
            E_i = 100000
    elif (p == '2'):#falta
        if (sys == '1'):
            print(df.iloc[86:102, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        elif (sys == '2'):
            print(df.iloc[21:38, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        else:
            E_i = 100000
    elif (p == '3'):
        if (sys == '1'):
            print(df.iloc[71:85, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        elif (sys == '2'):
            print(df.iloc[39:52, 1:4])
            print('Seleccione la designación de la viga: ')
            d = int(input())
            E_i = (df.iloc[d,4])*(df.iloc[d,5])

        else:
            E_i = 100000
    else:
        print('Por favor ingrese un numero valido')
        E_i = 1000000
    return E_i

def vector_puntual(coordenadas_p, ax,x,y,L):                    ##Gráfica_vector puntual

    for m in range(len(coordenadas_p)):
        if (coordenadas_p[m]>=0):
            ax.quiver(coordenadas_p[m]*10/L, 2, x, y, scale_units='xy', scale=1, color = "g")
