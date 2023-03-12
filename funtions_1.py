import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

"""def centroide (w, W_eq, l):
    y = w
    A = W_eq
    Qy = smp.integrate(x*y,(x,0,l)).simplify() #Coordenada en X
    X_dis = Qy/A
    print("La coordenada en x del centroide es {}". format(X_dis))
    return X_dis

def reaciones (carga_tot, xel,Xa,momentums,Xb):
    R2 = (carga_tot*(xel-Xa)-sum(momentums)) / (Xb-Xa)
    R1 = carga_tot - R2
    print("La reaccion en el primer apoyo es: {}".format(R1))
    print("La reaccion en el segundo apoyo es: {}".format(R2))
    return R1, R2
"""
"""def distribuidas(c, area, xa, xb, rec, dx ,x):
    #dx = (xb-xa)/rec
    #Qy = sum(c)*(xb-xa)
    #Qy = sum(c*dx)
    Qy = integrate.trapz(np.multiply(x,c), x)*1000
    return area , Qy"""

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

def get_interval():                                             ##Intervalo de la funcion            
    a = float(input("Limite inferior de la función: "))
    b = float(input("Limite superior de la función: "))
    return a, b

def Fun_Position():                                             ##Posicion de la funcion en la viga  
    xa = float(input("ingrese la posición donde inicia la carga"))
    xb = float(input("ingrese la posición donde finaliza la carga"))
    return xa,xb

def Cargas_puntuales_f(Vpc):                                    ##cargas y primer momento puntual    
    cordenada = float(input("ingrese la coordenada en X: "))
    load = float(input("ingrese la carga puntual: "))
    Qy = cordenada * load
    if cordenada == 0:
        Vpc[0] = Vpc[0] + load
    else:
        Vpc[int(cordenada*1000)-1] += load  
    return Qy , cordenada

def Cargas_distribuidas_f(elm,LR,Vpc , function_num = [0]):     ##cargas y primer momento distribuido
    
    function_num[0] +=1
    
    print(f"Function {function_num[0]}:")
    
    
    y2 = input("ingrese una funcion: ")
    a, b = get_interval()                          #Limites de la funcion                   
    xa, xb= Fun_Position()                         #Posicion de la funcion

    Rectangles = int(abs(LR[int(xb*elm)-1]*elm-LR[int(xa*elm)-1]*elm))

    if xa == 0:
        Rectangles = int(abs(LR[int(xb*elm)-1]*elm-LR[int(xa*elm)]*elm))

    dx = (b - a) / Rectangles
    x = np.linspace(a, b - dx, Rectangles+1)
    y = []
    y.append(eval(y2))                              #Evaluar la funcion


    if (type(y[0]) == int)|(type(y[0]) == float):
        y = np.repeat(y,Rectangles+1)
        Cargas_puntuales = y * dx
    else:
            Cargas_puntuales = y[0] * dx

    j = int((xa*elm))
    k = 0

    while k<len(Cargas_puntuales)-1:
        Vpc[j+k] = Vpc[j+k] + Cargas_puntuales[k]
        k+=1

    area = sum(Cargas_puntuales)
    Qy = integrate.trapz(np.multiply(x,Cargas_puntuales), x)*elm
    return Qy , area

def Momentos_f(elm,momentums):                                  ##Vector de momentos puntuales       s
    
    i = 0
    print('Desea ingresar un momento?')
    print('1. Si  |  2. No')
    j = int(input())
    while (i==0):
        
        if j == 1:
            m = float(input("ingrese el momento "))
            cordenada = float(input("ingrese la coordenada en X "))
            momentums[int(cordenada*elm)-1] += m

            print('Desea agregar otro momento?')
            print('1. Si  |  2. No')
            r = input()

            if (r=='1'):
                i = 0
            else:
                i = 1
        elif j == 2:
            i = 1
    return momentums

def integrar_num(V_in, V_e, tetha):
    if len(tetha) == 0:
        V_in = np.repeat(float(0),len(V_e))
        V_in[0] = V_e[0]
        for t in range(len(V_e)-1):
            V_in[t+1] = V_in[t] + V_e[t+1]
    else:
        V_in = np.repeat(float(0),len(V_e))
        for t in range(len(V_e)-1):
            V_in[t+1] = V_in[t] + tetha[t+1]


    return V_in