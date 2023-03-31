import numpy as np
import pandas as pd

df = pd.read_json("vigas_data.json")
print('Ingrese en que sistema va a trabajar')
print("1. SI    |   2. Ingles")
sys = input()

print('Seleccione el numero del material de la viga y tipo de perfil: ')
print('Acero:')
print('     1.W')
print('     2.S')
print('Aluminio:')
print('     3.I')
p = input()

if (p == '1'):
    if (sys == '1'):
        print(df.iloc[0:4, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    elif (sys == '2'):
        print(df.iloc[8:11, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    else:
        E_i = 100000
elif (p == '2'):#falta
    if (sys == '1'):
        print(df.iloc[0:4, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    elif (sys == '2'):
        print(df.iloc[8:11, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    else:
        E_i = 100000
elif (p == '3'):
    if (sys == '1'):
        print(df.iloc[11:14, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    elif (sys == '2'):
        print(df.iloc[4:8, 1:4])
        print('Seleccione la designación de la viga: ')
        d = int(input())
        E_i = (df.iloc[d,4])*(df.iloc[d,5])

    else:
        E_i = 100000
else:
    print('Por favor ingrese un numero valido')
    E_i = 1000000