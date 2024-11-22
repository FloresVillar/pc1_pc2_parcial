import numpy as np     
import math as mat 
from LDL_mi_version import isSymetric 
from LU_mi_version import sustitucion_regresiva
from LU_mi_version import sustitucion_progresiva

def cholesky(a,g):
    n = len(a)
    for j in range(n):
        suma = 0
        for k in range(j):
            suma = suma + g[j,k]**2
        g[j,j] = mat.pow(a[j,j]-suma,0.5)
        for i in range(j+1,n):
            suma = 0
            for k in range(j):
                suma = suma + g[i,k]*g[j,k]
            g[i,j] = (a[i,j]-suma)/g[j,j]

def l_inicial(g):
    n = len(g)
    for j in range(n):
        for i in range(n):
            if(i==j):
                g[i,j]=1.0

if __name__=='__main__':
    a = np.array([[12.0,1.0,4.0,4.0],[6.0,10.0,15.0,18.0],[4.0,5.0,8.0,7.0],[4.0,5.0,7.0,1.0]])
    b = np.array([[0.0],[20.0],[9.0],[50.0]])
    n = len(a)
    if(not isSymetric(a)):
        b = np.transpose(a)@b
        a = np.transpose(a)@a
    g = np.zeros((n,n),float)
    l_inicial(g)
    cholesky(a,g)
    z = sustitucion_progresiva(g,b) #gz= b' es inferior- progresiva
    x = np.zeros((n,1),float)       #g^t x = z   superior - regresiva    
    x = sustitucion_regresiva(np.transpose(g),z)
    print(x)
