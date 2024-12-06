import numpy as np

if __name__ == "__main__":
    L11 =np.array([[22802.67337, 0.00000, 0.00000 ,0.00000],
[17871.10999, 10086.02968, 0.00000 ,0.00000],
[13851.43965 ,8399.94470, 6795.69980 ,0.00000],
[16205.15437, 2135.01879, 6051.63983, 13238.08015]])
    A10=np.array([[368509378.00000 ,334329452.00000 ,262652770.00000 ,380720332.00000],
[397114485.00000, 350432462.00000 ,316364739.00000 ,375965118.00000], 
[466365863.00000, 426471740.00000 ,345866718.00000 ,432874355.00000],
[371443738.00000 ,371697971.00000 ,337512163.00000 ,386924395.00000]])
    L11t=np.transpose(L11)
    L11ti=np.linalg.inv(L11t)
    print(L11)
    print(L11t)
    print(L11ti)
    print(A10@L11ti)
    A = np.array([[1,2,4],[3,2,4],[6,7,8]])
    U = np.array([[3,5,4],[0,4,5],[0,0,9]])
    print(A)
    print(U)
    print(A@U)
    M =np.array([[-6395205761.41770800 ,-8066020962.38804100, -5984246602.65509000, -5401668591.06543600],
[-8066020962.38804100, -10055691338.84700800 ,-7453983359.02241200, -6761615613.37243500]
,[-5984246602.65509000 ,-7453983359.02241200 ,-5568834084.26806200, -5018381694.00170500]
,[-5401668591.06543600, -6761615613.37243500 ,-5018381694.00170500 ,-4544204069.29161600]
])
    print(M)
    l,v =np.linalg.eig(M)
    print(l)
    print(v)
"""//Imprimir("ABloques[1][0]",ABloques[1][0]);
        //Imprimir("ABloques[2][0]",ABloques[2][0]);
        //Imprimir("Transpuesta de L11",Transponer(L11));
       // Imprimir("L11 L11^t= ABloques[0][0] prodParalelo", ProdParalelo(L11, Transponer(L11)));
       // Imprimir("L11 L11^t= ABloques[0][0] prodTriangular", ProductoTriangular(L11, Transponer(L11)));
        //Imprimir("ABloques[0][0]", ABloques[0][0]);
        
        //Imprimir("L11^t L11^-t",ProdParalelo(Transponer(L11),InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("L11^t L11^-t",ProductoTriangular(Transponer(L11),InversaTriangularSuperior(Transponer(L11))));
       // Imprimir("L11^t",Transponer(L11));
        //Imprimir("L11^-t",InversaTriangularSuperior(Transponer(L11)));
        Imprimir("ABloques[2][1]", ABloques[1][0]);
        //Imprimir("ABloques[1][0] L11^-t productoTriangular",ProductoTriangular(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("ABloques[1][0] L11^-t productoParalelo",ProdParalelo(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));
        //Imprimir("ABloques[1][0] L11^-t productoSerial",ProdSerial(ABloques[1][0], InversaTriangularSuperior(Transponer(L11))));"""