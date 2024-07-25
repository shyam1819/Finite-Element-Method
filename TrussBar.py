import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt

def removnan(array1):
    nan_array = np. isnan(array1)
    not_nan_array = ~ nan_array
    array2 = array1[not_nan_array]
    return array2

#Function for assembling stiffness matrix
def assemble(k,B,N,n):
    N=2*N
    K=np.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            for l in range(n):
                if i == j:
                    if i+1 in B[l]:
                        m =list(B[l]).index((i+1))
                        K[i][j] += k[l][m][m]
                else:
                    if i+1 in B[l] and j+1 in B[l]:
                        m = list(B[l]).index((i+1))
                        o = list(B[l]).index((j+1))
                        K[i][j] += k[l][m][o]
            K[j][i]=K[i][j]
    return K
File=input("Enter the file name eg: bar.xlsx: ")
n=int(input("Enter no of elements: "))
N=int(input("Eneter no of nodes: "))
df=pd.read_excel(File) #Dataframe
# Reading values from the excel sheet.
L=np.zeros(n)
L=df.loc[:,"Length"].values
E=np.zeros(n)
E=df.loc[:,"E"].values
A=np.zeros(n)
A=df.loc[:,"Area"].values
theta=np.zeros(n)
theta=df.loc[:,"Angle"].values
theta=theta*np.pi/180
B=np.zeros((n,4))
for i in range(4):
    e="L"
    e+=str(i+1)
    B[:,i]=df.loc[:,e].values
B=B.astype("int64")
pn=np.zeros(n)
pn=df.loc[:,"Boundary Nodes"].values
pn=pn.astype("int64")
pv=np.zeros(n)
pv=df.loc[:,"Value"].values
sn=np.zeros(n)
sn=df.loc[:,"S Variable"].values
sn=sn.astype("int64")
sv=np.zeros(n)
sv=df.loc[:,"S Value"].values
# Computing local stiffness matrix and assemble them
k=np.zeros((n,4,4))
for i in range(n):
    klocal=np.array([[1,-1],[-1,1]])
    klocal=klocal.astype("float64")
    klocal*=A[i]*E[i]/L[i]
    t=np.array([[math.cos(theta[i]),math.sin(theta[i]),0,0],[0,0,math.cos(theta[i]),math.sin(theta[i])]])
    k[i]=(np.transpose(t).dot(klocal)).dot(t)
#Assembling the stiffness matrix using asemble function
K=assemble(k,B,N,n)
#Enforcing the boundary conditions
I=np.identity(N*2)
b=K.copy()
Q=np.zeros(N*2)
j=0
for i in pn:
    if i>=0:
        b[i-1]=I[i-1]
        Q[i-1]=pv[j]
        j+=1
j=0
for i in sn:
    if i>=0:
        Q[i-1]=sv[j]
        j+=1
# 1D Bar solver
ans=input("Solving for 1D bar Y/N")
if ans=="Y":
    a=np.arange(N*2)
    a=a[a%2==1]
    K=np.delete(K,a,axis=0)
    K=np.delete(K,a,axis=1)
    Q=np.delete(Q,a)
    b=K.copy()
    I=np.identity(N)
    j=0
    b[0]=I[0]
    U=np.dot(np.linalg.inv(b),Q)
    Q=np.dot(K,U)
    print("Displacements are ",U)
    print("#######################################################################################################")
    print("Secondary variables are ",Q)
    l=[0]
    for i in range(n):
        l.append(sum(L[i] for j in range(i+1)))
    on=np.ones(n+1)
    plt.plot(l+U,on)
    plt.scatter(l+U,on)
    plt.legend("Deformed bar")
    print("Nodes displacement")

#Truss solver      
else:
    #Solving KU=Q
    U=np.zeros(N*2)
    U=np.dot(np.linalg.inv(b),Q)  # U is the primary variable i.e displacements
    Q=np.dot(K,U)   #Q is the Secondary variable
    print("Displacements are ",U)
    print("########################################################################################################")
    print("Secondary variables are ",Q)

