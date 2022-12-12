#!/usr/bin/env python
# coding: utf-8

# In[17]:


#1D FEM code for Beam problems

import pandas as pd
import numpy as np
import  math
import sympy as sym
from matplotlib import pyplot as plt
from sympy.plotting import plot
def removnan(array1):
    nan_array = np. isnan(array1)
    not_nan_array = ~ nan_array
    array2 = array1[not_nan_array]
    return array2
# Function for assembling stiffness matrix
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
x=sym.Symbol("x")
File=input("Enter the file name: ")
df=pd.read_excel(File)    #Q507
n=int(input("Enter no of elements: ")) #2
N=int(input("Enter no of nodes: "))  #3
#Reading data from excel
L=np.zeros(n)
L=df.loc[:,"Length"].values
L=L.astype("float64")
EI=np.zeros(n)
EI=df.loc[:,"EI"].values
EI=EI.astype("float64")
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
K=np.zeros((2*N,2*N))
n1=np.zeros(n)
n1=df.loc[:,"Node no"].values
n1=n1.astype("int64")
ql=np.zeros(n)
ql=df.loc[:,"Uniform Load"].values
#Generation of local stiffness matrix and assembling them
for i in range(n):
    klocal=np.array([[6,-3*L[i],-6,-3*L[i]],[-3*L[i],2*L[i]**2,3*L[i],L[i]**2],[-6,3*L[i],6,3*L[i]],[-3*L[i],L[i]**2,3*L[i],2*L[i]**2]])
    klocal=klocal.astype("float64")
    klocal*=2*EI[i]/(L[i]**3)
    K[2*i:2*i+4,2*i:2*i+4]=K[2*i:2*i+4,2*i:2*i+4]+klocal
IF=np.zeros(2*N)
x1=0
#Generation internal force term
for i in range(n):
    h=L[i]
    phi_1 = 1-3*((x-x1)/h)**2 + 2*((x-x1)/h)**3
    phi_2 = -(x-x1)*(1-(x-x1)/h)**2
    phi_3 = 3*((x-x1)/h)**2 - 2*((x-x1)/h)**3
    phi_4 = -(x-x1)*(((x-x1)/h)**2 - (x-x1)/h)
    x2 = x1+L[i]
    F = np.zeros(4)
    F[0] = sym.integrate(ql[i]*phi_1,(x, x1, x2))
    F[1] = sym.integrate(ql[i]*phi_2,(x, x1, x2))
    F[2] = sym.integrate(ql[i]*phi_3,(x, x1, x2))
    F[3] = sym.integrate(ql[i]*phi_4,(x, x1, x2))
    IF[2*i:2*i+4]+=F
    x1 = x1+L[i]
I=np.identity(N*2)
b=K.copy()
Q=np.zeros(N*2)
#Enforcing boundary conditions
j=0
for i in sn:
    if i>=0:
        Q[i-1]=sv[j]
        j+=1
rhs=IF+Q
j=0
for i in pn:
    if i>=0:
        b[i-1]=I[i-1]
        rhs[i-1]=pv[j]
        j+=1
U=np.zeros(N*2)
U=np.dot(np.linalg.inv(b),rhs)  # U is the primary variable i.e displacements
rhs=np.dot(K,U)
Q=rhs-IF
print("The deflections and slopes at nodes are ",U)
print("######################################################################################################")
print("Boundary term ",Q)
print("######################################################################################################")
print("Internal force term ",IF)
#Plotting of the flexure beam
#l=[0]
#P=[]
#for i in range(n):
#    l.append(sum(L[i] for j in range(i+1)))
#for i in range(n):
#    h=L[i]
#    phi_1 = 1-3*((x-x1)/h)**2 + 2*((x-x1)/h)**3
#    #phi_2 = -(x-x1)*(1-(x-x1)/h)**2
#    phi_3 = 3*((x-x1)/h)**2 - 2*((x-x1)/h)**3
#    #phi_4 = -(x-x1)*(((x-x1)/h)**2 - (x-x1)/h)
#    x2 = x1+L[i]
#    flexture=U[2*i]*phi_1+U[2*i+2]*phi_3
#    P.append(plot(flexture,(x,l[i],l[i+1])))
#P.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




