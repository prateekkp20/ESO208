from copy import copy, deepcopy
import numpy
from math import *
n=int(input('Enter the number of equations:'))
print('Enter the matrix:')
a=[list(map(float, input().split())) for x in range(n)]
q=input('Please select your method:\na.GE without pivoting\nb.GE with pivoting\nc.GE with scaling and pivoting\nd.LU decomposition by using GE without pivoting\ne.LU decomposition by using GE with pivoting\nf.LU decomposition by using Crout method without pivoting\ng.Cholesky decomposition for symmetric positive definite matrix\n')

def GE_without_pivoting(mat,n):
    #Forward elimination
    for i in range(0,n-1):
        for j in range(i+1,n):
            alpha=-(mat[j][i])/mat[i][i]
            mat[j][i]=0
            for k in range(i+1,n+1):
                mat[j][k]=mat[j][k]+alpha*mat[i][k]

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=mat[n-1][n]/mat[n-1][n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+mat[i][j]*x[j]
        x[i]=(mat[i][n]-sumi)/mat[i][i]
    print(mat)
    print(x)

def GE_with_pivoting(mat,n):
    I=numpy.eye(n)
    for i in range(0,n-1):
        count=i
        for j in range(i+1,n):
            if mat[count][i]<mat[j][i]:

                count=j
        newpivot=mat[count]
        newpivoti=deepcopy(I[count])
        mat[count]=mat[i]
        I[count]=I[i]
        mat[i]=newpivot
        I[i]=newpivoti
    print(I)
    print(mat)
    GE_without_pivoting(mat,n)

def LU_decomposition_by_using_GE_without_pivoting(U,n):
    L= [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        L[i][i]=1
    for i in range(0,n-1):
        for j in range(i+1,n):
            alpha=-(U[j][i])/U[i][i]
            L[j][i]=-1*alpha
            U[j][i]=0
            for k in range(i+1,n):
                U[j][k]=U[j][k]+alpha*U[i][k]

    #Forward substitution
    y=[0]*n
    for i in range(0,n):
        if i==0:
            y[i]=U[0][n]
            continue
        sumi=0
        for j in range(0,i):
            sumi=sumi+L[i][j]*y[j]
        y[i]=(U[i][n]-sumi)

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=y[n-1]/U[n-1][n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+U[i][j]*x[j]
        x[i]=(y[i]-sumi)/U[i][i]
    print('X=\n',x)
    print('L=\n',L)
    print('U=\n',U)

def LU_decomposition_by_using_GE_with_pivoting(U,n):
    for i in range(0,n-1):
        count=i
        for j in range(i+1,n):
            if U[count][i]<U[j][i]:
                count=j
        newpivot=U[count]
        U[count]=U[i]
        U[i]=newpivot
    LU_decomposition_by_using_GE_without_pivoting(U,n)

def LU_decomposition_by_using_crout_without_pivoting(U,n):
    L= [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        L[i][i]=U[i][i]
        for t in range(0,n):
            U[i][t]=U[i][t]/L[i][i]
        for j in range(i+1,n):
            L[j][i]=U[j][i]
            U[j][i]=0
            for k in range(i+1,n):
                U[j][k]=U[j][k]-L[j][i]*U[i][k]
    print(L)
    print(U)

    #Forward substitution
    y=[0]*n
    for i in range(0,n):
        if i==0:
            y[i]=U[0][n]/L[0][0]
            continue
        sumi=0
        for j in range(0,i):
            sumi=sumi+L[i][j]*y[j]
        y[i]=(U[i][n]-sumi)/L[i][i]

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=y[n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+U[i][j]*x[j]
        x[i]=(y[i]-sumi)
    print(x)

def Cholesky_decomposition(U,n):
    L= [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        L[i][i]=sqrt(U[i][i])
        for k in range(n):
            U[i][k]=U[i][k]/L[i][i]
        for j in (i,n-1):
            L[j][i]=U[j][i]/U[i][i]
            for k in range(n):
                U[j][k]=U[j][k]-L[j][i]*U[i][k]
    print(L)
    print(U)

if q=='a':
    GE_without_pivoting(a,n)
if q=='b':
    GE_with_pivoting(a,n)
if q=='c':
    GE_with_pivoting(a,n)
if q=='d':
    LU_decomposition_by_using_GE_without_pivoting(a,n)
if q=='e':
    LU_decomposition_by_using_GE_with_pivoting(a,n)
if q=='f':
    LU_decomposition_by_using_crout_without_pivoting(a,n)
if q=='g':
    Cholesky_decomposition(a,n)