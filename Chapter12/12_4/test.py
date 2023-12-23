import numpy as np
from matplotlib import pyplot as plt
import random

n=4
p=0.1
amatrix = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i>j:
            x= random.random()
            if x< p:
                amatrix[i][j]= 1
                amatrix[j][i]=1
degree= np.zeros((n,1))
print(amatrix)
for i in range(n):
    degree[i,:] = np.sum(amatrix[i,:])

l = np.full((n,n),-1)
t=1
print(l)
def check_for_off_diagonal_terms (a):
    for i in range(n):
        for j in range(n):
            if i!=j:
                if a[i][j]==-1:
                    return True
                
    return False
def check_for_off_diagonal_terms1(a):
    n = len(a)  # Assuming 'n' is defined somewhere before this function is called
    for i in range(n):
        for j in range(n):
            if i != j and a[i][j] == -1:
                return True  # Return True if any off-diagonal element is -1
    return False
                
int_a = amatrix

# while check_for_off_diagonal_terms1(l):
for t in range(1,n**2-n):
    for i in range(n):
        for j in range(i+1,n):
            if amatrix[i,j]!=0:
                if l[i,j]==l[j,i]==-1:
                    l[i,j]=t
                    l[j,i]=t
                # else:
                #     l[i,j]==-l[i,j]
        
    print(t, l)
    amatrix=np.dot(amatrix,int_a)
    # t+=1
    # print(amatrix , t)

print(l)      
print(t)     
# print(amatrix)