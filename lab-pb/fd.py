import numpy as np
from math import pi
import scipy.sparse as sp

#def fd_sym(N, R, L, Q_pos, Q_f, eps):
N=100
R=1
L=10
Q_pos=10
Q_f=10
eps=1
rr = np.linspace(R, L, N)
dr = (L-R)/N

S_d=4*pi*R**2
I_d=4*pi
sigma_f=Q_f/S_d/eps

phi0=np.ones(N)
C0=Q_pos/(((np.exp(-phi0)*rr**2).sum()-np.exp(-phi0[0])*rr[0]**2/2-np.exp(-phi0[-1])*rr[-1]**2/2)*dr)/I_d

A=sp.diags([-1, 2, -1], [-1, 0, 1], shape=(N, N))/dr**2*eps
B=sp.diags([-np.concatenate(([0], 1/rr[1:-1])), np.concatenate((1/rr[1:-1],[0]))], [-1, 1])/dr*eps
A=A.tocsr()
A[ 0, 0] =  eps*2/dr**2
A[ 0, 1] = -eps*2/dr**2
A[-1,-1] =  eps*2/dr**2
A[-1,-2] = -eps*2/dr**2

tol=1e-7
k=0
ep=1
phi=phi0
C=C0
while(ep>tol):
    P=C*(np.exp(-phi)+np.exp(phi))
    F=C*(np.exp(-phi)-np.exp(phi)+(np.exp(-phi)+np.exp(phi))*phi)
    F[0]+=(2/dr-2/R)*sigma_f*eps
    Coeff=A+B+sp.diags(P,0)
    phi_num=sp.linalg.spsolve(Coeff, F)
    ep=np.abs(phi_num-phi).max()
    C=Q_pos/(((np.exp(-phi_num)*rr**2).sum()-np.exp(-phi_num[0])*rr[0]**2/2-np.exp(-phi_num[-1])*rr[-1]**2/2)*dr)/I_d
    phi=phi_num
    k+=1
    print(ep)

#   return (phi_num, C)