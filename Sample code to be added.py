# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:28:11 2021

@author: Nicholas Goh
"""
from numpy import *
#Cost of piping and ducting:

velocity(i,j) = 
rhog = 
C_ref = 
DIA_ref = 
n = 1.08
CEPCI_ref = 
MF = 
A_f = 
FX_rate = 
B = 
F = 
b = 


for i in units:
    for j in units:
        if j>i:
            
            DIA[i][j]= sqrt( (4*Q[i][j] / (velocity[i][j] * pi * rhog))
            C[i][j] = C_ref * (DIA[i][j]/DIA_ref)**n * (CEPCI_2021/CEPCI_ref) * MF * FX_rate * A_f
            PC[i][j] = B * (DIA[i][j]**n) * (CEPCI_2021/CEPCI_ref) *FX_rate
            C[i][j] = PC[i][j] * (1+F) * (A_f + b)
            
            #overall: inclusive of number of pipes in between units, eg. need two pipes if flow split into two. 
            CC = lpSum(D[i][j] * C[i][j] * n[i][j])
#Cost of pumping:
    k[i][j] = epsilon[i][j] / DIA[i][j]
    Re[j][j] = rhog * velocity[i][j] * DIA[i][j] / visc[i][j]
    A[i][j] = (k[i][j] ** 1.1098)/2.8257 + (7.149 / Re[i][j])**0.8981
    #friction factor
    f[i][j] = (1/ (-2 * log(k[i][j]/3.7065 - (5.0452/Re[i][j])*log(A[i][j]))) )**2
    #Pressure drop
    deltaPD[i][j] = 8 * f[i][j] * (D[i][j]/DIA[i][j]) *(phog/2) * velocity[i][j]**2
    #Pumping power (mechEffic = mechanical efficiency.)
    P[i][j] = Q[i][j] * deltaPD[i][j] / (rhog * mechEffic)
    #operating cost for the year
    POC[i][j] = C_elec * OH * P[i][j] * n[i][j]