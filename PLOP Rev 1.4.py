# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:02:59 2020

@author: sj1617
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 14:31:20 2020

@author: sj1617
"""

'testing'

"""

Continuous-domain representation of the process plant layout optimisation
problem with Dow F&EI constraints

Illustrative example adapted from D. I. PATSIATZIS et. al., 2004

Supplementary equations from Dow's F&EI Guide

Authors: Rey Tan, Chris Tighe 2020

Notes:
- Definition of units and pertinent_units sets needs to be easily changed in a GUI.
- Consider adding switches to toggle terms in objective function while constraints are still in problem

"""
#%% --------------Set-up--------------
# Import functions
from pulp import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatch
from scipy.stats import norm
import math
import csv
import time

# Create the model object
layout = pulp.LpProblem("Layout_Problem_Model",LpMinimize)

#%% --------------Errors--------------
class NoFEIError(Exception):
    """Exception raised when protection devices enabled without FEI constraints"""
pass #pass here as class has no definition so avoids error 

#%%--------------Switches--------------
# CPLEX free edition only supports up to 1000 variables. For larger land sizes, use CBC or increase coarseness (g)
# 1 = CPLEX, 2 = CBC
solver = 1

# Toggle constraints in layout problem (1 is on; 0 is off)

# Land Shape Constraints (1 for non-square polygon, 0 for square based on xmax and ymax)
SwitchLandShape = 1

# Toggle Minimum Separation Distances switch
SwitchMinSepDistance = 1


# FEI Constraints
SwitchFEI = 0

# Toggle protection devices (must have FEI enabled if 1)
SwitchProt = 0

# FEI Cost of Life Constraints (only works if SwitchFEI is on)
SwitchFEIVle = 1

# CEI Constraints
SwitchCEI = 0

# Check for errors if SwitchProt == 1 and SwitchFEI == 0:
#    raise NoFEIError("FEI cost of life constraints not allowed without FEI constraints")        
if SwitchProt == 1 and SwitchFEI == 0:
    raise NoFEIError("Protection devices not allowed without FEI constraints") 

#%% --------------Define Sets--------------
# Define the process units
units = ['reactor', 'hex1', 'eoabs', 'hex2', 'co2abs', 'flash', 'pump']
pertinent_units = ['reactor', 'eoabs', 'co2abs']
hazardous_chemicals = ['EO']
Nunits = len(units)


#%% --------------Define Parameters and Values--------------
# Base layout model

# M (m) is used in the contraints which set A or B, L or R, and Din or Dout.
# It should be big enough not to constrain the size of the plot, but not too big ##for use of big M method?
M = 1e3

#polygon layout:

X_begin = np.array([0,0,15,30,30])
X_end = np.array([0,15,30,30,0])
Ng = max(X_end)
XDiff = X_end - X_begin

Y_begin = np.array([0,21,35,21,0])
Y_end = np.array([21,35,21,0,0])
YDiff = Y_end - Y_begin

#for plotting:
X_beginplot = list(X_begin)
X_endplot = list(X_end)
Y_beginplot = list(Y_begin)
Y_endplot = list(Y_end)

#check for convex shape or not

#check for vertical lines
grad0_list = list(reversed(list(np.where(XDiff==0) [0])))

if len(grad0_list)>0:
    X_begin = list(X_begin)
    X_end = list(X_end)
    Y_begin = list(Y_begin)
    Y_end = list(Y_end)
    
    for i in grad0_list:
        del X_begin[i]
        del X_end[i]
        del Y_begin[i]
        del Y_end[i]
    
    X_begin = np.asarray(X_begin)
    X_end = np.asarray(X_end)
    Y_begin = np.asarray(Y_begin)
    Y_end = np.asarray(Y_end)

XDiff = X_end - X_begin
YDiff = Y_end - Y_begin      
grad = YDiff/XDiff
absgrad = abs(grad)
colin = Y_begin - (grad * X_begin)

# Converting to list for manipulation.
X_begin = list(X_begin)
X_end = list(X_end)
Y_begin = list(Y_begin)
Y_end = list(Y_end)
grad = list(grad)
absgrad = list(absgrad)
colin = list(colin)

#Area calculation if needed??
Area = np.trapz(Y_begin, X_begin)


# Dimensions of each unit (m)
alpha = dict.fromkeys(units)
beta = dict.fromkeys(units)
alpha['reactor'] = 10
alpha['hex1'] = 11.42
alpha['eoabs'] = 7.68
alpha['hex2'] = 8.48
alpha['co2abs'] = 5
alpha['flash'] = 2.68
alpha['pump'] = 2.40
beta['reactor'] = 5
beta['hex1'] = 11.42
beta['eoabs'] = 7.68
beta['hex2'] = 8.48
beta['co2abs'] = 5
beta['flash'] = 2.68
beta['pump'] = 2.40

# Purchase cost of each unit (dollars)
Cp = dict.fromkeys(units)
Cp['reactor'] = 335000
Cp['hex1'] = 11000
Cp['eoabs'] = 107000
Cp['hex2'] = 4000
Cp['co2abs'] = 81300
Cp['flash'] = 5000
Cp['pump'] = 1500

#Minimum separation distances
Demin = np.zeros((len(units), len(units)))
Demin[0][1] = 1
Demin[1][2] = 1
Demin[5][6] = 20

Demin = Demin + Demin.T - np.diag(Demin.diagonal())
Demin = makeDict([units,units],Demin,0)

## define the velocities 
velocity = np.zeros((len(units), len(units)))
velocity[0][1] = 1  # connection cost between reactor and hex1
velocity[0][4] = 1  # connection cost between reactor and co2abs
velocity[1][2] = 1  # connection cost between hex1 and eoabs
velocity[2][3] = 1  # connection cost between eoabs and hex2
velocity[3][4] = 1  # connection cost between hex2 and co2abs
velocity[4][5] = 1  # connection cost between co2abs and flash
velocity[4][6] = 1  # connection cost between co2abs and pump
velocity[5][6] = 1  # connection cost between flash and pump

velocity = velocity + velocity.T - np.diag(velocity.diagonal())#
velocity = makeDict([units,units],velocity,0)

## define the flowrates Q
Q = np.zeros((len(units), len(units)))
Q[0][1] = 100  # connection cost between reactor and hex1
Q[0][4] = 100  # connection cost between reactor and co2abs
Q[1][2] = 100  # connection cost between hex1 and eoabs
Q[2][3] = 100 # connection cost between eoabs and hex2
Q[3][4] = 100  # connection cost between hex2 and co2abs
Q[4][5] = 100  # connection cost between co2abs and flash
Q[4][6] = 100  # connection cost between co2abs and pump
Q[5][6] = 100  # connection cost between flash and pump

Q = Q + Q.T - np.diag(Q.diagonal())
Q = makeDict([units,units],Q,0)

## define epsilon 
epsilon = np.zeros((len(units), len(units)))
epsilon[0][1] = 98.4  # connection cost between reactor and hex1
epsilon[0][4] = 98.4  # connection cost between reactor and co2abs
epsilon[1][2] = 98.4  # connection cost between hex1 and eoabs
epsilon[2][3] = 98.4  # connection cost between eoabs and hex2
epsilon[3][4] = 98.4  # connection cost between hex2 and co2abs
epsilon[4][5] = 98.4  # connection cost between co2abs and flash
epsilon[4][6] = 98.4  # connection cost between co2abs and pump
epsilon[5][6] = 98.4  # connection cost between flash and pump

epsilon = epsilon + epsilon.T - np.diag(epsilon.diagonal())
epsilon = makeDict([units,units],epsilon,0)

## define viscosity 
visc = np.zeros((len(units), len(units)))
visc[0][1] = 98.4  # connection cost between reactor and hex1
visc[0][4] = 98.4  # connection cost between reactor and co2abs
visc[1][2] = 98.4  # connection cost between hex1 and eoabs
visc[2][3] = 98.4  # connection cost between eoabs and hex2
visc[3][4] = 98.4  # connection cost between hex2 and co2abs
visc[4][5] = 98.4  # connection cost between co2abs and flash
visc[4][6] = 98.4  # connection cost between co2abs and pump
visc[5][6] = 98.4  # connection cost between flash and pump

visc = visc + visc.T - np.diag(visc.diagonal())
visc = makeDict([units,units],visc,0)

rhog = np.zeros((len(units), len(units)))
rhog[0][1] = 98.4  # connection cost between reactor and hex1
rhog[0][4] = 98.4  # connection cost between reactor and co2abs
rhog[1][2] = 98.4  # connection cost between hex1 and eoabs
rhog[2][3] = 98.4  # connection cost between eoabs and hex2
rhog[3][4] = 98.4  # connection cost between hex2 and co2abs
rhog[4][5] = 98.4  # connection cost between co2abs and flash
rhog[4][6] = 98.4  # connection cost between co2abs and pump
rhog[5][6] = 98.4  # connection cost between flash and pump

rhog = rhog + rhog.T - np.diag(rhog.diagonal())
rhog = makeDict([units,units],rhog,0)

npp = np.zeros((len(units), len(units)))
npp[0][1] = 1  # connection cost between reactor and hex1
npp[0][4] = 1  # connection cost between reactor and co2abs
npp[1][2] = 1  # connection cost between hex1 and eoabs
npp[2][3] = 1  # connection cost between eoabs and hex2
npp[3][4] = 1  # connection cost between hex2 and co2abs
npp[4][5] = 1  # connection cost between co2abs and flash
npp[4][6] = 1  # connection cost between co2abs and pump
npp[5][6] = 1  # connection cost between flash and pump

npp = npp + npp.T - np.diag(npp.diagonal())
npp = makeDict([units,units],npp,0)

#define constants for piping cost 

C_ref = 1
n_1 = 1.08
n_2 = 1
CEPCI_ref  =1 
MF = 1
A_f = 11
FX_rate = 1
BB = 1
F = 1
bb = 1
CEPCI_2021 =1
DIA_ref = 1
mechEffic = 0.6
C_elec = 0.000045
OH = 8000

           
            # kr[i][j] += epsilon[i][j] / DIA[i][j]
            # Rey[i][j] += rhog[i][j] * velocity[i][j] * DIA[i][j] / visc[i][j]

            # AA[i][j] += (kr[i][j] ** 1.1098)/2.8257 + (7.149 / Rey[i][j])**0.8981
            # ff[i][j] += (1/ (-2 * math.log(kr[i][j]/3.7065 - (5.0452/Rey[i][j])*math.log(AA[i][j]))) )**2

            # deltaP[i][j] += 8 * ff[i][j] * (D[i][j]/DIA[i][j]) *(rhog/2) * velocity[i][j]**2
            # PP[i][j] += Q[i][j] * deltaP[i][j] / (rhog * mechEffic)
            # POC[i][j] += C_elec * OH * PP[i][j] * n[i][j]

#%% Land shape constraint: 1 if non-rectangular, 0 if rectangular.
if SwitchLandShape == 0: #sets default max available plot area.
    xmax = 40
    ymax = 40


#%% ----------F&EI model parameter input-----------
# Operating condidions
T = dict.fromkeys(pertinent_units)
Pg = dict.fromkeys(pertinent_units)
w = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
rhol = dict.fromkeys(pertinent_units)
MWunit = dict.fromkeys(pertinent_units)
INV = dict.fromkeys(pertinent_units)
# Assign values
T['reactor'] = 275  # degc
T['eoabs'] = 40
T['co2abs'] = 20
Pg['reactor'] = 2200  # kPa
Pg['eoabs'] = 2000
Pg['co2abs'] = 500
w['reactor']['EO'] = 0.94  # mass frac
w['eoabs']['EO'] = 0.995
w['co2abs']['EO'] = 0.05
rhol['reactor'] = 1000
rhol['eoabs'] = 1000
rhol['co2abs'] = 1000
MWunit['reactor'] = 44.05  # g/mol #EO=ethylene oxide
MWunit['eoabs'] = 44.05
MWunit['co2abs'] = 44.05
INV['reactor'] = 10000  # kg
INV['eoabs'] = 10000
INV['co2abs'] = 10000
#%% ----------- SwitchFEI---------
if SwitchFEI == 1:
    # F&EI factors
#    MF = dict.fromkeys(pertinent_units)
#    F1 = dict.fromkeys(pertinent_units)
#    F2 = dict.fromkeys(pertinent_units)
#    F = dict.fromkeys(pertinent_units)
    De = dict.fromkeys(pertinent_units)
    DF = dict.fromkeys(pertinent_units)
#    # assign values
#    MF['reactor'] = 29
#    MF['eoabs'] = 29
#    MF['co2abs'] = 24
#    F1['reactor'] = 2.2
#    F1['eoabs'] = 1.2
#    F1['co2abs'] = 1.2
#    F2['reactor'] = 2.45
#    F2['eoabs'] = 2.45
#    F2['co2abs'] = 2.45
    # compute factors
#    F3 = {i: F1[i]*F2[i] for i in pertinent_units}
#    F = {i: MF[i]*F3[i] for i in pertinent_units}
#    De = {i: 0.256*F[i] for i in pertinent_units}

    DF['reactor'] = 0.87
    DF['eoabs'] = 0.73
    DF['co2abs'] = 0.66
    De['reactor'] = 40
    De['eoabs'] = 21.8
    De['co2abs'] = 18.06

# damage factor equations based on value of MF, F3[i] = 8 used when it is greater than 8
#    def DFcalc(value):
#         if value == 1:
#             if F3[i] > 8:
#                 return 0.003907 + 0.002957*8 + 0.004031*8**2 - 0.00029*8**3
#             return 0.003907 + 0.002957*F3[i] + 0.004031*F3[i]**2 - 0.00029*F3[i]**3
#         elif value == 4:
#             if F3[i] > 8:
#                 return 0.025817 + 0.019071*8 -0.00081*8**2 + 0.000108*8**3
#             return 0.025817 + 0.019071*F3[i] -0.00081*F3[i]**2 + 0.000108*F3[i]**3
#         elif value == 10:
#             if F3[i] > 8:
#                 return 0.098582 + 0.017596*8 + 0.000809*8**2 - 0.000013*8**3
#             return 0.098582 + 0.017596*F3[i] + 0.000809*F3[i]**2 - 0.000013*F3[i]**3
#         elif value == 14:
#             if F3[i] > 8:
#                 return 0.20592 + 0.018938*8 + 0.007628*8**2 - 0.00057*8**3
#             return 0.20592 + 0.018938*F3[i] + 0.007628*F3[i]**2 - 0.00057*F3[i]**3
#         elif value == 16:
#             if F3[i] > 8:
#                 return 0.256741 + 0.019886*8 + 0.011055*8**2 - 0.00088*8**3
#             return 0.256741 + 0.019886*F3[i] + 0.011055*F3[i]**2 - 0.00088*F3[i]**3
#         elif value == 21:
#             if F3[i] >8:
#                 return 0.340314 + 0.076531*8 + 0.003912*8**2 - 0.00073*8**3
#             return 0.340314 + 0.076531*F3[i] + 0.003912*F3[i]**2 - 0.00073*F3[i]**3
#         elif value == 24:
#             if F3[i] > 8:
#                 return 0.395755 + 0.096443*8 - 0.00135*8**2 - 0.00038*8**3
#             return 0.395755 + 0.096443*F3[i] - 0.00135*F3[i]**2 - 0.00038*F3[i]**3
#         elif value == 29:
#             if F3[i] > 8:
#                 return 0.484766 + 0.094288*8 - 0.00216*8**2 - 0.00031*8**3
#             return 0.484766 + 0.094288*F3[i] - 0.00216*F3[i]**2 - 0.00031*F3[i]**3
#         elif value == 40:
#             if F3[i] > 8:
#                 return 0.554175 + 0.080772*8 + 0.000332*8**2 - 0.00044*8**3
#             return 0.554175 + 0.080772*F3[i] + 0.000332*F3[i]**2 - 0.00044*F3[i]**3
#         else:
#             raise ValueError(value)
#    for i in pertinent_units:
#         DF[i] = DFcalc(MF[i])

    # With literature values for DF and De

    # Upper bound for actual maximum probable property damage cost
    U = 1e8
#%%--------------- SwitchProt-----------------------------
    # Protection device model
    if SwitchProt == 1:
        # Define protection device configuration
        configurations = list(range(1, len(units)))
        # Loss control credit factor of protection device configuration k on item i
        gamma = np.zeros((len(pertinent_units), len(configurations)))
        gamma = makeDict([pertinent_units,configurations],gamma,0)
        # assign values
        gamma['reactor'][1] = 1
        gamma['reactor'][2] = 0.900
        gamma['reactor'][3] = 0.750
        gamma['reactor'][4] = 0.365
        gamma['reactor'][5] = 0.292
        gamma['reactor'][6] = 0.117
        gamma['eoabs'][1] = 1
        gamma['eoabs'][2] = 0.900
        gamma['eoabs'][3] = 0.760
        gamma['eoabs'][4] = 0.684
        gamma['eoabs'][5] = 0.612
        gamma['eoabs'][6] = 0.465
        gamma['co2abs'][1] = 1
        gamma['co2abs'][2] = 0.900
        gamma['co2abs'][3] = 0.760
        gamma['co2abs'][4] = 0.684
        gamma['co2abs'][5] = 0.612
        gamma['co2abs'][6] = 0.465
        # purchase cost of configuration k on unit i
        P = np.zeros((len(pertinent_units),len(configurations)))
        P = makeDict([pertinent_units,configurations],P,0)
        # assign values
        P['reactor'][1] = 0
        P['reactor'][2] = 5000
        P['reactor'][3] = 15000
        P['reactor'][4] = 40000
        P['reactor'][5] = 60000
        P['reactor'][6] = 125000
        P['eoabs'][1] = 0
        P['eoabs'][2] = 5000
        P['eoabs'][3] = 20000
        P['eoabs'][4] = 25000
        P['eoabs'][5] = 35000
        P['eoabs'][6] = 55000
        P['co2abs'][1] = 0
        P['co2abs'][2] = 5000
        P['co2abs'][3] = 20000
        P['co2abs'][4] = 25000
        P['co2abs'][5] = 35000
        P['co2abs'][6] = 55000

#%%--------------------CEI Factors parameters---------------------------
# # Initialise dictionaries
AQ = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced
AQf = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced by flash
AQp = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced by pool
Fv = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # fraction flashed
Pvap = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
CEI = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
Dc = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
Dh = dict.fromkeys(pertinent_units)
Liq = dict.fromkeys(pertinent_units)
Deltah = dict.fromkeys(pertinent_units)
WT = dict.fromkeys(pertinent_units)
WP = dict.fromkeys(pertinent_units)
AP = dict.fromkeys(pertinent_units)
heatratio = dict.fromkeys(hazardous_chemicals)
Tb = dict.fromkeys(hazardous_chemicals)
antoineA = dict.fromkeys(hazardous_chemicals)
antoineB = dict.fromkeys(hazardous_chemicals)
antoineC = dict.fromkeys(hazardous_chemicals)
MW = dict.fromkeys(hazardous_chemicals)
Pr = dict.fromkeys(hazardous_chemicals)
prob_death = dict.fromkeys(hazardous_chemicals)
CONC = dict.fromkeys(hazardous_chemicals)
# # Assign values
te = 10  # minutes of exposure time
prob_death['EO'] = 0.5  # probability of death for probit
Dh['reactor'] = 50.8  # mm
Dh['eoabs'] = 50.8  # hole diameter
Dh['co2abs'] = 50.8
Deltah['reactor'] = 1  # m
Deltah['eoabs'] = 1
Deltah['co2abs'] = 1
heatratio['EO'] = 0.00365  # 1/degc
Tb['EO'] = 10.5  # degc
antoineA['EO'] = 4.386
antoineB['EO'] = 1115.1
antoineC['EO'] = -29.015
MW['EO'] = 44.05

#%%-------------- Compute CEI factors -------------------------
#     # Probit equation
for h in hazardous_chemicals:
    Pr[h] = norm.ppf(prob_death[h]) + 5
    CONC[h] = math.sqrt(math.exp(Pr[h] + 17.5)/te)  # mg/m^3
for i in pertinent_units:
    for h in hazardous_chemicals:
        # Vapour pressure (kPa)
        Pvap[i][h] = 100 * 10**(antoineA[h] - antoineB[h]/(antoineC[h] + (Tb[h] + 273.15)))
for i in pertinent_units:
    # Liquid flowrate (kg/s)
    Liq[i] = 9.44E-07 * Dh[i]**2 * rhol[i] * math.sqrt(1000*Pg[i]/rhol[i] + 9.8*Deltah[i])
    if 300*Liq[i] >= INV[i]:
        Liq[i] = INV[i]/300
    else:
         Liq[i] = 9.44E-07 * Dh[i]**2 * rhol[i] * math.sqrt(1000*Pg[i]/rhol[i] + 9.8*Deltah[i])
    # Total liquid release (kg)
    if Liq[i] < INV[i]/900:
        WT[i] = 900*Liq[i]
    elif Liq[i] >= INV[i]/900:
        WT[i] = INV[i]
    else:
        pass
    # Fraction flashed
    for h in hazardous_chemicals:
        Fv[i][h] = w[i][h] * heatratio[h] * (T[i] - Tb[h])
        if Fv[i][h]/w[i][h] < 0.2:
            # Air quantity from flashed (kg/s)
            AQf[i][h] = 5*Fv[i][h]*Liq[i]
            # Pool size (kg)
            WP[i] = WT[i] * (1 - 5*Fv[i][h]/w[i][h])
        elif Fv[i][h]/w[i][h] >= 0.2:
            AQf[i][h] = w[i][h]*Liq[i]
            WP[i] = 0
        else:
            pass
        if AQf[i][h] == 0:
            WP[i] = WT[i]
        else:
            pass
    # Area of pool (m^2)
    AP[i] = 100 * WP[i]/rhol[i]
    for h in hazardous_chemicals:
        # Air quantity from evaporating pool (kg/s)
        AQp[i][h] = 9E-04 * AP[i]**0.95 * ((MW[h] * Pvap[i][h] * w[i][h])/(Tb[h] + 273))
        # Total air quantity
        AQ[i][h] = AQf[i][h] + AQp[i][h]
        # Chemical  exposure index
        CEI[i][h] = 655.1 * math.sqrt(AQ[i][h]/CONC[h])
        # Hazard distance (m)
        Dc[i][h] = 6551 * math.sqrt(AQ[i][h]/CONC[h])
        if Dc[i][h] > 10000:
            Dc[i][h] = 10000
# Find max chemical hazard distance
Dcvalues = []
for a, b in Dc.items():
    for c, d in b.items():
        Dcvalues.append(d)
maxDc = max(Dcvalues)

#%%----------------- SwitchFEIVle -----------------------
# # Occupancy calculations
# # Number of workers
# if SwitchFEIVle == 1:
Nw = dict.fromkeys(units)
Nw['reactor'] = 3
Nw['hex1'] = 1
Nw['eoabs'] = 3
Nw['hex2'] = 1
Nw['co2abs'] = 3
Nw['flash'] = 1
Nw['pump'] = 1
# Percentage of time present at unit
tw = dict.fromkeys(units)
tw['reactor'] = 0.10
tw['hex1'] = 0.05
tw['eoabs'] = 0.10
tw['hex2'] = 0.05
tw['co2abs'] = 0.10
tw['flash'] = 0.05
tw['pump'] = 0.05
# Occupancy
OCC = dict.fromkeys(units)
OCC = {i: Nw[i]*tw[i] for i in units}
# Cost of life (dollars)
Cl = dict.fromkeys(units)
for i in units:
    Cl[i] = 10e6

#%% --------------Define Variables--------------
# Base layout model

# 1 if length of item i is equal to alpha; 0 otherwise
O = LpVariable.dicts("O",(units),lowBound=0,upBound=1,cat="Integer")
# pair-wise values to ensure equipment items do not overlap
E1 = LpVariable.dicts("E1",(units,units),lowBound=0,upBound=1,cat="Integer")
E2 = LpVariable.dicts("E2",(units,units),lowBound=0,upBound=1,cat="Integer")
# is unit i to the right or above unit j
Wx = LpVariable.dicts("Wx",(units,units),lowBound=0,upBound=1,cat="Integer")
Wy = LpVariable.dicts("Wy",(units,units),lowBound=0,upBound=1,cat="Integer")

# Define continuous variables for base layout model
l = LpVariable.dicts("l",(units),lowBound=0,upBound=None,cat="Continuous")
    # breadth of item i
d = LpVariable.dicts("d",(units),lowBound=0,upBound=None,cat="Continuous")
    # relative distance in x coordinates between items i and j, if i is to the right of j
R = LpVariable.dicts("R",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # relative distance in x coordinates between items i and j, if i is to the left of j
L = LpVariable.dicts("L",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # relative distance in y coordinates between items i and j, if i is above j
A = LpVariable.dicts("A",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # relative distance in y coordinates between items i and j, if i is below j
B = LpVariable.dicts("B",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # total rectilinear distance between items i and j
D = LpVariable.dicts("D",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # coordinates of the geometrical centre of item i
x = LpVariable.dicts("X",(units),lowBound=None,upBound=None,cat="Continuous")
y = LpVariable.dicts("Y",(units),lowBound=None,upBound=None,cat="Continuous")

# Define continuousnterunit connection cost
CD = LpVariable.dicts("CD",(units,units),lowBound=0,upBound=None,cat="Continuous")
# Total connection cost
SumCD = LpVariable("SumCD",lowBound=0,upBound=None,cat="Continuous")

# pressure drop between pipes
deltaP = LpVariable.dicts("deltaP",(units,units),lowBound=0,upBound=None,cat="Continuous")
DIA = np.zeros((len(units), len(units)))
C = np.zeros((len(units), len(units)))
PC = np.zeros((len(units), len(units)))
C_annual = np.zeros((len(units), len(units)))
DIA = makeDict([units,units],DIA,0)
C = makeDict([units,units],C,0)
PC = makeDict([units,units],PC,0)
C_annual = makeDict([units,units],C_annual,0)

kr = np.zeros((len(units), len(units)))
kr = makeDict([units,units],kr,0)
Rey = np.zeros((len(units), len(units)))
Rey = makeDict([units,units],Rey,0)
AA = np.zeros((len(units), len(units)))
AA = makeDict([units,units],AA,0)
ff = np.zeros((len(units), len(units)))
ff = makeDict([units,units],ff,0)
# deltaP = np.zeros((len(units), len(units)))
# deltaP = makeDict([units,units],deltaP,0)
# # PP = np.zeros((len(units), len(units)))
# PP = makeDict([units,units],PP,0)
# # POC = np.zeros((len(units), len(units)))
# POC = makeDict([units,units],POC,0)

for idxj, j in enumerate(units):
    for idxi, i in enumerate(units):
        if idxj > idxi:
            DIA[i][j] += np.sqrt( (4*Q[i][j] / (velocity[i][j] * np.pi * rhog[i][j])))
            C[i][j] += C_ref*(DIA[i][j]*DIA_ref)**n_1 * (CEPCI_2021/CEPCI_ref) * MF * FX_rate * A_f
            PC[i][j] += BB * (DIA[i][j]**n_2) * (CEPCI_2021/CEPCI_ref) *FX_rate
            C_annual[i][j] += PC[i][j] * (1+F) * (A_f + bb)
            kr[i][j] += epsilon[i][j] / DIA[i][j]
            Rey[i][j] += rhog[i][j] * velocity[i][j] * DIA[i][j] / visc[i][j]

            AA[i][j] += (kr[i][j] ** 1.1098)/2.8257 + (7.149 / Rey[i][j])**0.8981
            ff[i][j] += (1/ (-2 * math.log(kr[i][j]/3.7065 - (5.0452/Rey[i][j])*math.log(AA[i][j]))) )**2





if SwitchFEI == 1:
    # 1 if j is allocated within the area of exposure if i; 0 otherwise
    Psie = LpVariable.dicts("Psie",(pertinent_units,units),lowBound=0,upBound=1,cat="Integer")
    # total rectilinear distance if D < De
    Din = LpVariable.dicts("Din",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
    # total rectilinear distance if D > De
    Dout = LpVariable.dicts("Dout",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
    # value of area of exposure of incident on i
    Ve = LpVariable.dicts("Ve",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")
    # second term for value of area of exposure of incident on i
    #Ve2 = LpVariable.dicts("Ve2",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")
    Ve2 = LpVariable.dicts("Ve2",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")

    # base maximum probable property damage cost for pertinent unit i
    Omega0 = LpVariable.dicts("Omega0",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")
    
# Total actual maximum probable property damage cost
SumOmega = LpVariable("SumOmega",lowBound=0,upBound=None,cat="Continuous")

if SwitchProt == 1:
    # 1 if protection device configuration k is installed on item i; 0 otherwise
    Z = LpVariable.dicts("Z",(pertinent_units,configurations),lowBound=0,upBound=1,cat="Integer")
    # actual maximum probable property damage cost for pertinent unit i
    Omega = LpVariable.dicts("Omega",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")
    # linearisation variable denoting product of Omega0 and Z
    linOmega0 = LpVariable.dicts("linOmega0",(pertinent_units,configurations),lowBound=0,upBound=None,cat="Continuous")
    # Total protection devices cost
    SumPZ = LpVariable("SumPZ",lowBound=0,upBound=None,cat="Continuous")
else:
    SumPZ = 0

#if SwitchFEIVle == 1:
#    # value of life in fire and explosion area of exposure
#    Vle = LpVariable.dicts("Vle",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")
#    Vle2 = LpVariable.dicts("Vle2",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
#    # Total cost of life due to fire and explosion
#    SumVle = LpVariable("SumVle",lowBound=0,upBound=None,cat="Continuous")
#    SumVle2 = LpVariable.dicts("SumVle2",(pertinent_units),lowBound=0,upBound=None,cat="Continuous")



if SwitchCEI == 1:
    # Currently unused variables for cost of life in FEI and CEI
    # 1 if j is allocated within chemical the area of exposure if i; 0 otherwise
    Psic = LpVariable.dicts("Psic",(pertinent_units,units,hazardous_chemicals),lowBound=0,upBound=1,cat="Integer")
    Dc_in = LpVariable.dicts("Dc_in",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
    Dc_out = LpVariable.dicts("Dc_out",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
    # value of life in chemical area of exposure
    Vlc = LpVariable.dicts("Vlc2",(pertinent_units,units),lowBound=0,upBound=None,cat="Continuous")
    Vlc2 = LpVariable.dicts("Vlc2",(pertinent_units,units,hazardous_chemicals),lowBound=0,upBound=None,cat="Continuous")
# Total cost of life due to chemical release
SumVlc = LpVariable("SumVlc",lowBound=0,upBound=None,cat="Continuous")

#%% --------------Define Objective Function--------------

obj_sumOmega = SwitchFEI*SumOmega
obj_PZ = SwitchProt*SumPZ
# obj_Vle = SwitchFEIVle*SumVle
obj_Vlc = SwitchCEI*SumVlc
obj_CD = SumCD
layout += obj_CD + obj_sumOmega + obj_PZ + obj_Vlc



#%% --------------Define Constraints and Objective Function Contributions--------------
#%% Base model constraints for all units i
for i in units:
    # Orientation constraints (1 - 2)
    layout += l[i] == alpha[i]*O[i] + beta[i]*(1 - O[i])
    layout += d[i] == alpha[i] + beta[i] - l[i]
    # Lower bounds of coordinates (19 - 22)
    layout += x[i] >= 0.5*l[i]
    layout += y[i] >= 0.5*d[i]
        

for idxj, j in enumerate(units):
    for idxi, i in enumerate(units):
        if idxj > idxi:
            
            
            # kr[i][j] += epsilon[i][j] / DIA[i][j]
            # Rey[i][j] += rhog[i][j] * velocity[i][j] * DIA[i][j] / visc[i][j]


            #layout += deltaP[i][j] == 8 * 7 * (D[i][j]/8) *(rhog/2) * 15**2
            # PP[i][j] += Q[i][j] * deltaP[i][j] / (rhog * mechEffic)
            # POC[i][j] += C_elec * OH * PP[i][j] * n[i][j]
            
            layout += CD[i][j] == C_annual[i][j] * D[i][j]*npp[i][j]
            # Distance calculation constraints (3 - 10)
            layout += R[i][j] - L[i][j] == x[i] - x[j]
            layout += A[i][j] - B[i][j] == y[i] - y[j]
            layout += R[i][j] <= M*Wx[i][j]
            layout += L[i][j] <= M*(1 - Wx[i][j])
            layout += A[i][j] <= M*Wy[i][j]
            layout += B[i][j] <= M*(1 - Wy[i][j])
            layout += D[i][j] == R[i][j] + L[i][j] + A[i][j] + B[i][j]
            layout += D[i][j] == D[j][i]
   
            
            # Nonoverlapping constraints (15 - 18)
            # Including switch for minimum separation distances
            if SwitchMinSepDistance == 1:
                layout += x[i] - x[j] + M*(E1[i][j] + E2[i][j]) >= Demin[i][j] + (l[i] + l[j])/2
                layout += x[j] - x[i] + M*(1 - E1[i][j] + E2[i][j]) >= Demin[i][j]+ (l[i] + l[j])/2
                layout += y[i] - y[j] + M*(1 + E1[i][j] - E2[i][j]) >= Demin[i][j]+ (d[i] + d[j])/2
                layout += y[j] - y[i] + M*(2 - E1[i][j] - E2[i][j]) >= Demin[i][j]+ (d[i] + d[j])/2
            else:
                layout += x[i] - x[j] + M*(E1[i][j] + E2[i][j]) >= (l[i] + l[j])/2
                layout += x[j] - x[i] + M*(1 - E1[i][j] + E2[i][j]) >= (l[i] + l[j])/2
                layout += y[i] - y[j] + M*(1 + E1[i][j] - E2[i][j]) >= (d[i] + d[j])/2
                layout += y[j] - y[i] + M*(2 - E1[i][j] - E2[i][j]) >=  (d[i] + d[j])/2
            
            # These constraints ensure consistency in interdependent variables
            layout += L[i][j] == R[j][i]
            layout += R[i][j] == L[j][i]
            layout += A[i][j] == B[j][i]
            layout += B[i][j] == A[j][i]
            layout += Wx[i][j] == 1 - Wx[j][i]
            layout += Wy[i][j] == 1 - Wy[j][i]

# Objective function contribution for base model
layout += SumCD == lpSum([CD[i][j] for i in units for j in units])        

#%% Land shape constraints (or set max plot size if not used)
if SwitchLandShape == 1:

    for i in units:
        #pentagon constraint: # no half so that no part of the unit is outside constraint
        ##need to account for the top corner, imagine sliding triangle along line connected to unit
        for k in colin:
           if k > 0:
               layout += y[i] + 0.5*d[i] + l[i]*absgrad[colin.index(k)]/2 <= grad[colin.index(k)]*x[i] + k
           elif k <= 0:
               layout += y[i] - 0.5*d[i] - l[i]*absgrad[colin.index(k)]/2 >= grad[colin.index(k)]*x[i] + k
        
        layout += x[i] + 0.5*l[i] <= Ng #CHANGE N TO PEAK OF PENTAGON
        layout += x[i] -0.5*l[i] >= 0
               
else:
    for i in units:
        layout += x[i] + 0.5*l[i] <= xmax
        layout += y[i] + 0.5*d[i] <= ymax
        
#%% F&EI Constraints
if SwitchFEI == 1:
    # for all i in pertinent units, j in units, j != i
    for i in units:
        if i in pertinent_units:
            for j in units:
                if i != j:
                    # Area of exposure constraints (26 - 28)
                    layout += Din[i][j] + Dout[i][j] == D[i][j]
                    layout += Din[i][j] <= De[i]*Psie[i][j]
                    layout += Dout[i][j] >= De[i]*(1 - Psie[i][j])
                    layout += Dout[i][j] <= M*(1 - Psie[i][j])
                    # Cost of life (37)
                    #Commented out as appeared to prevent Psie taking value of 1 - check later
                    #layout += Vle[i] == OCC[i]*Cl[i] + lpSum([OCC[j]*Cl[j]*Psie[i][j]])
                    # Summation term in value of area of exposure constraint (30i)                
                    #layout += Ve2[i] == lpSum([Cp[j]*Psie[i][j] - Cp[j]*Din[i][j]/De[i]])
                    if SwitchFEIVle == 1:
                        layout += Ve2[i][j] == Cp[j]*Psie[i][j] + OCC[j]*Cl[j]*Psie[i][j] - Cp[j]*Din[i][j]/De[i] - OCC[j]*Cl[j]*Din[i][j]/De[i]
                    else:
                        layout += Ve2[i][j] == Cp[j]*Psie[i][j] - Cp[j]*Din[i][j]/De[i]
                else:
                    layout += Ve2[i][j] == 0
    # for all i in pertinent units
    for i in pertinent_units:
        # Value of area of exposure constraint (30ii)
        if SwitchFEIVle == 1:
            layout += Ve[i] == Cp[i] + OCC[i]*Cl[i] + lpSum([Ve2[i][j] for j in units])
        else:
            layout += Ve[i] == Cp[i] + lpSum([Ve2[i][j] for j in units])
        # layout += Vle[i] == OCC[i]*Cl[i] + lpSum([Vle2[i][j] for j in units])
        # Base maximum probable property damage cost (31)
        layout += Omega0[i] == DF[i]*Ve[i]

#    if SwitchFEIVle == 1:
#       for i in units:
#            if i in pertinent_units: 
#                for j in units:
#                    if i != j:
#                        layout += Vle2[i][j] == OCC[j]*Cl[j]*Psie[i][j]
#                    else:
#                        layout += Vle2[i][j] == 0
#        for i in units:
#            if i in pertinent_units:
#                for j in units:
#                    layout += SumVle2[i] == lpSum([Vle2[i][j]])
#        for i in pertinent_units:
#            layout += Vle[i] == OCC[i]*Cl[i] #+ SumVle2[i]
#            layout += SumVle == lpSum(Vle[i] for i in pertinent_units)

    if SwitchProt == 1:
        # for all i in pertinent units,k in configurations
        for i in pertinent_units:
            for k in configurations:
                # Linearisation of product of base MPPD and configuration binary variable (34)
                layout += linOmega0[i][k] <= U*Z[i][k]

        # for all i in pertinent units
        for i in pertinent_units:

            # Only one configuration active (33)
            layout += lpSum([Z[i][k] for k in configurations]) == 1
            # Base maximum probable property damage cost is sum of linearisation term (35)
            layout += Omega0[i] == lpSum([linOmega0[i][k] for k in configurations])
            # Actual maximum probable property damage cost (36)
            layout += Omega[i] == lpSum([gamma[i][k]*linOmega0[i][k] for k in configurations])

            # Objective function contribution with protection devices
            layout += SumOmega == lpSum(Omega[i] for i in pertinent_units) 
            layout += SumPZ == lpSum(P[i][k]*Z[i][k] for i in pertinent_units for k in configurations)
            #layout += SumVle == lpSum(Vle[i] for i in pertinent_units)

    else:
        # Objective function contribution without protection devices
        layout += SumOmega == lpSum(Omega0[i] for i in pertinent_units)


#%% CEI Equations
# All commented out while under review - appears to use same Din and Dout as FEI which is not right - separate case
# Unused CEI constraints
# for i in units:
#     if i in pertinent_units:
#         for j in units:
#             if i != j:
#                 # (26)
#                 layout += Din[i][j] + Dout[i][j] == D[i][j]
#                 #for h in hazardous_chemicals:
#                 #    # (42)
#                 #    layout += Dout[i][j] <= M*(2 - Psie[i][j] - Psic[i][j][h])
if SwitchCEI == 1:
    # for all i in pertinent units, j in units, j != i
    for i in units:
        if i in pertinent_units:
            for j in units:
                if i != j:
                    for h in hazardous_chemicals:
                        #Area of exposure constraints (39 - 40)
                        layout += Dc_in[i][j] + Dc_out[i][j] == D[i][j]
                        layout += Dc_in[i][j] <= Dc[i][h]*Psic[i][j][h]
                        layout += Dc_out[i][j] >= Dc[i][h]*(1 - Psic[i][j][h])
                        layout += Dc_out[i][j] <= M*(1 - Psic[i][j][h])
                        layout += Vlc2[i][j][h] == OCC[j]*Cl[j]*Psic[i][j][h] - OCC[j]*Cl[j]*Dc_in[i][j]/Dc[i][h]
                else:
                    layout += Vlc2[i][j][h] == 0


    for i in pertinent_units:
        # Value of chemical area of exposure constraint (43)
        layout += Vlc[i] == OCC[i]*Cl[i] + lpSum([Vlc2[i][j][h] for j in units for h in hazardous_chemicals])

    # Objective function term
    layout += SumVlc == lpSum(Vlc[i] for i in pertinent_units)

# --------------Fixing Variable Values--------------
# Define function to fix value
# def fix_variable(variable, value):
#     variable.setInitialValue(value)
#     variable.fixValue()

# For the purpose of solving case with no protection devices, uncomment out this section if desired
# Fix Z_i,1 = 1, the rest = 0
# if SwitchProt != 1:
#    for k in configurations[1:]:
#        for i in pertinent_units:
#            fix_variable(Z[i][k], 0)
#    for i in pertinent_units:
#        fix_variable(Z[i][1], 1)
# --------------Initialise--------------
# Initialise with optimal solution from reference paper

# solution_x = {
#     ('reactor'): 30,#.39,
#     ('hex1'): 15.39,
#     ('eoabs'): 26.04,
#     ('hex2'): 4.24,
#     ('co2abs'): 4.24,
#     ('flash'): 5.58,
#     ('pump'): 2.9,
# }
# solution_y = {
#     ('reactor'): 30,#.39,
#     ('hex1'): 15.39,
#     ('eoabs'): 4.24,
#     ('hex2'): 4.24,
#     ('co2abs'): 26.54,
#     ('flash'): 9.82,
#     ('pump'): 9.82,
#}


# for i, v in solution_x.items():
#     x[i].setInitialValue(v)

# for i, v in solution_y.items():
#     y[i].setInitialValue(v)

# for i, v in solution_x.items():
#     fix_variable(x[i], v)

# for i, v in solution_y.items():
#     fix_variable(y[i], v)


# fix_variable(Psie['reactor']['hex1'], 1)
# fix_variable(Psie['reactor']['eoabs'], 1)
# fix_variable(Psie['reactor']['hex2'], 1)
# fix_variable(Psie['reactor']['co2abs'], 1)
# fix_variable(Psie['reactor']['flash'], 1)
# fix_variable(Psie['reactor']['pump'], 1)


# for i, v in solution_x.items():
#     fix_variable(x[i],v)

# for i, v in solution_y.items():
#     fix_variable(y[i],v)

#%% --------------Initiate Solve--------------
layout.writeLP("DowFEI.lp")
#CPLEXsolver = CPLEX_PY(msg=1, warmStart=1, gapRel=0, logPath='cplex.log')
CPLEXsolver = CPLEX_PY(msg=1, gapRel=0)
CBCsolver = PULP_CBC_CMD(msg=1)
starttime = time.perf_counter()
if solver == 1:
    layout.solve(CPLEXsolver)
elif solver == 2:
    layout.solve(CBCsolver)
totaltime = time.perf_counter() - starttime
# Print solver status
print("Status: ", LpStatus[layout.status])

#%%--------------Print Results--------------
# Print variable and objective function values
for v in layout.variables():
    print(v.name, "=", v.varValue)
print("Total cost of connections =", SumCD.varValue)

if SwitchFEI == 1:
    print("Total actual MPPD =", SumOmega.varValue)
if SwitchProt == 1:
    print("Total cost of protection devices =", SumPZ.varValue)
#    if SwitchFEIVle == 1:
#        print("Total cost of fatality due to fire/explosion =", SumVle.varValue)
if SwitchCEI == 1:
    print("Total cost of fatality due to chemical release =", SumVlc.varValue)
print ("Total cost of layout =", value(layout.objective))

print("Elapsed time =", totaltime)

#%%--------------Export Results--------------
# filename = 'Optimisation_Plot.csv'
# with open(filename, 'w', newline='') as file:
#     # Write objective function
#     writer = csv.writer(file)
#     writer.writerow(['Objective', value(layout.objective)])
    
#     # Write coordinates
#     #fieldnames = ['unit','x', 'y', 'l', 'd', 'De', 'Dc']
#     fieldnames = ['unit','x', 'y', 'l', 'd', 'De', 'Dc']
#     writer = csv.DictWriter(file, fieldnames=fieldnames)    
#     writer.writeheader()
#     for i in units:
#         #writer.writerow({'unit': i, 'x': x[i].varValue, 'y': y[i].varValue, 'l': l[i].varValue, 'd': d[i].varValue, 'De': De.get(i), 'Dc': Dc.get(i)})
#         writer.writerow({'unit': i, 'x': x[i].varValue, 'y': y[i].varValue, 'l': l[i].varValue, 'd': d[i].varValue, 'De': De.get(i), 'Dc': Dc.get(i)})

# filename = 'Optimisation_Results.csv'
# with open(filename, 'w', newline='') as file:
#     # Write value of all variables
#     writer = csv.writer(file, delimiter=',')    
#     for v in layout.variables():
#         writer.writerow([v.name, v.varValue])

#%%--------------Plot Results--------------
xpos, ypos = [], []
for i in units:
    xpos.append(x[i].varValue)
    ypos.append(y[i].varValue)
    


# Plot invisible scatter
fig, ax = plt.subplots()
ax.scatter(xpos,ypos,alpha=0)
#### LINE THING IS X1,X2 AND Y1,Y2 WTF

#Pentagon:
if SwitchLandShape == 1:
    for k in X_begin:
        line = plt.Line2D((k, X_end[X_begin.index(k)]), (Y_begin[X_begin.index(k)], Y_end[X_begin.index(k)]), lw = 1.5)
        plt.gca().add_line(line)
    if len(grad0_list)>0:  
        for i in grad0_list:
            line = plt.Line2D((X_beginplot[i], X_endplot[i]), (Y_beginplot[i], Y_endplot[i]), lw = 1.5)
            plt.gca().add_line(line)
            
else:
    line = plt.Line2D((0, xmax), (ymax,ymax))
    plt.gca().add_line(line)
    line = plt.Line2D((xmax, xmax), (0, ymax))
    plt.gca().add_line(line)
    
# Set bounds of axis
plt.axis('square')
ax.set_xlim(0,max(xpos+ypos)+15)
ax.set_ylim(0,max(xpos+ypos)+15)
# Place unit number at each scatter point
numbers = list(range(1,len(xpos)+1))
for i,txt in enumerate(numbers):
    ax.annotate(txt, (xpos[i]-.5,ypos[i]-.5))
# Value of objective function as title of plot
ax.set_title(('Cost = '+str(round(value(layout.objective)))),loc='center')

# Build rectangles with respective unit sizes around units
# Create dictionary for coordinates of bottom left of rectangle of each unit
xrectdict = dict.fromkeys(units)
yrectdict = dict.fromkeys(units)
# Subtract half the dimension based on orientation and store into dictionary
for i in units:
    xrectdict[i] = x[i].varValue - l[i].varValue/2
    yrectdict[i] = y[i].varValue - d[i].varValue/2
# Convert dictionary to array
xrect, yrect = [],[]
for i in units:
    xrect.append(xrectdict[i])
    yrect.append(yrectdict[i])
# Extract length and depth data of units
length,depth = [],[]
for i in units:
    length.append(l[i].varValue)
    depth.append(d[i].varValue)
# Plot rectangles
for i in range(len(numbers)):
    rect = mpatch.Rectangle((xrect[i], yrect[i]), length[i], depth[i], fill=None, edgecolor="black")
    ax.add_patch(rect)

    
    
#%%
