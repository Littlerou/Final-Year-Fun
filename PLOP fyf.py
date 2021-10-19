## testing git comment
<<<<<<< HEAD
#dgdrgdg 


#dfgdf 
#suckmadeek python ya mong  


=======
#testing testing
>>>>>>> 5bd8217b87a7c26d062826bd23f2542124f657cb

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:21:50 2021

@author: Nicholas Goh

PLOP Edited for new input

**IMPORTANT:
1) Ctrl+F and search the keyword 'UNIT' to look for places where data needs to be inputted.
2) Ctrl+F and search 'REVIEW' to look for sections where code needs to be understood/removed/edited.


"""
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

# --------------Set-up--------------
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

# --------------Errors--------------
class NoFEIError(Exception):
    """Exception raised when protection devices enabled without FEI constraints"""
pass

# --------------Switches--------------
# CPLEX free edition only supports up to 1000 variables. For larger land sizes, use CBC or increase coarseness (g)
# 1 = CPLEX, 2 = CBC
solver = 1

# Toggle constraints in layout problem (1 is on; 0 is off)
# Land Use Constraints
SwitchLandUse = 1
# FEI Constraints
SwitchFEI = 1
# Toggle protection devices (must have FEI enabled if 1)
SwitchProt = 1
# FEI Cost of Life Constraints (only works if SwitchFEI is on)
SwitchFEIVle = 1
# CEI Constraints
SwitchCEI = 1

# Check for errors if SwitchProt == 1 and SwitchFEI == 0:
#    raise NoFEIError("FEI cost of life constraints not allowed without FEI constraints")        
if SwitchProt == 1 and SwitchFEI == 0:
    raise NoFEIError("Protection devices not allowed without FEI constraints") 

# --------------Define Sets--------------
# Define the process units
units = []
pertinent_units = []
hazardous_chemicals = []
Nunits = len(units)

# --------------Define Parameters and Values--------------
# Base layout model

# M (m) is used in the contraints which set A or B, L or R, and Din or Dout.
# It should be big enough not to constrain the size of the plot, but not too big
M = 1e3

# Dimensions of each unit (m)
alpha = dict.fromkeys(units)
beta = dict.fromkeys(units)
alpha['UNIT'] = 

beta['UNIT'] =

# Purchase cost of each unit (dollars)
Cp = dict.fromkeys(units)
Cp['UNIT'] = 

# Connection/piping costs (cost per unit length)
C = np.zeros((len(units), len(units)))  # square matrix with elements unit,unit
# assign values where flowsheet connection
C['UNIT']['UNIT2'] =  # connection cost between UNIT and UNIT

# symmetrise the matrix
C = C + C.T - np.diag(C.diagonal())
# Cost data is made into a dictionary
C = makeDict([units,units],C,0)

# Land use model
if SwitchLandUse == 1:
    # Land cost (per unit distance squared)
    LC = 3e3
    # Number of grid points in square plot
    N = 250
    # Size of one grid square (m)
    g = 2
    # Define the set for binary variable Gn
    gridsize = list(range(1,N))
else:
    # Maximum plot size if land use model not switched on
    xmax = 1000
    ymax = 1000

# F&EI model
if SwitchFEI == 1:
    # Operating condidions
    T = dict.fromkeys(pertinent_units)
    Pg = dict.fromkeys(pertinent_units)
    w = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
    rhol = dict.fromkeys(pertinent_units)
    MWunit = dict.fromkeys(pertinent_units)
    INV = dict.fromkeys(pertinent_units)
    # Assign values
    T['UNIT'] = 275  # Temperature, degc
    
    Pg['reactor'] = 2200  #Pressure (gauge???), kPa

    w['UNIT']['UNITCONTENT'] = 0.94  # mass frac

    rhol['UNIT'] =       #to identify what rhol is.
    
    MWunit['UNIT'] = 44.05  # Molecular weight, g/mol

    INV['UNIT'] = 10000  # kg


'REVIEW'
    # F&EI factors
#    MF = dict.fromkeys(pertinent_units)
#    F1 = dict.fromkeys(pertinent_units)
#    F2 = dict.fromkeys(pertinent_units)
#    F = dict.fromkeys(pertinent_units)
    De = dict.fromkeys(pertinent_units)
    DF = dict.fromkeys(pertinent_units)
#    # assign values
#    MF['UNIT']

#    F1['UNIT'] = 

#    F2['UNIT'] = 

    # compute factors
#    F3 = {i: F1[i]*F2[i] for i in pertinent_units}
#    F = {i: MF[i]*F3[i] for i in pertinent_units}
#    De = {i: 0.256*F[i] for i in pertinent_units}

    DF['UNIT'] = 
    
    De['UNIT'] =

'REVIEW'
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

    # Upper bound for actual maximum probable property damage cost  'REVIEW'
    U = 1e8 # To be edited maybe

    # Protection device model
    if SwitchProt == 1:
        # Define protection device configuration
        configurations = list(range(1, len(units)))
        # Loss control credit factor of protection device configuration k on item i
        gamma = np.zeros((len(pertinent_units), len(configurations)))
        gamma = makeDict([pertinent_units,configurations],gamma,0)
        # assign values
        gamma['UNIT'][1] = 1
        gamma['UNIT'][2] = 0.900
        etc.

        gamma['UNIT2'][1] = 1
        gamma['UNIT2'][2] = 0.900
        etc.

        # purchase cost of configuration k on unit i
        P = np.zeros((len(pertinent_units),len(configurations)))
        P = makeDict([pertinent_units,configurations],P,0)
        # assign values
        P['UNIT'][1] = 0
        P['UNIT'][2] = 5000

        P['UNIT2'][1] = 0
        P['UNIT2'][2] = 5000

# ## CEI Factors
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
te = 10  # minutes of exposure time 'REVIEW'
prob_death['EO'] = 0.5  # probability of death for probit 'REVIEW'
Dh['UNIT'] = 50.8  # hole diameter, mm

Deltah['UNIT'] = 1  # m

heatratio['UNIT'] = 0.00365  # 1/degc

Tb['UNITCONTENT'] = 10.5  # degc

antoineA['UNITCONTENT'] = 4.386

antoineB['UNITCONTENT'] = 1115.1

antoineC['UNITCONTENT'] = -29.015

MW['UNITCONTENT'] = 44.05

# # Compute CEI factors ' BIG REVIEW'
#     # Probit equation  'REVIEW' CHECK probit equation linearity
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


# # Occupancy calculations
# # Number of workers
if SwitchFEIVle == 1:
    Nw = dict.fromkeys(units)
    Nw['UNIT'] = 3

    # Percentage of time present at unit
    tw = dict.fromkeys(units)
    tw['UNIT'] = 0.10

    # Occupancy
    OCC = dict.fromkeys(units)
    OCC = {i: Nw[i]*tw[i] for i in units}
    # Cost of life (dollars) 'REVIEW'
    Cl = dict.fromkeys(units)
    for i in units:
        Cl[i] = 10e6

# --------------Define Variables--------------
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

if SwitchLandUse == 1:
    # N binary variables representing plot grid
    Gn = LpVariable.dicts("Gn",(gridsize),lowBound=0,upBound=1,cat="Integer")
    # Total land cost
    TLC = LpVariable("TLC",lowBound=0,upBound=None,cat="Continuous")

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

# --------------Define Objective Function--------------  'REVIEW' - has objective function been simplified to line 484?
# if SwitchFEI == 1 and SwitchProt == 1 and SwitchLandUse == 1:
#     layout += SumCD + TLC + SumOmega + SumPZ
# elif SwitchFEI == 1 and SwitchProt == 1 and SwitchLandUse == 0:
#     layout += SumCD + SumOmega + SumPZ
# elif SwitchFEI == 1 and SwitchProt == 0 and SwitchLandUse == 0:
#     layout += SumCD + SumOmega
# elif SwitchFEI == 1 and SwitchProt == 0 and SwitchLandUse == 1:
#     layout += SumCD + TLC + SumOmega
# else:
#     layout += SumCD

# layout += SumCD + SumOmega + SumPZ + SumVlc + TLC

# commented out as it is currently under review 'REVIEW' - to be used as switches.
 obj_sumOmega = SwitchFEI*SumOmega
 obj_PZ = SwitchProt*SumPZ
 obj_Vle = SwitchFEIVle*SumVle
 obj_TLC = SwitchLandUse*TLC
 obj_Vlc = SwitchCEI*Vlc
 obj_CD = SumCD
 layout += obj_CD + obj_sumOmega + obj_PZ + obj_Vle + obj_TLC + obj_Vlc



# --------------Define Constraints and Objective Function Contributions--------------
# Base model constraints for all units i
for i in units:
    # Orientation constraints (1 - 2)
    layout += l[i] == alpha[i]*O[i] + beta[i]*(1 - O[i])
    layout += d[i] == alpha[i] + beta[i] - l[i]
    # Lower and upper bounds of coordinates (19 - 22)
    layout += x[i] >= 0.5*l[i]
    layout += y[i] >= 0.5*d[i]

for idxj, j in enumerate(units):
    for idxi, i in enumerate(units):
        if idxj > idxi:
            layout += CD[i][j] == C[i][j] * D[i][j]
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
            layout += x[i] - x[j] + M*(E1[i][j] + E2[i][j]) >= (l[i] + l[j])/2
            layout += x[j] - x[i] + M*(1 - E1[i][j] + E2[i][j]) >= (l[i] + l[j])/2
            layout += y[i] - y[j] + M*(1 + E1[i][j] - E2[i][j]) >= (d[i] + d[j])/2
            layout += y[j] - y[i] + M*(2 - E1[i][j] - E2[i][j]) >= (d[i] + d[j])/2
            # These constraints ensure consistency in interdependent variables
            layout += L[i][j] == R[j][i]
            layout += R[i][j] == L[j][i]
            layout += A[i][j] == B[j][i]
            layout += B[i][j] == A[j][i]
            layout += Wx[i][j] == 1 - Wx[j][i]
            layout += Wy[i][j] == 1 - Wy[j][i]

# Objective function contribution for base model
layout += SumCD == lpSum([CD[i][j] for i in units for j in units])

# Land use constraints (or set max plot size if not used)
if SwitchLandUse == 1:
    for i in units:
        # Land area approximation constraints for square plot (24 - 25)
        layout += x[i] + 0.5*l[i] <= lpSum(n*g*Gn[n] for n in range(1,N))
        layout += y[i] + 0.5*d[i] <= lpSum(n*g*Gn[n] for n in range(1,N))
        # for rectangular plot
        # layout += x[i] + 0.5*l[i] <= lpSum(n1*g*Gn1n2[n1][n2] for n1 in range(1,N1) for n2 in range(1,N2))
        # layout += y[i] + 0.5*d[i] <= lpSum(n2*g*Gn1n2[n1][n2] for n1 in range(1,N1) for n2 in range(1,N2))

    # Only 1 grid size selected (23)
    layout += lpSum(Gn[n] for n in range(1,N)) == 1

    # Objective function contribution for land use model
    layout += TLC == LC*lpSum(Gn[n]*(n*g)**2 for n in range(1,N))
else:
    layout += x[i] + 0.5*l[i] <= xmax
    layout += y[i] + 0.5*d[i] <= ymax

# F&EI Constraints
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
                    #Commented out as appeared to prevent Psie taking value of 1 - check later 'REVIEW'
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


# CEI Equations 'REVIEW'
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

# --------------Fixing Variable Values-------------- 'REVIEW'
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
#     ('UNIT'): 30,#.39,
#     ('UNIT2'): 15.39,
#     ('UNIT3'): 26.04,
#     ('UNIT4'): 4.24,
#     ('UNIT5'): 4.24,
#     ('UNIT6'): 5.58,
#     ('UNIT7'): 2.9,
# }
# solution_y = {
#     ('UNIT'): 30,#.39,
#     ('UNIT2'): 15.39,
#     ('UNIT3'): 4.24,
#     ('UNIT4'): 4.24,
#     ('UNIT5'): 26.54,
#     ('UNIT6'): 9.82,
#     ('UNIT7'): 9.82,
#}


# for i, v in solution_x.items():
#     x[i].setInitialValue(v)

# for i, v in solution_y.items():
#     y[i].setInitialValue(v)

# for i, v in solution_x.items():
#     fix_variable(x[i], v)

# for i, v in solution_y.items():
#     fix_variable(y[i], v)

'REVIEW'
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

# --------------Initiate Solve--------------
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

#--------------Print Results--------------
# Print variable and objective function values
for v in layout.variables():
    print(v.name, "=", v.varValue)
print("Total cost of connections =", SumCD.varValue)
if SwitchLandUse == 1:
    for n in range(1, N):
        if Gn[n].varValue == 1:
            print("Size of land area =", n*g, "metres square")
    print("Total cost of land =", TLC.varValue)
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

#--------------Export Results--------------
filename = 'Optimisation_Plot.csv'
with open(filename, 'w', newline='') as file:
    # Write objective function
    writer = csv.writer(file)
    writer.writerow(['Objective', value(layout.objective)])
    
    # Write coordinates
    #fieldnames = ['unit','x', 'y', 'l', 'd', 'De', 'Dc']
    fieldnames = ['unit','x', 'y', 'l', 'd', 'De', 'Dc']
    writer = csv.DictWriter(file, fieldnames=fieldnames)    
    writer.writeheader()
    for i in units:
        #writer.writerow({'unit': i, 'x': x[i].varValue, 'y': y[i].varValue, 'l': l[i].varValue, 'd': d[i].varValue, 'De': De.get(i), 'Dc': Dc.get(i)})
        writer.writerow({'unit': i, 'x': x[i].varValue, 'y': y[i].varValue, 'l': l[i].varValue, 'd': d[i].varValue, 'De': De.get(i), 'Dc': Dc.get(i)})

filename = 'Optimisation_Results.csv'
with open(filename, 'w', newline='') as file:
    # Write value of all variables
    writer = csv.writer(file, delimiter=',')    
    for v in layout.variables():
        writer.writerow([v.name, v.varValue])

#--------------Plot Results--------------
xpos, ypos = [], []
for i in units:
    xpos.append(x[i].varValue)
    ypos.append(y[i].varValue)

# Plot invisible scatter
fig, ax = plt.subplots()
ax.scatter(xpos,ypos,alpha=0)
# Set bounds of axis
plt.axis('square')
ax.set_xlim(0,max(xpos+ypos)+5)
ax.set_ylim(0,max(xpos+ypos)+5)
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

