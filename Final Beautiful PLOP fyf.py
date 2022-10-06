# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:49:20 2021

@author: Nicholas Goh and Isobel Melvin
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


Continuous-domain representation of the process plant layout optimisation
problem with Dow F&EI constraints

Illustrative example adapted from D. I. PATSIATZIS et. al., 2004

Supplementary equations from Dow's F&EI Guide

Authors: Rey Tan, Chris Tighe 2020
"""
#%% --------------Library set-up--------------
# Import functions
from pulp import *
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatch
# from scipy.stats import norm
import math
import csv
import time

# Create the model object
layout = pulp.LpProblem("Layout_Problem_Model",LpMinimize)

#%% --------------Errors--------------
class LandAreaError(Exception):
    """Exception raised when Land Feasible region does not tally with grid size definition."""
pass #pass here as class has no definition so avoids error 

class LandShapeError(Exception):
    """Exception raised when multiple land shape definitions are selected"""
pass #pass here as class has no definition so avoids error 

class ConvexShapeError(Exception):
    """Exception raised when the 'tail' of non-convex shape is placed out of bounds"""
pass #pass here as class has no definition so avoids error 
#%%--------------Switches--------------
# CPLEX free edition only supports up to 1000 variables. For larger land sizes, use CBC or increase coarseness (g)
# 1 = CPLEX, 2 = CBC
solver = 1

# Toggle constraints in layout problem (1 is on; 0 is off)

# Land Shape Constraints (1 for convex non-square polygon, 0 for rectangles or square)
SwitchLandShape = 0

# Land Shape Constraints (1 for Non-convex shape (T or L))
SwitchNonConvex = 1

# Land Cost constraints (1 is on; 0 is off)
SwitchLandCost = 0

# Fixed aspect ratio grids or variable aspect ratio grids
SwitchAspRatio = 0

# Toggle Minimum Separation Distances switch
SwitchMinSepDistance = 0

# FEI Constraints
SwitchFEI = 1

# FEI Cost of Life Constraints (only works if SwitchFEI is on)
SwitchFEIVle = 1

#%% --------------Define Sets, Parameters and other Input Specifications --------------
#%% 
### 1. DEFINE THE PROCESS UNITS
units =  ['furnace', 'reactor', 'flash','comp','distil', 'store','ctrlroom'] #3distill comlmns in total
pertinent_units = ['furnace', 'reactor', 'distil', 'store']

# Dimensions of each unit (m)
alpha = dict.fromkeys(units)
alpha['furnace'] = 7
alpha['reactor'] = 3.6 + 2 + 2 #8.6
alpha['flash'] = 4
alpha['comp'] = 14
alpha['distil'] = 8
alpha['store'] = 22 
alpha['ctrlroom'] = 15

beta = dict.fromkeys(units)
beta['furnace'] = 16
beta['reactor'] = 1.8 + 2 + 2 #5.8 #reactor diameter + 2m allowance for diameter of shell.
beta['flash'] = 4
beta['comp'] = 14
beta['distil'] = 8
beta['store'] = 12
beta['ctrlroom'] = 15

# Purchase cost of each unit (dollars)
Cp = dict.fromkeys(units)
Cp['furnace'] = 778820
Cp['reactor'] = 1243620
Cp['flash'] = 763840
Cp['comp'] = 3426640
Cp['distil'] = 4308640
Cp['store'] = 5151768     
Cp['ctrlroom'] = 1227010 

### 2. DEFINE THE SHAPE OF PLANT
if SwitchNonConvex == 0:
    # For right-angled quadrilaterals definitions only 
    if SwitchLandShape == 0:
        # Dimensions of rectangular plot area
        xmax = 150
        ymax = 150
        
    #For all convex polygons
    elif SwitchLandShape == 1: 
        # Coordinates of vertices of convex polygon.
        X_coord = [0, 0, 25, 50, 50, 0]
        Y_coord = [0, 50, 70, 50, 0, 0]
 
# Inputs for L-shape or T-shape plot areas.   
elif SwitchNonConvex == 1:
    d_nc = 50 # The width/height of a square (m)
    l_nc = 20 # The length of 'tail' from square (m)
    t_nc = 20# The thickness of the 'tail'
    h_nc = 20# The height at the bottom of the tail from y=0

### 3. DEFINE PIPE CONNECTIONS AND PIPE FLOW CONTENTS
# Flow velocity in pipe [m/s] (3 for liquids and 30 for gas)
velocity = np.zeros((len(units), len(units)))
velocity[0][1] = 30  # velocity between furnace and reactor
velocity[1][2] = 30  # velocity between reactor and flash
velocity[2][3] = 30  # velocity between flash and comp
velocity[2][4] = 3  # velocity between flash and distil
velocity[3][0] = 30  # velocity between comp and furnace
velocity[4][5] = 3  # velocity between distil and store
velocity[5][0] = 3  # velocity between store and furnace
velocity[4][0] = 3 # velocity between distil and furnace

velocity = velocity + velocity.T - np.diag(velocity.diagonal())#
velocity = makeDict([units,units],velocity,0)

## Mass flowrate in pipe [kg/s]
Q = np.zeros((len(units), len(units)))
Q[0][1] = 15.0778  # flowrate between furnace and reactor
Q[1][2] = 15.0778   # flowrate between reactor and flash
Q[2][3] = 6.2015 # flowrate between flash and comp
Q[2][4] = 8.87631  # flowrate between flash and distil
Q[3][0] = 4.34102  # flowrate between comp and furnace
Q[4][5] = 5.94936  # flowrate between distil and store
Q[5][0] = 7.66106  # flowrate between store and furnace
Q[4][0] = 2.44706      # flowrate between distil and furnace

Q = Q + Q.T - np.diag(Q.diagonal())
Q = makeDict([units,units],Q,0)

# Flow viscosity in pipe [Pa.s]
visc = np.zeros((len(units), len(units)))
visc[0][1] = 0.0000244021  # viscosity between furnace and reactor
visc[1][2] = 0.0000210 # viscosity between reactor and flash
visc[2][3] = 0.0000115938 # viscosity between flash and comp
visc[2][4] = 0.00039666  # viscosity between flash and distil
visc[3][0] = 0.0000122  # viscosity between comp and furnace
visc[4][5] = 0.000592162597018692  # viscosity between distil and store
visc[5][0] = 0.0005699812  # viscosity between store and furnace
visc[4][0] = 0.000238829053918871   # viscosity between distil and furnace

visc = visc + visc.T - np.diag(visc.diagonal())
visc = makeDict([units,units],visc,0)

# Density of fluid in pipe [kg/cum]
rhog = np.zeros((len(units), len(units)))
rhog[0][1] = 8.01272  # density between furnace and reactor
rhog[1][2] = 10.2485  # density between reactor and flash
rhog[2][3] = 10.0488  # density between flash and comp
rhog[2][4] = 842.6054  # density between flash and distil
rhog[3][0] = 11.1785  # density between comp and furnace
rhog[4][5] = 871.310417768016  # density between distil and store
rhog[5][0] = 866.8986 # density between store and furnace
rhog[4][0] = 778.577117350617  # density between distil and furnace

rhog = rhog + rhog.T - np.diag(rhog.diagonal())
rhog = makeDict([units,units],rhog,0)

# Number of pipes between units
npp = np.zeros((len(units), len(units)))
npp[0][1] = 1  # number of pipes between furnace and reactor
npp[1][2] = 1  # number of pipes between reactor and flash
npp[2][3] = 1  # number of pipes between flash and comp
npp[2][4] = 1  # number of pipes between flash and distil
npp[3][0] = 1  # number of pipes between comp and furnace
npp[4][5] = 1  # number of pipes between distil and store
npp[5][0] = 1  # number of pipes between store and furnace
npp[4][0] = 1  # number of pipes between distil and furnace

npp = npp + npp.T - np.diag(npp.diagonal())
npp = makeDict([units,units],npp,0)

### 4. DEFINE SAFETY-RELATED PARAMETERS
# Inputs for minimum separation distance between units (m)
Demin = np.zeros((len(units), len(units)))
Demin[0][1] = 15.24 #50 ft
Demin[0][2] = 15.24
Demin[0][3] = 15.24
Demin[0][4] = 15.24
Demin[0][5] = 15.24
Demin[0][6] = 60.96 #200 ft
Demin[1][2] = 4.572 #15 feet
Demin[1][3] = 4.572
Demin[1][4] = 4.572
Demin[1][5] = 15.24
Demin[1][6] = 60.96
Demin[2][3] = 4.572
Demin[2][4] = 4.572
Demin[2][5] = 15.24
Demin[2][6] = 60.96
Demin[3][4] = 4.572
Demin[3][5] = 15.24
Demin[3][6] = 60.96
Demin[4][5] = 15.24
Demin[4][6] = 60.96
Demin[5][6] = 60.96

Demin = Demin + Demin.T - np.diag(Demin.diagonal())
Demin = makeDict([units,units],Demin,0)

# Inputs for F&EI and MPPD calculations
if SwitchFEI == 1:
    # Exposure radius of pertinent units due to fire and explosion (m)
    De = dict.fromkeys(pertinent_units)
    De['furnace'] =  28.8 * np.sqrt(2)
    De['reactor'] = 40.4 * np.sqrt(2)
    De['distil'] = 35.4* np.sqrt(2)
    De['store'] = 45.0* np.sqrt(2)

    # Damage factor of pertinent units
    DF = dict.fromkeys(pertinent_units)
    DF['furnace'] =  0.75
    DF['reactor'] = 0.82
    DF['distil'] = 0.77
    DF['store'] = 0.82
    
    # Frequency/probability of event occuring at a pertinent unit (/year)
    freq_event = dict.fromkeys(units)
    freq_event['furnace'] = 1e-3
    freq_event['reactor'] = 1e-3
    freq_event['distil']  = 1e-3
    freq_event['store']   = 1e-3
    
    if SwitchFEIVle == 1:
        # Number of workers
        Nw = dict.fromkeys(units)
        Nw['furnace'] = 2
        Nw['reactor'] = 3
        Nw['flash'] = 2
        Nw['comp'] = 3
        Nw['distil'] = 4 
        Nw['store'] = 2
        Nw['ctrlroom'] = 10

        # Percentage of time orkers are present at unit
        tw = dict.fromkeys(units)
        tw['furnace'] = 0.1
        tw['reactor'] = 0.1
        tw['flash'] = 0.1
        tw['comp'] = 0.1
        tw['distil'] = 0.1
        tw['store'] = 0.1
        tw['ctrlroom'] = 0.375
        
### 5. DEFINE DEGREE OF FINENESS FOR AREA CALCULATION
if SwitchLandCost == 1:
    if SwitchAspRatio == 0:
        # For fixed aspect ratios - N determines the dimensions of each individual grid relative to xmax and ymax
        N = 80
        
    elif SwitchAspRatio == 1: 
        # For free aspect ratios - N1 determines the number of grids in x-direction, N2 determines number of grids in y-directions.
        #Number of grids in each direction
        N1 = 80
        N2 = 80
        #Length and width of one grid square: each grid square is g * g in dimension
        #!! Ensure that g*N1 = xmax, g*N2 = ymax
        g = xmax/N1
    
#%% --------------Default Parameters---------------------------
# M (m) is used in the contraints which set A or B, L or R, and Din or Dout.
# It should be big enough not to constrain the size of the plot, but not too big
M = 1e3

# Default parameters for Connection Costs
C_ref = 7000  # Reference installed cost of duct of $700 per metre
n_1 = 1.08 # Parameter for installed cost of ducts 
n_2 = 0.74 # Paramter for type of material 
CEPCI_1959  = 100 # Accounts for rise of inflation 1959
CEPCI_2006 = 499.6 # Accounts for rise of inflation 2006
CEPCI_2021 = 677.7 # Accounts for rise of inflation 2021
MF = 1 # Material factor (different for each component)
A_f = 0.1102  #Annualised factor - multiplied to Land cost, Road cost, Piping cost and Pumping cost
FX_rate = 0.75  # USD to GBP exchange rate as of 1 Dec 2021
BB = 880  # Parameter dependent on material
F = 1.5 # Parameter dependent on material 
DIA_ref = 2.3 # Diameter of pipe
epsilon = 0.0046 # Absolute roughness of steel pipe/dict [m]. check chin chern paper for reference.
mechEffic = 0.6  # Mechanical efficiency 
C_elec = 0.000045 # Wholesale cost of electircity 
OH = 8000  # Operating hours

# Annualised land cost per squared distance (£/m^2/year)
LC = 61.78 * A_f 

# Annualised cost of building paths/roads (£/m^2/year)
RC = 35 * A_f 
# Width of Road
road_width = 4

# Cost of life (dollars)
Cl = dict.fromkeys(units)
for i in units:
    Cl[i] = 10e6

#%% Error checking
if SwitchLandShape + SwitchNonConvex >= 2:
    raise LandShapeError("Please select just one land shape definition")

if SwitchNonConvex == 1:
    if h_nc + t_nc > d_nc :
        raise ConvexShapeError("Please check d_nc >= t_nc + h_nc")
    
if SwitchAspRatio == 1:
    #Check for whether N1*g and N2*g is = xmax,ymax
    if N1*g != xmax or N2*g != ymax:
        raise LandAreaError("Please check N1, N2 and g value to ensure N1*g = xmax, N2*g = ymax")
    
    
#%% --------------PRE-CONSTRAINTS SETUP-----------------------
### 1. CONVEX POLYGON FORMATION
if SwitchLandShape == 1:
    X_begin = np.array(X_coord[:-1])
    X_end = np.array(X_coord[1:])
    X_endlim = max(X_end)
    XDiff = X_end - X_begin
    
    Y_begin = np.array(Y_coord[:-1])
    Y_end = np.array(Y_coord[1:])
    YDiff = Y_end - Y_begin
    
    # Definitions for helping line  plotting:
    X_beginplot = list(X_begin)
    X_endplot = list(X_end)
    Y_beginplot = list(Y_begin)
    Y_endplot = list(Y_end)
    
    # Check for infinite gradients
    grad0_list = list(reversed(list(np.where(XDiff==0) [0])))
    
    # Compensates for inability to plot infinite gradients
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
    
    # Construction of line constraints
    XDiff = X_end - X_begin
    YDiff = Y_end - Y_begin      
    grad = YDiff/XDiff
    absgrad = abs(grad)
    y_intercept = Y_begin - (grad * X_begin)
    
    # Converting to list for manipulation.
    X_begin = list(X_begin)
    X_end = list(X_end)
    Y_begin = list(Y_begin)
    Y_end = list(Y_end)
    grad = list(grad)
    absgrad = list(absgrad)
    y_intercept = list(y_intercept)

### 2. OCCUPANCY CALCULATIONS
OCC = dict.fromkeys(units)
OCC = {i: Nw[i]*tw[i] for i in units}

### 3. LAND COST SETUP
if SwitchLandCost == 1:    
    if SwitchAspRatio == 0:
        # Length and width of one grid square (m)
        g_x = xmax/N
        g_y = ymax/N       
        gridsize = list(range(1,N))
        
    elif SwitchAspRatio == 1:
        # Defining set for binary variable Gn and Gn1n2
        gridsizen1 = list(range(1,N1+1))
        gridsizen2 = list(range(1,N2+1))
        
### 4. CONNECTION COST SETUP
# Initialising dictionaries for pressure drop calculations
DIA = np.zeros((len(units), len(units)))
DIA = makeDict([units,units],DIA,0)
C = np.zeros((len(units), len(units)))
C = makeDict([units,units],C,0)
PC = np.zeros((len(units), len(units)))
PC = makeDict([units,units],PC,0)
C_annual = np.zeros((len(units), len(units)))
C_annual = makeDict([units,units],C_annual,0)

# Initialiseing dictionary for pumping costs
kr = np.zeros((len(units), len(units)))
kr = makeDict([units,units],kr,0)
Rey = np.zeros((len(units), len(units)))
Rey = makeDict([units,units],Rey,0)
AA = np.zeros((len(units), len(units)))
AA = makeDict([units,units],AA,0)
ff = np.zeros((len(units), len(units)))
ff = makeDict([units,units],ff,0)
PC_unit = np.zeros((len(units), len(units))) #Pipe cost per unit length
PC_unit = makeDict([units,units],PC_unit,0)
delP_unit = np.zeros((len(units), len(units))) #pressure drop per unit length
delP_unit = makeDict([units,units],delP_unit,0)
pPump_unit = np.zeros((len(units), len(units))) #pump power per unit length
pPump_unit = makeDict([units,units],pPump_unit,0)
pipetotal_unit = np.zeros((len(units), len(units))) #pump power per unit length
pipetotal_unit = makeDict([units,units],pipetotal_unit,0)

# Calculations for pipe/pumping costs
for idxj, j in enumerate(units):
    for idxi, i in enumerate(units):
        if idxj > idxi:
            #Equations for piping cost
            DIA[i][j] = np.sqrt( (4*Q[i][j] / (velocity[i][j] * np.pi * rhog[i][j])))
            C[i][j] = C_ref*(DIA[i][j]/DIA_ref)**n_1 * (CEPCI_2021/CEPCI_1959) * MF * FX_rate * A_f
            PC[i][j] = BB * (DIA[i][j]**n_2) * (CEPCI_2021/CEPCI_2006) *FX_rate
            C_annual[i][j] = PC[i][j] * (1+F) * (A_f)
            
            #Equations for pumping cost
            kr[i][j] = epsilon/ DIA[i][j]
            Rey[i][j] = rhog[i][j] * velocity[i][j] * DIA[i][j] / visc[i][j]
            AA[i][j] = (kr[i][j] ** 1.1098)/2.8257 + (7.149 / Rey[i][j])**0.8981
            ff[i][j] = (1/ (-2 * math.log(kr[i][j]/3.7065 - (5.0452/Rey[i][j])*math.log(AA[i][j]))) )**2
            delP_unit[i][j] = 8 * ff[i][j] * rhog[i][j] * (velocity[i][j] **2) / (2 * DIA[i][j])
            pPump_unit[i][j] = Q[i][j] * delP_unit[i][j] /(rhog[i][j]*mechEffic)
            PC_unit[i][j] = C_elec * OH * pPump_unit[i][j]*npp[i][j]
            
            # Annualised Piping/Pumping cost (£/m/year)
            pipetotal_unit[i][j] = PC_unit[i][j]+ C_annual[i][j] *npp[i][j]


#%%---------------------CONSTRAINTS VARIABLE DEFINITIONS----------------
# 1 if length of item i is equal to alpha; 0 otherwise
O = LpVariable.dicts("O",(units),lowBound=0,upBound=1,cat="Integer")
# pair-wise values to ensure equipment items do not overlap
E1 = LpVariable.dicts("E1",(units,units),lowBound=0,upBound=1,cat="Integer")
E2 = LpVariable.dicts("E2",(units,units),lowBound=0,upBound=1,cat="Integer")
# Binary viriables determining whether unit i to the right/left or above/below unit j
Wx = LpVariable.dicts("Wx",(units,units),lowBound=0,upBound=1,cat="Integer")
Wy = LpVariable.dicts("Wy",(units,units),lowBound=0,upBound=1,cat="Integer")

# Binary variable for non-convex shapes 
G1 = LpVariable.dicts("G1",(units),lowBound=0,upBound=1,cat="Integer")

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

# Total rectilinear distance between items i and j
D = LpVariable.dicts("D",(units,units),lowBound=0,upBound=None,cat="Continuous")
    # coordinates of the geometrical centre of item i
x = LpVariable.dicts("X",(units),lowBound=None,upBound=None,cat="Continuous")
y = LpVariable.dicts("Y",(units),lowBound=None,upBound=None,cat="Continuous")

# Define connection cost
CD = LpVariable.dicts("CD",(units,units),lowBound=0,upBound=None,cat="Continuous")
# Total connection cost
SumCD = LpVariable("SumCD",lowBound=0,upBound=None,cat="Continuous")

# Define road construction cost
roadCost = LpVariable("roadCost",lowBound=0,upBound=None,cat="Continuous")
D_road =  LpVariable.dicts("D_road",(units, units),lowBound=0,upBound=None,cat="Continuous")


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

if SwitchLandCost == 1:
    # N binary variables representing plot grid
    if SwitchAspRatio == 1: 
        Gn1n2 = LpVariable.dicts("Gn1n2", (gridsizen1, gridsizen2), lowBound = 0, upBound = 1, cat = "Integer")# For g_x[i] and g_y[i], select the larger as g    
    else:
        Gn = LpVariable.dicts("Gn",(gridsize), lowBound=0, upBound=1, cat="Integer")
    # Total Land Cost
    TLC = LpVariable("TLC",lowBound=0,upBound=None,cat="Continuous")

elif SwitchLandCost == 0:
    TLC = 0


#%% --------------DEFINE OBJECTIVE FUNCTION--------------
obj_sumOmega = SwitchFEI*SumOmega
obj_CD = SumCD
obj_TLC = SwitchLandCost *TLC

layout += obj_CD + obj_sumOmega + obj_TLC

#%% --------------DEFINE CONSTRAINTS AND OBJECITVE FUNCTION CONTRIBUTION--------------
### 1. BASE MODEL CONSTRAINTS FOR ALL UNITS i
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
            # Distance calculation constraints (3 - 10)
            layout += R[i][j] - L[i][j] == x[i] - x[j]
            layout += A[i][j] - B[i][j] == y[i] - y[j]
            layout += R[i][j] <= M*Wx[i][j]
            layout += L[i][j] <= M*(1 - Wx[i][j])
            layout += A[i][j] <= M*Wy[i][j]
            layout += B[i][j] <= M*(1 - Wy[i][j])
            layout += D[i][j] == R[i][j] + L[i][j] + A[i][j] + B[i][j]
            layout += D[i][j] == D[j][i]
            
            # These constraints ensure consistency in interdependent variables
            layout += L[i][j] == R[j][i]
            layout += R[i][j] == L[j][i]
            layout += A[i][j] == B[j][i]
            layout += B[i][j] == A[j][i]
            layout += Wx[i][j] == 1 - Wx[j][i]
            layout += Wy[i][j] == 1 - Wy[j][i]
            
            # Defines cost of pumping and piping
            layout += CD[i][j] == (pipetotal_unit[i][j]) * D[i][j] * npp[i][j]

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

# Defines cost of road from control room to other units.
for idxj, j in enumerate(units):
    if idxj<6:
        layout += D_road["ctrlroom"][j] == D["ctrlroom"][j]   
# Cost of Road
layout += roadCost ==  RC * road_width * lpSum([D_road["ctrlroom"][k] for k in units])

# Objective function contribution for base model
layout += SumCD == (roadCost + lpSum([CD[i][j] for i in units for j in units]))

### 2. LAND SHAPE CONSTRAINTS
if SwitchLandShape == 1:
    for i in units:
        for k in y_intercept:
           if k > 0:
               layout += y[i] + 0.5*d[i] + l[i]*absgrad[y_intercept.index(k)]/2 <= grad[y_intercept.index(k)]*x[i] + k
           elif k <= 0:
               layout += y[i] - 0.5*d[i] - l[i]*absgrad[y_intercept.index(k)]/2 >= grad[y_intercept.index(k)]*x[i] + k
        
        layout += x[i] + 0.5*l[i] <= X_endlim
        layout += x[i] -0.5*l[i] >= 0
              
elif SwitchLandShape == 0:
    if SwitchNonConvex == 0:
        for i in units:
            layout += x[i] + 0.5*l[i] <= xmax
            layout += y[i] + 0.5*d[i] <= ymax
    
    elif SwitchNonConvex == 1:
        for i in units:
            layout += y[i] + 0.5*d[i] <= d_nc*G1[i] + (h_nc + t_nc) * (1-G1[i])
            layout += x[i] + 0.5*l[i] <= d_nc*G1[i] + (d_nc + l_nc) *(1-G1[i])
            layout += y[i] - 0.5*d[i] >= h_nc*(1-G1[i])
            layout += x[i] - 0.5*l[i] >= 0
        
    
### 3. LAND COST CONSTRAINTS
if SwitchLandShape == 0 and SwitchNonConvex == 0:
    if SwitchLandCost == 1:
        if SwitchAspRatio == 1:    
            for i in units:              
                # Variable aspect ratio up to limits N1 and N2 * g            
                layout += x[i] + 0.5 * l[i] <= lpSum(n1 * g * Gn1n2[n1][n2] for n1 in range(1,N1) for n2 in range(1, N2))
                layout += y[i] + 0.5 * d[i] <= lpSum(n2 * g * Gn1n2[n1][n2] for n1 in range(1,N1) for n2 in range(1, N2))
                                                             
                # Only 1 grid size selected (rectangle)
                layout += lpSum((Gn1n2[n1][n2]) for n1 in range(1,N1) for n2 in range(1,N2)) == 1
    
                # Objective function contribution for variable aspect ratio land use model
                layout += TLC == LC*lpSum((Gn1n2[n1][n2] * n1*g * n2*g) for n1 in range(1,N1) for n2 in range(1,N2))    
        
        elif SwitchAspRatio == 0:
            for i in units:
                # Fixed aspect ratio land area approximation constraints for plot
                layout += x[i] + 0.5*l[i] <= lpSum(n*g_x*Gn[n] for n in range(1,N))
                layout += y[i] + 0.5*d[i] <= lpSum(n*g_y*Gn[n] for n in range(1,N))
            
                # Only 1 grid size selected (fixed aspect ratio)
                layout+= lpSum(Gn[n] for n in range(1,N)) == 1
            
                # Objective function contribution for fixed aspect ratio land use model
                layout += TLC == LC*lpSum(Gn[n]*n*g_x*n*g_y for n in range(1,N))
    
    
### 4. F&EI CONSTRAINTS
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

                    if SwitchFEIVle == 1:
                        layout += Ve2[i][j] == Cp[j]*Psie[i][j] + OCC[j]*Cl[j]*Psie[i][j] - Cp[j]*Din[i][j]/De[i] - OCC[j]*Cl[j]*Din[i][j]/De[i]
                    else:
                        layout += Ve2[i][j] == Cp[j]*Psie[i][j] - Cp[j]*Din[i][j]/De[i]
                else:
                    layout += Ve2[i][j] == 0
    # for all i in pertinent units
    for i in pertinent_units:
        # Value of area of exposure constraint (30ii)
        if SwitchFEIVle == 1: # Cost of replacing unit that explodes is not included, since its just a constant and doesnt actually tell us anything useful.
            layout += Ve[i] == freq_event[i] * (Cp[i] + OCC[i]*Cl[i] + lpSum([Ve2[i][j] for j in units]))
        else:
            layout += Ve[i] == freq_event[i]*(Cp[i] + lpSum([Ve2[i][j] for j in units]))

        # Base maximum probable property damage cost (31)
        layout += Omega0[i] == DF[i]*Ve[i]


    # F&EI objective function contribution
    layout += SumOmega == lpSum(Omega0[i] for i in pertinent_units)


### 5. FIXING BLOCKS TO CERTAIN AREAS
def fix_variable(variable, value):
    variable.setInitialValue(value)
    variable.fixValue()

x_store = alpha["store"]/2
y_store = beta["store"]/2

solution_x = {('store'): x_store}
solution_y = {('store'): y_store}

for i, v in solution_x.items():
    x[i].setInitialValue(v)

for i, v in solution_y.items():
    y[i].setInitialValue(v)

for i, v in solution_x.items():
    fix_variable(x[i], v)

for i, v in solution_y.items():
    fix_variable(y[i], v)
    
#%% -------------- INITIATE SOLVE --------------
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

if SwitchLandCost == 1:
    if SwitchAspRatio == 1:
        for n1 in range(1,N1):
            for n2 in range(1,N2):
                if Gn1n2[n1][n2].varValue ==1:
                    xfinalaxis = n1
                    yfinalaxis = n2                    
        print("Number of grids", xfinalaxis, "x", yfinalaxis, ". Size of land area =", (xfinalaxis*yfinalaxis*g*g), "metres square")
    else:
        for n in range(1, N):
            if Gn[n].varValue == 1:
                xyfinalaxis = n
                print("Number of grids", n, "Size of land area =", (n*n*g_x*g_y), "metres square")

if SwitchFEI == 1:
    print("Total actual MPPD =", SumOmega.varValue)
else:
    print("Total actual MPPD =", SumOmega.varValue)

    
if SwitchLandCost == 1:
    print("Total Land Use =", TLC.varValue)

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

# Plotting Convex Polygons
if SwitchLandShape == 1 and SwitchNonConvex == 0:
    for k in X_begin:
        line = plt.Line2D((k, X_end[X_begin.index(k)]), (Y_begin[X_begin.index(k)], Y_end[X_begin.index(k)]), lw = 1.5)
        plt.gca().add_line(line)
    if len(grad0_list)>0:  
        for i in grad0_list:
            line = plt.Line2D((X_beginplot[i], X_endplot[i]), (Y_beginplot[i], Y_endplot[i]), lw = 1.5)
            plt.gca().add_line(line)
            
# Plotting right angled quadrilaterals
if SwitchLandShape == 0 and SwitchNonConvex == 0:
    line = plt.Line2D((0, xmax), (ymax,ymax))
    plt.gca().add_line(line)
    line = plt.Line2D((xmax, xmax), (0, ymax))
    plt.gca().add_line(line)
    
#Plotting Non-convex shapes:
if SwitchNonConvex == 1:
    
    line1 = plt.Line2D((0,d_nc),(d_nc,d_nc))
    plt.gca().add_line(line1)
    line2 = plt.Line2D((d_nc,d_nc),(d_nc,h_nc + t_nc))
    plt.gca().add_line(line2)
    line3 = plt.Line2D((d_nc,d_nc + l_nc),(h_nc + t_nc,h_nc + t_nc))
    plt.gca().add_line(line3)
    line4 = plt.Line2D((d_nc + l_nc,d_nc + l_nc),(h_nc + t_nc,h_nc))
    plt.gca().add_line(line4)
    line5 = plt.Line2D((d_nc + l_nc,d_nc),(h_nc,h_nc))
    plt.gca().add_line(line5)
    line6 = plt.Line2D((d_nc,d_nc),(h_nc,0))
    plt.gca().add_line(line6)
    
# Set bounds of axis
plt.axis('square')

# For Non Convex Shapes
if SwitchNonConvex == 1:
    ax.set_xlim(0,d_nc+l_nc+10)
    ax.set_ylim(0,d_nc+10)

# For convex polygons:
if SwitchLandShape == 1:
    ax.set_xlim(0,max(X_beginplot))
    ax.set_ylim(0,max(Y_beginplot))

# For right angled quadrilateral shapes
if SwitchLandShape == 0 and SwitchNonConvex == 0:
    if SwitchLandCost == 1:
        # For Variable Aspect ratios:
            if SwitchAspRatio == 1:
                ax.set_xlim(0,xfinalaxis*g)
                ax.set_ylim(0,yfinalaxis*g)      
        # For Fixed aspect ratios:
            elif SwitchAspRatio == 0:
                ax.set_xlim(0,xyfinalaxis*g_x)
                ax.set_ylim(0,xyfinalaxis*g_y)     
    # Display if no land cost considered
    elif SwitchLandCost == 0:   
        ax.set_xlim(0,max(xpos)+0.5*l[units[xpos.index(max(xpos))]].varValue)
        ax.set_ylim(0,max(ypos)+0.5*d[units[ypos.index(max(ypos))]].varValue)   
    
# Place unit number at each scatter point
numbers = list(range(1,len(xpos)+1))
sizes = [3,3,3,3,3,3]
for i,txt in enumerate(numbers):
    ax.annotate(txt, (xpos[i]-.5,ypos[i]-.5), fontsize = 8)
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

    
    