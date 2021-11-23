# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 13:28:37 2021

@author: Nicholas Goh
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
# from scipy.stats import norm
import math
import csv
import time

# Create the model object
layout = pulp.LpProblem("Layout_Problem_Model",LpMinimize)

#%% --------------Errors--------------
# class NoFEIError(Exception):
#     """Exception raised when protection devices enabled without FEI constraints"""
# pass #pass here as class has no definition so avoids error 

#%%--------------Switches--------------
# CPLEX free edition only supports up to 1000 variables. For larger land sizes, use CBC or increase coarseness (g)
# 1 = CPLEX, 2 = CBC
solver = 1

# Toggle constraints in layout problem (1 is on; 0 is off)

# Land Shape Constraints (1 for non-square polygon, 0 for square based on xmax and ymax)
SwitchLandShape = 0

# Land Cost constraints (1 is on; 0 is off)
SwitchLandUse = 0

# Toggle Minimum Separation Distances switch
SwitchMinSepDistance = 0

# FEI Constraints
SwitchFEI = 0

# Toggle protection devices (must have FEI enabled if 1)
SwitchProt = 0

# FEI Cost of Life Constraints (only works if SwitchFEI is on)
SwitchFEIVle = 1

# CEI Constraints
SwitchCEI = 0

#%% Case selection:
Casee = 1

if Casee == 1:
     SwitchFEI = 1
     SwitchLandUse = 0
elif Casee == 2:    
    SwitchMinSepDistance = 1
    SwitchLandUse = 1
    SwitchFEI = 1
elif Casee == 3:
    SwitchLandShape = 1
    SwitchFEI = 1
    

# Check for errors if SwitchProt == 1 and SwitchFEI == 0:
#    raise NoFEIError("FEI cost of life constraints not allowed without FEI constraints")        
if SwitchProt == 1 and SwitchFEI == 0:
    raise NoFEIError("Protection devices not allowed without FEI constraints") 

#%% --------------Define Sets--------------
# Define the process units
# furnace = furnace, reactor = FEHE + reactor, Flash = C1 + flash, Comp = Comp, Distil = RECCOL + STAB + PRODCOL + C2 + P2, Store = storage tanks + P1
units =  ['furnace', 'reactor', 'flash','comp','distil', 'store'] #3distill comlmns in total
pertinent_units = ['furnace', 'reactor', 'distil', 'store']
hazardous_chemicals = ['tol','benz','meth','h2','diph']
Nunits = len(units)


#%% --------------Define Parameters and Values--------------
# Base layout model

# M (m) is used in the contraints which set A or B, L or R, and Din or Dout.
# It should be big enough not to constrain the size of the plot, but not too big ##for use of big M method?
M = 1e3

#polygon layout:
X_begin = np.array([0,0,25,50,50])
X_end = np.array([0,25,50,50,0])
Ng = max(X_end)
XDiff = X_end - X_begin

Y_begin = np.array([0,40,60,40,0])
Y_end = np.array([40,60,40,0,0])
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

#%%----------------- INPUT SPECIFICATIONS-------------------------
# Dimensions of each unit (m)
alpha = dict.fromkeys(units)
beta = dict.fromkeys(units)
alpha['furnace'] = 1
alpha['reactor'] = 9 #to be checked reactor itself is vertical column, 10 feet or 3m. HEX has got 6.1m length, approximately diameter = 
alpha['flash'] = 2.3
alpha['comp'] = 12
alpha['distil'] = 5.5
alpha['store'] = 24 # Squeezing them into a rectangle


beta['furnace'] = 1
beta['reactor'] = 5 #reactor diameter + 2m allowance for diameter of shell.
beta['flash'] = 2.3
beta['comp'] = 12
beta['distil'] = 5.5
beta['store'] = 33


# Purchase cost of each unit (dollars)
Cp = dict.fromkeys(units)
Cp['furnace'] = 556300
Cp['reactor'] = 888300
Cp['flash'] = 545600
Cp['comp'] = 2447600
Cp['distil'] = 3077600
Cp['store'] = 1890074 # to be checked

#Minimum separation distances
Demin = np.zeros((len(units), len(units)))
Demin[0][1] = 15
Demin[0][2] = 15
Demin[0][3] = 15
Demin[0][4] = 15
Demin[0][5] = 15
Demin[1][2] = 15
Demin[1][3] = 15
Demin[1][4] = 15
Demin[1][5] = 15
Demin[2][3] = 15
Demin[2][4] = 15
Demin[2][5] = 15
Demin[3][4] = 15
Demin[3][5] = 15
Demin[4][5] = 15

Demin = Demin + Demin.T - np.diag(Demin.diagonal())
Demin = makeDict([units,units],Demin,0)

## define the velocities [m/s]
velocity = np.zeros((len(units), len(units)))
velocity[0][1] = 30  # velocity between furnace and reactor
velocity[1][2] = 30  # velocity between reactor and flash
velocity[2][3] = 30  # velocity between flash and comp
velocity[2][4] = 3  # velocity between flash and distil
velocity[3][0] = 30  # velocity between comp and furnace
velocity[4][5] = 30  # velocity between distil and store
velocity[5][0] = 3  # velocity between store and furnace

velocity = velocity + velocity.T - np.diag(velocity.diagonal())#
velocity = makeDict([units,units],velocity,0)

## define the flowrates Q [kg/s]
Q = np.zeros((len(units), len(units)))
Q[0][1] = 15.0778  # flowrate between furnace and reactor
Q[1][2] = 15.0778   # flowrate between reactor and flash
Q[2][3] = 4.34105  # flowrate between flash and comp
Q[2][4] = 8.87631  # flowrate between flash and distil
Q[3][0] = 15.0778  # flowrate between comp and furnace
Q[4][5] = 5.94936  # flowrate between distil and store
Q[5][0] = 7.66106  # flowrate between store and furnace

Q = Q + Q.T - np.diag(Q.diagonal())
Q = makeDict([units,units],Q,0)

## define viscosity #STREAM TABLES [Pa.s]
visc = np.zeros((len(units), len(units)))
visc[0][1] = 0.0000244021  # viscosity between furnace and reactor
visc[1][2] = 0.0000115938  # viscosity between reactor and flash
visc[2][3] = 0.0000115938 # viscosity between flash and comp
visc[2][4] = 0.00039666  # viscosity between flash and distil
visc[3][0] = 0.000017046  # viscosity between comp and furnace
visc[4][5] = 0.00000900998  # viscosity between distil and store
visc[5][0] = 0.000587848  # viscosity between store and furnace

visc = visc + visc.T - np.diag(visc.diagonal())
visc = makeDict([units,units],visc,0)

#density of the stream fluid #from aspen [kg/cum]
rhog = np.zeros((len(units), len(units)))
rhog[0][1] = 8.01272  # density between furnace and reactor
rhog[1][2] = 10.0488  # density between reactor and flash
rhog[2][3] = 10.0488  # density between flash and comp
rhog[2][4] = 842.605  # density between flash and distil
rhog[3][0] = 13.5901  # density between comp and furnace
rhog[4][5] = 2.73253  # density between distil and store
rhog[5][0] = 869.273  # density between store and furnace

rhog = rhog + rhog.T - np.diag(rhog.diagonal())
rhog = makeDict([units,units],rhog,0)

# number of pipes 
npp = np.zeros((len(units), len(units)))
npp[0][1] = 1  # number of pipes between furnace and reactor
npp[1][2] = 1  # number of pipes between reactor and flash
npp[2][3] = 1  # number of pipes between flash and comp
npp[2][4] = 1  # number of pipes between flash and distil
npp[3][0] = 1  # number of pipes between comp and furnace
npp[4][5] = 1  # number of pipes between distil and store
npp[5][0] = 1  # number of pipes between store and furnace

npp = npp + npp.T - np.diag(npp.diagonal())
npp = makeDict([units,units],npp,0)

#define constants for piping cost 

C_ref = 7000  # reference installed cost of duct of $700 per metre
n_1 = 1.08# parameter for installed cost of ducts 
n_2 = 0.74 # paramter for type of material 
CEPCI_1959  = 100 # accounts for rise of inflation 1959
CEPCI_2006 = 499.6 #accounts for rise of inflation 2006
CEPCI_2021 = 677.7   #accounts for rise of inflation 2021
MF = 1 # material factor (different fro each component)
A_f = 0.1102  # check report
FX_rate = 0.64  # exchange rate
BB = 880  # paramter dependent on material
F = 1.5 # parameter dependent on material 
# bb = 0.05 # 5% of installed cost accounting fro maintanance 
DIA_ref = 2.3 # diameter of pipe
epsilon = 0.0046 #absolute roughness of steel pipe/dict [m]. check chin chern paper for reference.
mechEffic = 0.6  #mechanical efficiency 
C_elec = 0.000045 # wholesale cost of electircity 
OH = 8000  # operating hours

#%% Land shape constraint: 1 if non-rectangular, 0 if rectangular.
if SwitchLandShape == 0: #sets default max available plot area.
    xmax = 60
    ymax = 60
    #%% Land use constraint: 
    if SwitchLandUse == 1:
        # Land cost per squared distance (m^2)
        LC = 125
        # Number of grid points in square plot
        N = 60
        # Length and width of one grid square (m)
        g_x = xmax/N
        g_y = ymax/N
        # g = xmax/N
        # Defining set for binary variable Gn
        gridsize = list(range(1,N))
        

#%% ----------- SwitchFEI---------
if SwitchFEI == 1:
    # F&EI factors
    De = dict.fromkeys(pertinent_units)
    DF = dict.fromkeys(pertinent_units)

    De['furnace'] =  28.8
    De['reactor'] = 23.9
    De['distil'] = 35.4
    De['store'] = 45.0
    # De['furnace'] =  1
    # De['reactor'] = 1
    # De['distil'] = 11
    # De['store'] = 2

    DF['furnace'] =  0.75
    DF['reactor'] = 0.66
    DF['distil'] = 0.77
    DF['store'] = 0.82
    
    # DF['furnace'] =  0.5
    # DF['reactor'] = 0.4
    # DF['distil'] = 0.4
    # DF['store'] = 0.4


    # Upper bound for actual maximum probable property damage cost
    U = 1e8
    
# #%%--------------- SwitchProt-----------------------------
#     # Protection device model
#     if SwitchProt == 1:
#         # Define protection device configuration
#         configurations = list(range(1, len(units)))
#         # Loss control credit factor of protection device configuration k on item i
#         gamma = np.zeros((len(pertinent_units), len(configurations)))
#         gamma = makeDict([pertinent_units,configurations],gamma,0)
#         # assign values
#         gamma['reactor'][1] = 1
#         gamma['reactor'][2] = 0.900
#         gamma['reactor'][3] = 0.750
#         gamma['reactor'][4] = 0.365
#         gamma['reactor'][5] = 0.292
#         gamma['reactor'][6] = 0.117
#         gamma['eoabs'][1] = 1
#         gamma['eoabs'][2] = 0.900
#         gamma['eoabs'][3] = 0.760
#         gamma['eoabs'][4] = 0.684
#         gamma['eoabs'][5] = 0.612
#         gamma['eoabs'][6] = 0.465
#         gamma['co2abs'][1] = 1
#         gamma['co2abs'][2] = 0.900
#         gamma['co2abs'][3] = 0.760
#         gamma['co2abs'][4] = 0.684
#         gamma['co2abs'][5] = 0.612
#         gamma['co2abs'][6] = 0.465
#         # purchase cost of configuration k on unit i
#         P = np.zeros((len(pertinent_units),len(configurations)))
#         P = makeDict([pertinent_units,configurations],P,0)
#         # assign values
#         P['reactor'][1] = 0
#         P['reactor'][2] = 5000
#         P['reactor'][3] = 15000
#         P['reactor'][4] = 40000
#         P['reactor'][5] = 60000
#         P['reactor'][6] = 125000
#         P['eoabs'][1] = 0
#         P['eoabs'][2] = 5000
#         P['eoabs'][3] = 20000
#         P['eoabs'][4] = 25000
#         P['eoabs'][5] = 35000
#         P['eoabs'][6] = 55000
#         P['co2abs'][1] = 0
#         P['co2abs'][2] = 5000
#         P['co2abs'][3] = 20000
#         P['co2abs'][4] = 25000
#         P['co2abs'][5] = 35000
#         P['co2abs'][6] = 55000

#%%--------------------CEI Factors parameters---------------------------
# # # Initialise dictionaries
# AQ = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced
# AQf = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced by flash
# AQp = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # airborne quantity produced by pool
# Fv = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}  # fraction flashed
# Pvap = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
# CEI = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
# Dc = {i: dict.fromkeys(hazardous_chemicals) for i in pertinent_units}
# Dh = dict.fromkeys(pertinent_units)
# Liq = dict.fromkeys(pertinent_units)
# Deltah = dict.fromkeys(pertinent_units)
# WT = dict.fromkeys(pertinent_units)
# WP = dict.fromkeys(pertinent_units)
# AP = dict.fromkeys(pertinent_units)
# heatratio = dict.fromkeys(hazardous_chemicals)
# Tb = dict.fromkeys(hazardous_chemicals)
# antoineA = dict.fromkeys(hazardous_chemicals)
# antoineB = dict.fromkeys(hazardous_chemicals)
# antoineC = dict.fromkeys(hazardous_chemicals)
# MW = dict.fromkeys(hazardous_chemicals)
# Pr = dict.fromkeys(hazardous_chemicals)
# prob_death = dict.fromkeys(hazardous_chemicals)
# CONC = dict.fromkeys(hazardous_chemicals)
# # # Assign values
# te = 10  # minutes of exposure time
# prob_death['EO'] = 0.5  # probability of death for probit
# Dh['reactor'] = 50.8  # mm
# Dh['eoabs'] = 50.8  # hole diameter
# Dh['co2abs'] = 50.8
# Deltah['reactor'] = 1  # m
# Deltah['eoabs'] = 1
# Deltah['co2abs'] = 1
# heatratio['EO'] = 0.00365  # 1/degc
# Tb['EO'] = 10.5  # degc
# antoineA['EO'] = 4.386
# antoineB['EO'] = 1115.1
# antoineC['EO'] = -29.015
# MW['EO'] = 44.05

#%%----------------- SwitchFEIVle -----------------------
# # Occupancy calculations
# # Number of workers
# if SwitchFEIVle == 1:
Nw = dict.fromkeys(units)
Nw['furnace'] = 2
Nw['reactor'] = 3
Nw['flash'] = 2
Nw['comp'] = 3
Nw['distil'] = 4 
Nw['store'] = 2
# Percentage of time present at unit
tw = dict.fromkeys(units)
tw['furnace'] = 0.1
tw['reactor'] = 0.1
tw['flash'] = 0.1
tw['comp'] = 0.1
tw['distil'] = 0.1
tw['store'] = 0.1

# Occupancy
OCC = dict.fromkeys(units)
OCC = {i: Nw[i]*tw[i] for i in units}

# Cost of life (dollars)
Cl = dict.fromkeys(units)
for i in units:
    Cl[i] = 1e6

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

#initialise dictionary for pumping costs
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

for idxj, j in enumerate(units):
    for idxi, i in enumerate(units):
        if idxj > idxi:
            #Equations for piping cost
            DIA[i][j] = np.sqrt( (4*Q[i][j] / (velocity[i][j] * np.pi * rhog[i][j])))
            C[i][j] = C_ref*(DIA[i][j]/DIA_ref)**n_1 * (CEPCI_2021/CEPCI_1959) * MF * FX_rate * A_f
            PC[i][j] = BB * (DIA[i][j]**n_2) * (CEPCI_2021/CEPCI_2006) *FX_rate
            # bb = 0.05 * Cp[i]
            C_annual[i][j] = PC[i][j] * (1+F) * (A_f)
            
            #Equations for pumping cost
            kr[i][j] = epsilon/ DIA[i][j]
            Rey[i][j] = rhog[i][j] * velocity[i][j] * DIA[i][j] / visc[i][j]
            AA[i][j] = (kr[i][j] ** 1.1098)/2.8257 + (7.149 / Rey[i][j])**0.8981
            ff[i][j] = (1/ (-2 * math.log(kr[i][j]/3.7065 - (5.0452/Rey[i][j])*math.log(AA[i][j]))) )**2
            delP_unit[i][j] = 8 * ff[i][j] * rhog[i][j] * (velocity[i][j] **2) / (2 * DIA[i][j])
            pPump_unit[i][j] = Q[i][j] * delP_unit[i][j] /(rhog[i][j]*mechEffic)
            PC_unit[i][j] = C_elec * OH * pPump_unit[i][j]*npp[i][j]
            
            pipetotal_unit[i][j] = PC_unit[i][j]+ C_annual[i][j] *npp[i][j]

if SwitchLandUse == 1:
    # N binary variables representing plot grid
    Gn = LpVariable.dicts("Gn",(gridsize), lowBound=0, upBound=1, cat="Integer")
    # For g_x[i] and g_y[i], select the larger as g    
    # Total Land Cost
    TLC = LpVariable("TLC",lowBound=0,upBound=None,cat="Continuous")
else:
    TLC = 0
    
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


#%% --------------Define Objective Function--------------

obj_sumOmega = SwitchFEI*SumOmega
obj_PZ = SwitchProt*SumPZ
obj_CD = SumCD
obj_TLC = SwitchLandUse *TLC
if Casee == 2:
    layout += obj_CD + obj_TLC
else:
    layout += obj_CD + obj_sumOmega + obj_PZ + obj_TLC


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
            layout += CD[i][j] == (pipetotal_unit[i][j]) * D[i][j] * npp[i][j]
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
    if SwitchLandUse == 1:
        for i in units:
            # Land area approximation constraints for plot
            layout += x[i] + 0.5*l[i] <= lpSum(n*g_x*Gn[n] for n in range(1,N))
            layout += y[i] + 0.5*d[i] <= lpSum(n*g_y*Gn[n] for n in range(1,N))
        
            # Only 1 grid size selected
            layout+= lpSum(Gn[n] for n in range(1,N)) == 1
        
            # Objective function contribution for land use model
            layout += TLC == LC*lpSum(Gn[n]*n*g_x*n*g_y for n in range(1,N))
    else:
        for i in units:
            layout += x[i] + 0.5*l[i] <= xmax
            layout += y[i] + 0.5*d[i] <= ymax
        
#%% Land use constraints
# if SwitchLandUse == 1:
#     # Land area approximation constraints for plot
#     layout += x[i] + 0.5*l[i] <= lpSum(n*g*Gn[n] for n in range(1,N))
#     layout += y[i] + 0.5*d[i] <= lpSum(n*g*Gn[n] for n in range(1,N))
    
#     # Only 1 grid size selected
#     layout+= lpSum(Gn[n] for n in range(1,N)) == 1
    
#     # Objective function contribution for land use model
#     layout += TLC == LC*lpSum(Gn[n]*(n*g)**2 for n in range(1,N))    
    
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

        # Base maximum probable property damage cost (31)
        layout += Omega0[i] == DF[i]*Ve[i]

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


#%% Fixing blocks to certain areas
def fix_variable(variable, value):
    variable.setInitialValue(value)
    variable.fixValue()

x_store = alpha["store"]/2
y_store = beta["store"]/2

# solution_x = {('store'): x_store,
#               ('furnace'): 90,
#               ('reactor'): 9,
#               ('comp'): 88}
# solution_y = {('store'): y_store,
#               ('furnace'):90,
#               ('reactor'):90,
#               ('comp'):12}

solution_x = {('store'): x_store}
              # ('comp'): 45}
solution_y = {('store'): y_store}
              # ('comp'): 37}

for i, v in solution_x.items():
    x[i].setInitialValue(v)

for i, v in solution_y.items():
    y[i].setInitialValue(v)

for i, v in solution_x.items():
    fix_variable(x[i], v)

for i, v in solution_y.items():
    fix_variable(y[i], v)
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

if SwitchLandUse == 1:
    for n in range(1, N):
        if Gn[n].varValue == 1:
            print("Number of grids", n, "Size of land area =", (n*n*g_x*g_y), "metres square")

if SwitchFEI == 1:
    print("Total actual MPPD =", SumOmega.varValue)
else:
    print("Total actual MPPD =", SumOmega.varValue)

    
if SwitchLandUse == 1:
    print("Total Land Use =", TLC.varValue)
# if SwitchProt == 1:
#     print("Total cost of protection devices =", SumPZ.varValue)
#    if SwitchFEIVle == 1:
#        print("Total cost of fatality due to fire/explosion =", SumVle.varValue)
# if SwitchCEI == 1:
#     print("Total cost of fatality due to chemical release =", SumVlc.varValue)
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
# ax.set_xlim(0,max(xpos+ypos)+15)
# ax.set_ylim(0,max(xpos+ypos)+15)
if SwitchLandShape == 1:
    ax.set_xlim(0,max(X_beginplot))
    ax.set_ylim(0,max(Y_beginplot))
else:
    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)

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

    
    
#%% Data processing
## Actual area used
#Collect all vertices of units
# xvertices = []
# for i in units:
#     xvertices.append(x[units].varValue + x[units].varValue)
#     xvertices.append(x[units].varValue - x[units].varValue)
#     yvertices.append(y[units].varValue + y[units]).varValue)
#     yvertices.append(y[units].varValue + y[units].varValue)
#create convex hull

#find area of convex hull