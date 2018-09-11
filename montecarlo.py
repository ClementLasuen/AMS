#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import random
import math
import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages

# variables
numrep=1
dt=0.01
beta=3.67  # tester pour 300K ou alors beta = 6.67 ou 1.67
nkill=1
# initial points
xi=-0.7
yi=0.

######## FUNCTIONS ########
# x=pos[0] | y=pos[1]
# V(x,y)=0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*exp(-x**2-(y-1.0/3)**2) - 3*exp(-x**2-(y-5.0/3)**2) - 5*exp(-(x-1)**2-y**2) - 5*exp(-(x+1)**2-y**2

def V(x,y):
    m= 4 + 0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*np.exp(-x**2-(y-1.0/3)**2) - 3*np.exp(-x**2-(y-5.0/3)**2) - 5*np.exp(-(x-1)**2-y**2) - 5*np.exp(-(x+1)**2-y**2);
    return(m)

# reaction coordinate: returns it for the last point
def reccoord(ind):
    x=rep[ind][-1][0];
    y=rep[ind][-1][1];
    d = math.sqrt( (x-1)**2 + y**2 );
    return 1/d ;# put here the reaction coordinate equation

# returns the zone where the last point of the replica is
def zone(ind):
    # val=-1 if in A
    # rep[ind][-1] est la dernière position de la particule d'indice ind
    val=0
    x=rep[ind][-1][0];
    y=rep[ind][-1][1];
	
    #si on est en A
    if V(x,y)<0.5 and math.sqrt( (x+1)**2 + y**2 ) <= 0.25:
        val = -1;
    elif V(x,y)<0.5 and math.sqrt( (x-1)**2 + y**2 ) <= 0.25:
        val = 1;
    else :
        val=0;
    # val=1 if in B
    # val=0 if otherwise
    # on definit A et B comme des zones ou le niveau est < quelque chose et proche d'un certain point (pour ne pas confondre A et B)
    return val

# potential derivate
def dPot(pos):
    x=pos[0]
    y=pos[1]
    return [0.8*x**3-6*x*math.exp(-x**2-(y-1./3)**2) +
            6*x*math.exp(-x**2-(y-5./3)**2) +
            10*(x-1)*math.exp(-(x-1)**2-y**2) +
            10*(x+1)*math.exp(-(x+1)**2-y**2) ,
            0.8*(y-1./3)**3-6*(y-1./3)*math.exp(-x**2-(y-1./3)**2) +
            6*(y-5./3)*math.exp(-x**2-(y-5./3)**2) +
            10*y*math.exp(-(x-1)**2-y**2) +
            10*y*math.exp(-(x+1)**2-y**2)] # put here the potential gradient

# evolution function
def run(ind):
    global rep
    global size
    # evolution for replica ind, return reaction coordinate
    u=random.random()
    v=random.random()
    tgx=math.sqrt(2.*dt/beta)*math.sqrt(-2*math.log(u))*math.cos(2*math.acos(-1)*v)
    tgy=math.sqrt(2.*dt/beta)*math.sqrt(-2*math.log(u))*math.sin(2*math.acos(-1)*v)
    # ref: https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
    tPot=dPot(rep[ind][size[ind]])
    rep[ind].append([rep[ind][size[ind]][0]-dt*tPot[0]+tgx,rep[ind][size[ind]][1]-dt*tPot[1]+tgy])
    size[ind]=size[ind]+1
    return reccoord(ind)

##################################################################
# building the replicas structure
    
#%%Montecarlo

num_test=10000
montecarlo = 0.  
modulo=100

X = np.linspace(-1.5,1.5,100)
Y = np.linspace(-1.,2.,100)

X,Y = np.meshgrid(X,Y)

Z = V(X,Y)

LISTE_X=[]
LISTE_Y=[]

for j in range(num_test):
    if j%modulo==0:
        print(j)
        sys.stdout.flush()

    rep=[]
    size=[]
    level=[]
    for i in range(numrep):
        rep.append([[xi,yi]])
        size.append(0)
        level.append(reccoord(i))

    # initialization of the replicas

    
    for i in range(numrep):
        while zone(i) == 0:
            where=run(i)
            if where > level[i]:
                level[i] = where
                
    if j%modulo==0:
        liste_x=[rep[0][k][0] for k in range(len(rep[0]))]
        liste_y=[rep[0][k][1] for k in range(len(rep[0]))]
        LISTE_X.append(liste_x)
        LISTE_Y.append(liste_y)
    if zone(i)==1:
        montecarlo+=1
        
montecarlo/=num_test
print(montecarlo) # Proba finale


# Sauvegarde de la figure au format pdf

filename = 'FIGURE_MONTECARLO_BETA-{}_dt-{}_NUMTEST-{}_NUMTRACES-{}.pdf'.format(beta,dt,num_test,int(num_test/modulo))
with PdfPages(filename) as pdf :
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title(r'Répliques, $\beta$ =' + str(beta) + r' , $dt$ = ' + str(dt) + r' , proba = ' + str(montecarlo))
    for i in range(len(LISTE_X)):
        ax.plot(LISTE_X[i],LISTE_Y[i])
    plt.contour(X,Y,Z)
    pdf.savefig()
    plt.close()      
