#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt

# variables
numrep=5
dt=0.1
beta=1.0  # tester pour 300K ou alors beta = 6.67 ou 1.67
nkill=1
# initial points
xi=1.
yi=1.

#%%####### FUNCTIONS ########
# x=pos[0] | y=pos[1]
# V(x,y)=0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*exp(-x**2-(y-1.0/3)**2) - 3*exp(-x**2-(y-5.0/3)**2) - 5*exp(-(x-1)**2-y**2) - 5*exp(-(x+1)**2-y**2

def V(x,y):
	return 0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*np.exp(-x**2-(y-1.0/3)**2) - 3*np.exp(-x**2-(y-5.0/3)**2) - 5*np.exp(-(x-1)**2-y**2) - 5*np.exp(-(x+1)**2-y**2) +4;

# reaction coordinate: returns it for the last point
def reccoord(ind, endroit):
	d=0;
	x=rep[ind][endroit][0];
	y=rep[ind][endroit][1];
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
	if V(x,y)<0.5 and np.sqrt( (x+1)**2 + y**2 ) <= 0.25:
		val = -1;
	elif V(x,y)<0.5 and np.sqrt( (x-1)**2 + y**2 ) <= 0.25:
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
	return reccoord(ind,-1)

#%%#################################################################
# building the replicas structure
rep=[]
size=[]
level=[]
for i in range(numrep):
	rep.append([[xi,yi]])
	size.append(0)
	level.append(reccoord(i,-1))

# initialization of the replicas
for i in range(numrep):
	while zone(i) == 0:
		where=run(i)
		if where > level[i]:
			level[i]=where
			
def replication(indice,alive_replicas, killing_level): 
	#indice pointe dans rep au replica qui doit être supprimé
	#killing_level est la valeur en-dessous de laquelle on supprime les réplicas
	#ind est l'indice du replica à copier
	#alive_replicas est une liste des replicas 
	global rep
	global level
	rep[indice]=[]
	ind=random.randint(0,len(alive_replicas)-1)
	ind=alive_replicas[ind]
	level[indice]=killing_level-1 #initialisé pour être inférieur à killing_level
	k=0
	while level[indice]<=killing_level:
		rep[indice].append(rep[ind][k])
		level[indice]=reccoord(ind,k)
		k+=1
	rep[indice].append(rep[ind][k])#point dont la reccoord est > à killing_level
	#on met à jour level
	level[indice]=reccoord(ind,k)
	#on poursuit la trajectoire de la particule correspondant au nouveau replica
	while zone(indice) == 0:
		where=run(indice)
		if where > level[indice]:
			level[i]=where

#AMS
zmax=4.255
while(True):
	alive_replicas=[]
	replicas_to_kill=[]
	liste_triee=copy.deepcopy(level);
	liste_triee.sort()
	killing_level=liste_triee[nkill-1]
	
	if killing_level>=zmax:
		break;
	#on remplit les listes des replicas à tuer ou à conserver
	for l in range(len(level)):
		if level[l]<killing_level:
			replicas_to_kill.append(l)
		else :
			alive_replicas.append(l)
	if len(replicas_to_kill)==numrep :
		break;
	for indice in replicas_to_kill :
		replication(indice,alive_replicas, killing_level)
        
        
#%% Courbes

X = np.linspace(-1.5,1.5,100)
Y = np.linspace(-1.,2.,100)

X,Y = np.meshgrid(X,Y)
Z = V(X,Y) 

plt.figure()

plt.contour(X,Y,Z) 

liste_X = [[] for i in range(numrep)]
liste_Y = [[] for i in range(numrep)]

for i in range(numrep):
    for j in range(len(rep[i])):
        liste_X[i].append(rep[i][j][0])
        liste_Y[i].append(rep[i][j][1])


for i in range(numrep):
    plt.plot(liste_X[i],liste_Y[i])

plt.show()
