#!/usr/bin/env python3
import random
import math
import numpy as np
# variables
numrep=5
dt=0.1
beta=1.0  # tester pour 300K ou alors beta = 6.67 ou 1.67
nkill=1
# initial points
xi=1.
yi=1.

######## FUNCTIONS ########
# x=pos[0] | y=pos[1]
# V(x,y)=0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*exp(-x**2-(y-1.0/3)**2) - 3*exp(-x**2-(y-5.0/3)**2) - 5*exp(-(x-1)**2-y**2) - 5*exp(-(x+1)**2-y**2

def V(x,y):
	return 0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*exp(-x**2-(y-1.0/3)**2) - 3*exp(-x**2-(y-5.0/3)**2) - 5*exp(-(x-1)**2-y**2) - 5*exp(-(x+1)**2-y**2;

# reaction coordinate: returns it for the last point
def reccoord(ind):
	rep[ind][-1][0]=x;
	rep[ind][-1][1]=y;
	d = np.sqrt( (x**2 + (y-1)**2 );
	return 1/d ;# put here the reaction coordinate equation

# returns the zone where the last point of the replica is
def zone(ind):
	# val=-1 if in A
	# rep[ind][-1] est la dernière position de la particule d'indice ind
	rep[ind][-1][0]=x;
	rep[ind][-1][1]=y;
	
	#si on est en A
	if V(x,y)<2 and np.sqrt( (x**2 + (y+1)**2 ) <= 0.25:
		val = -1;
	if V(x,y)<2 and np.sqrt( (x**2 + (y-1)**2 ) <= 0.25:
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
    return [0.8*x**4-6*x*math.exp(-x**2-(y-1./3)**2) +
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
			level[i]=where

# for AMS: we can introduce the loop here!
	# step 1: define the killing level
	# step 2: check the stop criterion
	# step 3: kill and replicate
	#
	# SOME IDEAS:
	# introduce a variable for the killing level
	# introduce a list with indices of replicas to kill
	# introduce a list with indices of alive replicas
	# introduce a function that replicates a replica (make a hard copy until the first point after the killing level and then evolve using the run function)
	# This last function will change the replica, the size and the level lists!
	# rajouter le nombre de répliques qui rentrent dans B
	
	# TO PRINT:
	# killing level elovution
	# the final replicas
	# the final probability

# for DNS: really easy!
	# Evolve until zone != 0
	# count number of times reached A or B
	# restart from the initial point
	# Note that some parts of the program are the same, but there is no need to keep the replicas in memory!
	#
	# TO PRINT:
	# the replicas that reached B
	# the final probability
	# number of killed replicas

print(rep)
print(size)

