#!/usr/bin/env python3
import random
import math
import numpy as np
import copy

# variables
numrep=5
dt=0.1
beta=1.0
nkill=1
# initial points
xi=1.
yi=1.

######## FUNCTIONS ########
# x=pos[0] | y=pos[1]
# V(x,y)=0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*exp(-x**2-(y-1.0/3)**2) - 3*exp(-x**2-(y-5.0/3)**2) - 5*exp(-(x-1)**2-y**2) - 5*exp(-(x+1)**2-y**2
def V(x,y):
	return 0.2*x**4 + 0.2*(y-1.0/3)**4 + 3*np.exp(-x**2-(y-1.0/3)**2) - 3*np.exp(-x**2-(y-5.0/3)**2) - 5*np.exp(-(x-1)**2-y**2) - 5*np.exp(-(x+1)**2-y**2)
	
# reaction coordinate: returns it for the last point
# def reccoord(ind):
# 	return 1 # put here the reaction coordinate equation
def reccoord(ind,endroit):
	rep[ind][endroit][0]=x;
	rep[ind][endroit][1]=y;
	d = np.sqrt( (x-1)**2 + y**2 );
	return 1/d # put here the reaction coordinate equation

# returns the zone where the last point of the replica is
# def zone(ind):
# 	# val=-1 if in A
# 	# val=1 if in B
# 	# val=0 if otherwise
# 	return val
def zone(ind):
	val=0
	# val=-1 if in A
	# rep[ind][-1] est la dernière position de la particule d'indice ind
	rep[ind][-1][0]=x;
	rep[ind][-1][1]=y;
	
	#si on est en A
	if V(x,y)<2 and np.sqrt( (x**2 + (y+1)**2 )) <= 0.25:
		val = -1;
	if V(x,y)<2 and np.sqrt( (x**2 + (y-1)**2 )) <= 0.25:
		val = 1;
	else :
		val=0;
	# val=1 if in B
	# val=0 if otherwise
	# on definit A et B comme des zones ou le niveau est < quelque chose et proche d'un certain point (pour ne pas confondre A et B)
	return val

# potential derivate
def dPot(pos):
	return [1,1] # put here the potential gradient

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
	rep[ind].append([rep[ind][size[ind]][0]-dt*tPot[0]+tgx,rep[ind][size[ind]][1]-dt*tPot[1]+tgy]) #rep(ind) est la liste avec tous les points, rep(ind)(ind)(ind) ça nous donne soit x soit y
	size[ind]=size[ind]+1
	return reccoord(ind,-1)

##################################################################
# building the replicas structure
rep=[]
size=[]
level=[]
for i in range(numrep):
	rep.append([[xi,yi]])
	size.append(0)
	level.append(reccoord(i))

# initialization of the replicas ; where stocke la valeur de la coordonné de réaction
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
	#
	# TO PRINT:
	# killing level evolution
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
# 
# print(rep)
# print(size)

#faire une deepcopy lorsqu'on copie le début d'une trajectoire
#à chaque itération, nombre de répliques tuées
#dans un fichier séparé mettre les résultats pour faire des probas
#variable nombre de répliques dans B
#1.67 et 6.67 => calculer les valeurs des températures correspondantes

# AMS
#function that replicates a replica
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
	killing_level=copy.deepcopy(level).sort()[k-1]
	if killing_level>=zmax:
		break;
	#on remplit les listes des replicas à tuer ou à conserver
	for l in range(level):
		if level[l]<killing_level:
			replicas_to_kill.append(l)
		else :
			alive_replicas.append(l)
	if len(replicas_to_kill)==numrep :
		break;
	for indice in replicas_to_kill :
		replication(indice,alive_replicas, killing_level)
			
	
	
	


	
	
	
	
	
	
	
	