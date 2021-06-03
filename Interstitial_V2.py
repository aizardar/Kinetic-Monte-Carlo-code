import numpy as np
import numpy.random as random
import copy


# SYSTEM SIZE
L = 10

LatticeF = np.chararray((L*L),itemsize = 2)

Al_number = int(0.5*(L**2))

LatticeF[:] = "Fe"

LatticeF[0:Al_number] = "Al"

np.random.shuffle(LatticeF)

Lattice = np.reshape(LatticeF,(L,L))



E_H_4Fe = 0.0
dE_Al_Fe = 0.22
Interstitial = np.chararray(((L-1),(L-1)),itemsize = 2)
Coordination = {}


for i in range(0,len(Interstitial)):
  for j in range(0,len(Interstitial)):
      
	 Coordination[i,j] = [Lattice[i][j],Lattice[i+1][j],Lattice[i][j+1],Lattice[i+1][j+1]]


Initial_position = [random.randint(len(Interstitial)-1),random.randint(len(Interstitial)-1)] # Start always here???? random.randint (41)
print "Initial position =", Initial_position  


Conf_0 = Coordination[Initial_position[0],Initial_position[1]]
print Conf_0


def NeighCount(Conf):

    N_Al = 0
    N_Fe = 0
    p = 0
    for p in Conf:
      if(p == "Fe"):
	N_Fe = N_Fe + 1
      
      elif(p == "Al"):
	N_Al =N_Al + 1 
	
    return {'N_Fe':N_Fe, 'N_Al':N_Al}

   
   
N_Fe = NeighCount(Conf_0)['N_Fe']
N_Al = NeighCount(Conf_0)['N_Al']


def E_Conf(N_Fe,N_Al):
  
  return E_H_4Fe + N_Al* dE_Al_Fe 

print E_Conf(N_Fe,N_Al)