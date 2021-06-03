import numpy as np
import numpy.random as random


# To read the parameter file
def read_file(filename):
	with open(filename) as f_in:
        	data = filter(None, (line.rstrip() for line in f_in)) # to store data line by line without empty lines
	SystemSize= float(data[0].split()[0])    # System size
	Num_H = int(data[1].split()[0])    	# Number of hydrogen atoms
	Temp = float(data[2].split()[0])	# Temperature
	E1 = float(data[3].split()[0])  	# Reference Energy (All Fe)
	E2 = float(data[4].split()[0])		# Activation Energy for jump from all Fe to all Fe Configuration
	E3 = float(data[5].split()[0])		# Energy difference between 5Fe1Al and 6Fe configuration
	E4 = float(data[6].split()[0])		#  Energy difference when an Al atom is in the 2nd nearest neighbour
	
	return SystemSize,Num_H,Temp,E1,E2,E3,E4

L,H,T,E_H_Fe,dE_H_FeFe,dE_H_Al,dE_H_Al_2 = read_file('param.dat')

Interstitial = np.chararray((L,L,L),itemsize = 2)
Coordination_NN = {}        # To store the Coordination of first nearest neighbours 
Coordination_2NN = {}		# To store the Coordination of second nearest neighbours
E_Conf = {}			# To store the Configuration energies  
events = {}			# To store the events
barrier = {}			# To store the barriers
pos = {}
read_pos = {}


# To calculate number of Al atoms 

def NeighCount(Coordination):

	N_Al = 0

	for p in Coordination:
	  if(p == "Al"):
	    N_Al = N_Al + 1
	  


	return float(N_Al)



# To calculate the barrier

def Calc_barrier(x,y,z,xp,yp,zp):

	n = read_pos[x,y,z]
	nprime = read_pos[xp,yp,zp]

	return  E_Conf[n] - E_Conf[nprime] + dE_H_FeFe + NeighCount(Coordination_NN[n])*dE_H_Al + NeighCount(Coordination_NN[nprime])*dE_H_Al
	

# To calculate the Coordinations, Energies, events, and barriers


def Calc_Coordination(Al_Concentration):

	Al_number = int((Al_Concentration*(L*L*L))/100.0)
	LatticeF = np.chararray((L*L*L),itemsize = 2)  #   Simple Cubic Superlattice
	TotSites = len(LatticeF)
	LatticeF[:] = "Fe"   #  Fe = Iron
	LatticeF[0:Al_number] = "Al"  # Al = Aluminium
	np.random.shuffle(LatticeF)
	Lattice = np.reshape(LatticeF,(L,L,L))
	
	
	n = 0 

	for x in range(len(Interstitial)):
		for y in range(len(Interstitial)):
			for z in range(len(Interstitial)):
				pos[n] = [x,y,z]
				read_pos[x,y,z] = n 
				n = n + 1 
				 
		


	for n in range(TotSites):
		
		x,y,z = pos[n]

		Coordination_NN[n] = [Lattice[np.mod(x,L)][np.mod(y-1,L)][np.mod(z,L)],Lattice[np.mod(x-1,L)][np.mod(y,L)][np.mod(z,L)],Lattice[np.mod(x,L)][np.mod(y+1,L)][np.mod(z,L)],Lattice[np.mod(x+1,L)][np.mod(y,L)][np.mod(z,L)],Lattice[np.mod(x,L)][np.mod(y,L)][np.mod(z-1,L)],Lattice[np.mod(x,L)][np.mod(y,L)][np.mod(z+1,L)]]
		
	for n in range(TotSites):
		
		x,y,z = pos[n]

		Coordination_2NN[n] = Coordination_NN[read_pos[np.mod(x+2,L),y,z]] + Coordination_NN[read_pos[np.mod(x-2,L),y,z]] + Coordination_NN[read_pos[x,np.mod(y+2,L),z]] + Coordination_NN[read_pos[x,np.mod(y-2,L),z]] + Coordination_NN[read_pos[x,y,np.mod(z+2,L)]] + Coordination_NN[read_pos[x,y,np.mod(z-2,L)]] 
	

	for n in range(TotSites):
		
		x,y,z = pos[n]

		E_Conf[n] = E_H_Fe + float(NeighCount(Coordination_NN[n]) * dE_H_Al) + float(NeighCount(Coordination_2NN[n]) * dE_H_Al_2)
		


	
	for n in range(TotSites):
		
		x,y,z = pos[n]  

		Ev_1 = [0,-1,-1]
		Ev_2 = [-1,0,-1]
		Ev_3 = [0,1,-1]
		Ev_4 = [1,0,-1]

		Ev_5 = [-1,-1,0]
		Ev_6 = [-1,1,0]
		Ev_7 = [1,1,0]
		Ev_8 = [1,-1,0]

		Ev_9 = [0,-1,1]
		Ev_10 = [-1,0,1]
		Ev_11 = [0,1,1]
		Ev_12 = [1,0,1]



		events[n] = [Ev_1, Ev_2, Ev_3, Ev_4, Ev_5, Ev_6, Ev_7, Ev_8, Ev_9, Ev_10, Ev_11, Ev_12]

		barrier_1 = Calc_barrier(x,y,z,x,np.mod(y-1,L),np.mod(z-1,L))		
		barrier_2 = Calc_barrier(x,y,z,np.mod(x-1,L),y,np.mod(z-1,L))
		barrier_3 = Calc_barrier(x,y,z,x,np.mod(y+1,L),np.mod(z-1,L))
		barrier_4 = Calc_barrier(x,y,z,np.mod(x+1,L),y,np.mod(z-1,L))
		barrier_5 = Calc_barrier(x,y,z,np.mod(x-1,L),np.mod(y-1,L),z)
		barrier_6 = Calc_barrier(x,y,z,np.mod(x-1,L),np.mod(y+1,L),z)
		barrier_7 = Calc_barrier(x,y,z,np.mod(x+1,L),np.mod(y+1,L),z)
		barrier_8 = Calc_barrier(x,y,z,np.mod(x+1,L),np.mod(y-1,L),z)
		barrier_9 = Calc_barrier(x,y,z,x,np.mod(y-1,L),np.mod(z+1,L))
		barrier_10 = Calc_barrier(x,y,z,np.mod(x-1,L),y,np.mod(z+1,L))
		barrier_11 = Calc_barrier(x,y,z,x,np.mod(y+1,L),np.mod(z+1,L))
		barrier_12 = Calc_barrier(x,y,z,np.mod(x+1,L),y,np.mod(z+1,L))

		barrier[n] = [barrier_1,barrier_2,barrier_3,barrier_4,barrier_5,barrier_6,barrier_7,barrier_8,barrier_9,barrier_10,barrier_11,barrier_12]
		

	return Coordination_NN,Coordination_2NN,E_Conf,events,barrier

