import numpy as np
import numpy.random as random

# system parameters
L =10     		# number of unit cells
number_of_H = 1000	# Number of hydrogen atoms
max_steps = 10000	# maximum kmc steps per atom.
T = 300     		# Temperature(T)
E_H_Fe = 0.0     	# Reference Energy (All Fe)
dE_H_FeFe = 0.55    	# Activation Energy for jump from all Fe to all Fe Configuration
dE_H_Al = 0.2     	# Energy difference between 5Fe1Al and 6Fe configuration
dE_H_Al_2 = -0.03   	# Energy of the interstitial site when 1 Al atom is 2nd nearest neighbour


#Interstitial = np.chararray((L,L,L),itemsize = 2)
Coordination_NN = {}        # To store the Coordination of first nearest neighbours 
Coordination_2NN = {}       # To store the Coordination of second nearest neighbours
E_Conf = {}		# To store the Configuration energies  
events = {}         	# To store the events
pos = {}
read_pos = {}
LatticeF = []



# To calculate number of Al atoms 

def NeighCount(Coordination):

    N_Al = 0

    for p in Coordination:
      if(p == "Al"):

        N_Al = N_Al + 1
      
    return N_Al


# To calculate the barrier

#def Calc_barrier(x,y,z,xp,yp,zp):
def Calc_barrier(n, nprime):
        Delta_E = E_Conf[n] - E_Conf[nprime]
        if (Delta_E<0):
            return -Delta_E + dE_H_FeFe  + NeighCount(Coordination_NN[n])*dE_H_Al + NeighCount(Coordination_NN[nprime])*dE_H_Al
        else:
            return  dE_H_FeFe+ NeighCount(Coordination_NN[n])*dE_H_Al + NeighCount(Coordination_NN[nprime])*dE_H_Al
    

# To calculate the Coordinations, Energies, events, and barriers


def Calc_Coordination(Al_Concentration):
    TotSites = (4 * L**3)
    LatticeF = np.chararray((TotSites),itemsize = 2)  # 4 atoms per unit cell in the FCC lattice
    Al_number = int((Al_Concentration*(TotSites))/100.0) 
    LatticeF[:] = "Fe"   #  Fe = Iron
    LatticeF[0:Al_number] = "Al"  # Al = Aluminium
    np.random.shuffle(LatticeF) # random redistribution of Al
    
    print 'Al_concentration', Al_Concentration
    print 'System size', 4  * L**3
    print 'Number of Al atoms', Al_number

    Lat_site = {}; Lat_Index = {}; Inter_site = {}; Inter_Index = {}

    barrier = {}; event_index = {}      # To store the barriers
#   Lattice = np.reshape(LatticeF,(L,L,L))
    
    
    n = 0 

    for x in range (L):
      for y in range (L):
        for z in range (L):
          # 4 lattice sites first
          Lat_site[n] = [2 * x, 2 * y, 2 * z]
          Lat_site[n+1] = [2 * x+ 1, 2 * y + 1, 2 * z]
          Lat_site[n+2] = [2 * x+1, 2 * y , 2 * z + 1]
          Lat_site[n+3] = [2 * x, 2 * y +1 , 2 * z + 1]
          
          Lat_Index[2 * x, 2 * y, 2 * z] = n
          Lat_Index[2 * x + 1, 2 * y+1, 2 * z] = n+1
          Lat_Index[2 * x +1, 2 * y, 2 * z+1] = n+2
          Lat_Index[2 * x, 2 * y+1, 2 * z+1] = n+3
          
          # Octahedral sites next
          
          Inter_site[n] = [2 * x+1, 2 * y, 2 * z]
          Inter_site[n+1] = [2 * x, 2 * y + 1, 2 * z]      
          Inter_site[n+2] = [2 * x, 2 * y , 2 * z + 1]
          Inter_site[n+3] = [2 * x+1, 2 * y +1 , 2 * z + 1]
          
          Inter_Index[2 * x+1, 2 * y, 2 * z] = n
          Inter_Index[2 * x, 2 * y+1, 2 * z] = n+1
          Inter_Index[2 * x, 2 * y, 2 * z+1] = n+2
          Inter_Index[2 * x+1, 2 * y+1, 2 * z+1] = n+3

          n += 4

        # Setting up configuration
    for n in range(TotSites):
        x,y,z = Inter_site[n]

#        Coordination_NN[n] = [Lat_Index[np.mod(x+1, 2 * L), y, z], Lat_Index[np.mod (x-1, 2* L), y, z], Lat_Index [x, np.mod (y+1, 2 * L), z], Lat_Index [x, np.mod (y-1, 2 * L), z], Lat_Index[x, y, np.mod (z+1, 2 * L)], Lat_Index[x,y, np.mod (z-1, 2 * L)]]
        Coordination_NN[n] = [LatticeF[Lat_Index[np.mod(x+1, 2 * L), y, z]], LatticeF[Lat_Index[np.mod (x-1, 2* L), y, z]], LatticeF[Lat_Index [x, np.mod (y+1, 2 * L), z]], LatticeF[Lat_Index [x, np.mod (y-1, 2 * L), z]], LatticeF[Lat_Index[x, y, np.mod (z+1, 2 * L)]], LatticeF[Lat_Index[x,y, np.mod (z-1, 2 * L)]]]
        E_Conf[n] = E_H_Fe + NeighCount(Coordination_NN[n]) * dE_H_Al
#        print 'Energy', n, E_Conf[n]

    for n in range (TotSites):  
        x, y, z = Inter_site[n]
        events[n] = [[1,1,0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, 1, -1], [0,-1, 1], [0,-1, -1]]

        barrier_list = []
        event_list = []
        for event in (events[n]):
             x_prime, y_prime, z_prime =  np.mod (x + event[0], 2 * L), np.mod (y + event[1], 2 * L), np.mod (z + event[2], 2 * L)
             barrier_list += [Calc_barrier (n, Inter_Index[ x_prime, y_prime, z_prime]) ] 
             event_list += [ Inter_Index [x_prime, y_prime, z_prime]]
           
        barrier[n] = barrier_list
        event_index[n] = event_list
    return TotSites, Inter_site, events, event_index,  barrier
    
