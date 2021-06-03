import numpy as np
import numpy.random as random
import copy
import numpy.linalg as lg
import matplotlib.pyplot as plt
import time
# kmc code in 2D (Cannot be generalised !!!!!)

startTime = time.time()
Atoms = {}
barriers = {}
Atoms [('Al', 'Al')] = 'A'
Atoms [('Fe', 'Fe')] = 'A'
Atoms [('Fe', 'Al')] = 'B'
Atoms [('Al', 'Fe')] = 'B'


barriers [('A', 'A')] = 0.25
barriers [('A', 'B')] = 0.275
barriers [('B', 'A')] = 0.275
barriers [('B', 'B')] = 0.3

# INITIALISE SYSTEM HERE
L = 100
kB = 8.617E-5
#density = np.round(np.linspace(0,1,10),1)
density = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
Al_number = 0;
store_D = [] #To store the diffusivities

for Al_density in density:

    Al_number = int(Al_density*(L**2))
    OneD = np.chararray((L*L))
    OneD[:] = 'Fe'
    OneD[0:Al_number] = 'Al'
    np.random.shuffle(OneD)
    #twoD = np.reshape(OneD,(L,L))
    #OneD = np.ndarray.flatten(symbols)
    print "##########################"
    print ' Al density = ', Al_density
    print "---------------------------------------------"
    Number_of_H = 1000

    T = 300
    D = 0
    t = 0
    temp = 0 
    
    for i in range (Number_of_H):
	    cummulative_Delta_t = 0
	    steps = 0
	    nu = 1.0E13 # Attempt frequency
	    Initial_position = [random.randint(L),random.randint(L)] # Start always here???? random.randint (41)
	    print "Initial position =", Initial_position  
	    cell_position = copy.copy(Initial_position)  
	    real_position = copy.copy(Initial_position)
	    while (steps<1000): # Keeps on simulating till given steps
	      # GET RATES:
			
		    
		    x,y = cell_position   
                    
                    
		    if (x<0):
			  x += L
			  
		    if (y<0):
			 y = y + L     
			  
		    if (x == L):
			  x = 0
			  
		    if (y == L):
			  y = 0	    
			  
		    #x,y = cell_position
		   
		   
		    d = x + y*L
						
		  					
		    if (y==0):
			rate_up =  np.exp ( - barriers [(OneD[d],OneD[d + L * y])]/(kB * T))     
		    else:
			rate_up =  np.exp ( - barriers [(OneD[d], OneD[d - L])]/(kB * T)) 
		    if (y==(L-1)):
			rate_down =  np.exp ( -barriers [(OneD[d], OneD[d - L * y])]/(kB * T))
		    else:
			rate_down =  np.exp ( -barriers [(OneD[d], OneD[d + L])]/(kB * T))
		    if (x==0):
			rate_left = np.exp ( -barriers [(OneD[d], OneD[d + (L - 1)])]/(kB * T))
		    else:
			rate_left =  np.exp ( - barriers [(OneD[d], OneD[d - 1])]/(kB * T))
		    if (x==(L-1)):
			rate_right =  np.exp ( -barriers [(OneD[d], OneD[d - (L-1) ])]/(kB * T))
		    else:
			rate_right =  np.exp ( -barriers [(OneD[d], OneD[d + 1])]/(kB * T))

	    
		    R = rate_up + rate_down + rate_left + rate_right # Total rate

		    rho_1 = random.random_sample()

		    if ((rho_1 * R) < rate_up): # move up
			#d = d - L
			y = y - 1 
			real_position[1] += 1
		    elif ((rho_1 * R) < (rate_up + rate_down)): # move down
			#d = d + L
			y = y + 1 
			real_position[1] -= 1
		    elif ((rho_1 * R) < (rate_up + rate_down + rate_left)): # move left
			#d = d -1
			x = x - 1
			real_position[0] -= 1
		    else: # move right
			#d = d + 1
			x = x + 1 
			real_position[0] += 1

	      
		    rho_2 = random.random_sample ()

		    Delta_t = -np.log (rho_2)/ (R * nu)	   
		    cummulative_Delta_t += Delta_t			# Increment time
		    steps += 1
		    cell_position = x,y  
		    
	    
	    print "Final Position = ", real_position 	    
	    D_i = ((real_position[0]-Initial_position[0])**2 + (real_position[1]-Initial_position[1])**2)*(2.5e-10)**2 / (2*cummulative_Delta_t)     #  Diffusion Constant for each Jump	    
	    print "#",i+1,"Diffusivity(m2/s) = ", D_i
	    temp += D_i*cummulative_Delta_t
	    t += cummulative_Delta_t

    D = temp/(2*t)
    print"Overall Diffusivity = ", D	
    store_D.append(D)


elapsedTime = time.time() - startTime
print "Elapsed Time = ",elapsedTime
#Writing to a file 

np.savetxt('data_bdensity_300K.txt', zip(density, store_D))

plt.plot(density,store_D, label = 'T = 300K')
plt.xlabel('B Concentration')
plt.ylabel('Diffusivity')
plt.legend()
plt.savefig('Di_Vs_Bconc_300K.png')
plt.show()
