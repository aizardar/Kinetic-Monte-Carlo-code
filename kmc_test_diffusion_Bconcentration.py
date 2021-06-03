import numpy as np
import numpy.random as random
import copy
import numpy.linalg as lg
import matplotlib.pyplot as plt

# Very simple kmc code

#  Fe = Iron , Al = Aluminium


barriers [('i', 'i')] = 0.25
barriers [('i', 'j')] = 0.275
barriers [('j', 'i')] = 0.275
barriers [('j', 'j')] = 0.3



# INITIALISE SYSTEM HERE
L = 10
kB = 8.617E-5
#density = np.round(np.linspace(0,1,10),1)
density = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
Al_number = 0;
store_D = [] #To store the diffusivities
for Al_density in density:

    Al_number = int(Al_density*(L**2))
    LatticeF = np.chararray((L*L),itemsize = 2)
    LatticeF[:] = "Fe"
    LatticeF[0:Al_number] = "Al"
    np.random.shuffle(LatticeF)
    Lattice = np.reshape(LatticeF,(L,L))

    print "##########################"
    print ' Al density = ', Al_density
    print "---------------------------------------------"
    Number_of_H = 1000

    T = 1000
    D = 0
    
    t = 0
    temp = 0 
    
    for i in range (Number_of_H):
	    cummulative_Delta_t = 0
	    
	    steps = 0
	    nu = 1.0E13 # Attempt frequency
	    Initial_position = [random.randint(L),random.randint(L) ] # Start always here???? random.randint (41)
	    print "Initial position =", Initial_position  
	    cell_position = copy.copy(Initial_position)  
	    real_position = copy.copy(Initial_position)
	    while (steps<1000): # Keeps on simulating till given steps
	      # GET RATES:
			
		    
		    xx,yy = cell_position   

	      
		    if (cell_position[0]<0):
			  cell_position[0] += L
			  
		    if (cell_position[1]<0):
			  cell_position[1] += L			
		    		     
		    if (cell_position[0] == L):
			  cell_position = [0,yy]
			  
		    if (cell_position[1] == L):
			  cell_position = [xx,0]	    
			  
		    x,y = cell_position
					
		    if (y==0):
		      
			
			rate_down =  np.exp ( - barriers [(symbols[x][y], symbols[x][L-1])]/(kB * T))
		    else:
			rate_down =  np.exp ( - barriers [(symbols[x][y], symbols[x][y-1])]/(kB * T)) 
		    if (y==(L-1)):
			rate_up =  np.exp ( -barriers [(symbols[x][y], symbols[x][0])]/(kB * T))
		    else:
			rate_up =  np.exp ( -barriers [(symbols[x][y], symbols [x][y+1])]/(kB * T))
		    if (x==0):
			rate_left = np.exp ( -barriers [(symbols[x][y], symbols[L-1][y])]/(kB * T))
		    else:
			rate_left =  np.exp ( - barriers [(symbols[x][y], symbols[x-1][y])]/(kB * T))
		    if (x==(L-1)):
			rate_right =  np.exp ( -barriers [(symbols[x][y], symbols[0][y])]/(kB * T))
		    else:
			rate_right =  np.exp ( -barriers [(symbols[x][y], symbols[x+1][y])]/(kB * T))

	    
		    R = rate_up + rate_down + rate_left + rate_right # Total rate

		    rho_1 = random.random_sample()

		    if ((rho_1 * R) < rate_up): # move up
			cell_position[1] +=1
			real_position[1] +=1
		    elif ((rho_1 * R) < (rate_up + rate_down)): # move down
			cell_position[1] -= 1
			real_position[1] -= 1
		    elif ((rho_1 * R) < (rate_up + rate_down + rate_left)): # move left
			cell_position[0] -= 1
			real_position[0] -= 1
		    else: # move right
			cell_position[0] += 1
			real_position[0] += 1

	      
		    rho_2 = random.random_sample ()

		    Delta_t = -np.log (rho_2)/ (R * nu)	   
		    cummulative_Delta_t += Delta_t			# Increment time
		    steps += 1
	    
	    print "Final Position = ", real_position 	    
	    D_i = ((real_position[0]-Initial_position[0])**2 + (real_position[1]-Initial_position[1])**2)*(2.5e-10)**2 / (2*cummulative_Delta_t)     #  Diffusion Constant for each Jump	    
	    print "#",i+1,"Diffusivity(m2/s) = ", D_i
	    temp += D_i*cummulative_Delta_t
	    t += cummulative_Delta_t

    D = temp/(2*t) 
    print"Overall Diffusivity = ", D	
    store_D.append(D)

#Writing to a file 

np.savetxt('data_bdensity_1000K.txt', zip(density, store_D))

plt.plot(density,store_D, label = 'T = 1000 K')
plt.xlabel('B Concentration')
plt.ylabel('Diffusivity')
plt.legend()
plt.savefig('Di_Vs_Bconc_1000K.png')
plt.show()
