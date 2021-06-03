import numpy as np
import numpy.random as random
import copy
import numpy.linalg as lg
import matplotlib.pyplot as plt
from utilities import *


kB = 8.617E-5 # Boltzmann Constant

Concentration = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0] # Concentration of Aluminium (in Percentage)

Al_number = 0

store_D = [] #To store the diffusivities


for Al_Concentration in Concentration:

	
      Coordination_NN,Coordination_2NN,E_Conf,events,barrier = Calc_Coordination(Al_Concentration)

      
      print "##########################"
      print ' Al Concentration(%) = ', Al_Concentration
      print "---------------------------------------------"

      D = 0
      t = 0
      temp = 0 
      
      for i in range (H):
	      cummulative_Delta_t = 0
	      nu = 1.0E13 # Attempt frequency
	      steps = 0
	       
	      Initial_position = np.array([random.randint(len(Interstitial)-1),random.randint(len(Interstitial)-1),random.randint(len(Interstitial)-1)]) 
	      print "Initial position =", Initial_position  
	      cell_position = copy.copy(Initial_position) 
	      real_position = copy.copy(Initial_position)
	      while (steps<1000): # Keeps on simulating till given steps
		# GET RATES:
			  
		    		    
			if (cell_position[0] < 0):
			    cell_position[0] = L-1
			    
			if (cell_position[1] < 0):
			    cell_position[1] = L-1

			if (cell_position[2] < 0):
			    cell_position[2] = L-1			
				      
			if (cell_position[0] > L-1):
			    cell_position[0] = 0

			if (cell_position[1] > L-1):
			    cell_position[1] = 0

			if (cell_position[2] > L-1):
			    cell_position[2] = 0
			
			x,y,z = cell_position
			n = read_pos[x,y,z]
		      
		      
		      	R_total = 0
		      	R = []
			 
		      	for j in range(len(events[n])):

 			
				Barriers = barrier[n][j]		
			  	R.append(np.exp(-Barriers/(kB*T)))
				
			
			R_total = np.sum(R)   # Sum of all rates
		
						  			      
			rho_1 = random.random_sample()
			cumm_R = 0.0

			for k in range(len(events[n])):

				cumm_R += R[k]
				
				if ((rho_1 * R_total) < cumm_R):
					cell_position += events[n][k]
					real_position += events[n][k]
					break
								
				
					
			rho_2 = random.random_sample ()

			Delta_t = -np.log (rho_2)/ (cumm_R * nu)	   
			cummulative_Delta_t += Delta_t			# Increment time
			steps += 1
	      
	      print "Final Position = ", real_position 	    
	      D_i = ((real_position[0]-Initial_position[0])**2 + (real_position[1]-Initial_position[1])**2 + (real_position[2]-Initial_position[2])**2)*(2.5e-10)**2 / (2*cummulative_Delta_t*3)     #  Diffusion Constant for each Jump	    
	      print "#",i+1,"Diffusivity(m2/s) = ", D_i
	      temp += D_i*cummulative_Delta_t
	      t += cummulative_Delta_t

      #D = temp/(2*t)
      D = temp/t
      print"Overall Diffusivity = ", D	
      store_D.append(D)

     



#Writing to a file 
np.savetxt('data_bdensity_500K.txt', zip(Concentration, store_D))

plt.plot(Concentration,store_D, label = 'T = 500 K')
plt.xlabel('B Concentration')
plt.ylabel('Diffusivity')
plt.legend()
plt.savefig('Di_Vs_Bconc_500K.png')
plt.show()
