import numpy as np
import numpy.random as random
import copy
import numpy.linalg as lg
import matplotlib.pyplot as plt

# INITIALISE SYSTEM HERE


L = 100
kB = 8.617E-5

Concentration = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] # Concentration of Aluminium

Al_number = 0

E_H_4Fe = 0

dE_Al_Fe = 0.3

store_D = [] #To store the diffusivities

Interstitial = np.chararray(((L-1),(L-1)),itemsize = 2)

Coordination = {}


for Al_Concentration in Concentration:



      Al_number = int(Al_Concentration*(L**2))

      LatticeF = np.chararray((L*L),itemsize = 2)

      LatticeF[:] = "Fe"   #  Fe = Iron

      LatticeF[0:Al_number] = "Al"  # Al = Aluminium

      np.random.shuffle(LatticeF)

      Lattice = np.reshape(LatticeF,(L,L))

      for i in range(len(Interstitial)):
	for j in range(len(Interstitial)):
	    
	      Coordination[i,j] = [Lattice[i][j],Lattice[i+1][j],Lattice[i][j+1],Lattice[i+1][j+1]]

      print "##########################"
      print ' Al density = ', Al_density
      print "---------------------------------------------"

      Number_of_H = 1000

      T = 500
      D = 0

      t = 0
      temp = 0 



	      
      for i in range (Number_of_H):
	      cummulative_Delta_t = 0
	      
	      steps = 0
	      nu = 1.0E13 # Attempt frequency
	      Initial_position = [random.randint(len(Interstitial)-1),random.randint(len(Interstitial)-1)] # Start always here???? random.randint (41)
	      print "Initial position =", Initial_position  
	      cell_position = copy.copy(Initial_position) 
	      real_position = copy.copy(Initial_position)
	      while (steps<1000): # Keeps on simulating till given steps
		# GET RATES:
			  
		      
		      xx,yy = cell_position   

		
		      if (cell_position[0] < 0):
			    cell_position[0] += L-2
			    
		      if (cell_position[1]<0):
			    cell_position[1] += L-2			
				      
		      if (cell_position[0] == L-2):
			    cell_position = [0,yy]
			    
		      if (cell_position[1] == L-2):
			    cell_position = [xx,0]	    
			    
		      x,y = cell_position
		      
		      N_Al_Conf0 = 0
		      N_Fe_Conf0 = 0
		      
		      Conf_0 = Coordination[x,y]
		     
		     # To calculate Different number of Neighbours 
		     
		     def NeighCount(Conf):

			N_Al = 0
			N_Fe = 0
			p = 0
			for p in Conf:
			  if(p == "Fe"):
			    N_Fe = N_Fe + 1
			  
			  elif(p == "Al"):
			    N_Al = N_Al + 1 

			return {'N_Fe':N_Fe, 'N_Al':N_Al}
		     
		     
		     N_Fe_Conf0 = NeighCount(Conf_0)['N_Fe']
		     N_Al_Conf0 = NeighCount(Conf_0)['N_Al']
			  
		      # To calculate Energy of Configurations
		     		      		      
		      def E_Conf(N_Fe,N_Al):
			
			return E_H_4Fe + N_Al * dE_Al_Fe 
		      
		      #Rate Down
		      
		      if (y==0):	
			  
			  Conf_1 = Coordination[x,L-2]
			  N_Fe = NeighCount(Conf_1)['N_Fe']
			  N_Al = NeighCount(Conf_1)['N_Al']
			  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_down =  np.exp ( - barrier/(kB * T))
		      else:

			  Conf_2 = Coordination[x,y-1]
			  N_Fe = NeighCount(Conf_2)['N_Fe']
			  N_Al = NeighCount(Conf_2)['N_Al']
			  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_down =  np.exp ( - barrier/(kB * T))
			  
		     #Rate Up
		     
		     
		      if (y==(L-2)):

			  Conf_3 = Coordination[x,0]
			  N_Fe = NeighCount(Conf_3)['N_Fe']
			  N_Al = NeighCount(Conf_3)['N_Al']
			  		  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_up =  np.exp ( - barrier/(kB * T)) 
		      
		      else:

			  Conf_4 = Coordination[x,y + 1]
			  N_Fe = NeighCount(Conf_4)['N_Fe']
			  N_Al = NeighCount(Conf_4)['N_Al']
			      

			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_up =  np.exp ( - barrier/(kB * T))
		     
		     
		     #Rate Left
		     
		      if (x==0):
			  
			  Conf_5 = Coordination[L-2,y]
			  N_Fe = NeighCount(Conf_5)['N_Fe']
			  N_Al = NeighCount(Conf_5)['N_Al']
			  			      
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_left =  np.exp ( - barrier/(kB * T))
		      else:
			  
			  Conf_6 = Coordination[x-1,y]
			  N_Fe = NeighCount(Conf_6)['N_Fe']
			  N_Al = NeighCount(Conf_6)['N_Al']
			  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  
			  rate_left =  np.exp ( - barrier/(kB * T))
			  
		      
		     #Rate Right	      
		      
		      
		      if (x==(L-2)):
			  
			  Conf_7 = Coordination[0,y]
			  N_Fe = NeighCount(Conf_7)['N_Fe']
			  N_Al = NeighCount(Conf_7)['N_Al']
			  
			  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  rate_right =  np.exp ( - barrier/(kB * T))
		      else:
			  
			  Conf_8 = Coordination[x+1,y]
			  N_Fe = NeighCount(Conf_8)['N_Fe']
			  N_Al = NeighCount(Conf_8)['N_Al']
			  
			  
			  barrier = E_Conf(N_Fe_Conf0,N_Al_Conf0) - E_Conf(N_Fe,N_Al) + E_H_4Fe + (N_Al_Conf0 + N_Al)*dE_Al_Fe*0.5
			  rate_right =  np.exp ( - barrier/(kB * T))
			  

	      
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

np.savetxt('data_bdensity_500K.txt', zip(density, store_D))

plt.plot(density,store_D, label = 'T = 500 K')
plt.xlabel('B Concentration')
plt.ylabel('Diffusivity')
plt.legend()
plt.savefig('Di_Vs_Bconc_500K.png')
plt.show()

