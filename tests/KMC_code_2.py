import numpy as np
import numpy.random as random
import copy
import numpy.linalg as lg
from utilities import *


kB = 8.617E-5 # Boltzmann Constant

#Concentration = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0] # Concentration of Aluminium (in Percentage)
Concentration = np.linspace (0.0,  10.0, 41)
print 'List of concentrations', Concentration
store_D = [] #To store the diffusivities

f = open ('D_vs_Al_concentration_300K.dat', 'w')
print >> f, '# Al conc (at.%) -- D (m^2/s)'
f.close()
for Al_Concentration in Concentration:

    
      N_sites, positions,  events, event_list, barrier = Calc_Coordination(Al_Concentration) # combine events, event_list and barrier later
      print "##########################"
      print ' Al Concentration(%) = ', Al_Concentration
      print "---------------------------------------------"

      D = 0.0
      t = 0.0
      temp = 0.0
      temp_average = 0.0
      D_median = 0.0
     
      for i in range (number_of_H):
          cumulative_Delta_t = 0
          nu = 1.0E13 # Attempt frequency # should be elsewhere in the code!!
          steps = 0
           
          site = random.randint (N_sites)
          r_0 = copy.copy(positions [site])
          print "Initial position =", r_0
#         r_cell = r_0
          r_real = copy.copy (r_0)
          #print "Initial position = ", r_real
          while (steps<max_steps): # Keeps on simulating till given steps - should be a parameter of the model!!!
            R = []           
            for j in range(len(events[site])):
               Barriers = barrier[site][j]     
               R += [np.exp(-Barriers/(kB*T))]               
            
            R_total = sum(R)   # Sum of all rates
        
                                          
            rho_1 = random.random_sample()
            cum_R = 0.0

            for k in range(len(events[site])):

                cum_R += R[k]
                
                if ((rho_1 * R_total) < cum_R):
##                    print 'Event', steps, events[site][k]
##                    norm_delta = np.max (delta)
##                    if (norm_delta>2):
##                          print 'Error?'
##                    print steps, site, positions[site], k, events[site][k], r_real
                    for coord in range (3):                         
                        r_real[coord] += events[site][k][coord]
                    site = event_list[site][k]

                    break                               
            rho_2 = random.random_sample ()

            Delta_t = -np.log (rho_2)/ (R_total * nu)      
            cumulative_Delta_t += Delta_t           # Increment time
            steps += 1
          
          print "Final Position = ", r_real     
         # print "Distance travelled in real space",  np.sqrt ((r_real[0]-r_0[0])**2 + (r_real[1]-r_0[1])**2 + (r_real[2]-r_0[2])**2)
          print "Time taken", cumulative_Delta_t
          D_i = ((r_real[0]-r_0[0])**2 + (r_real[1]-r_0[1])**2 + (r_real[2]-r_0[2])**2)*(2.5e-10)**2 / (2*cumulative_Delta_t*3)     #  Diffusion Constant for each Jump    
          print "#",i+1,"Diffusivity(m2/s) = ", D_i
          
          temp += D_i *cumulative_Delta_t
          temp_average += D_i
          t += cumulative_Delta_t
          

      D = temp/(2*t)
      print"Overall Diffusivity = ", D
      
      f = open ('D_vs_Al_concentration_300K.dat', 'a')
      print >> f, Al_Concentration, D, temp_average/number_of_H
      print Al_Concentration, D, temp_average/number_of_H, 
      f.close()


     

