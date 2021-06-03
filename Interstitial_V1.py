import numpy as np
import numpy.random as random
import copy


# INITIALISE SYSTEM HERE
L = 10

Al_number = int(0.5*(L**2))
Lattice = np.chararray((L*L),itemsize = 2)
Interstitial = np.chararray((L*L)/4)
Lattice[:] = "Fe"
Lattice[:Al_number] = "Al"
#np.random.shuffle(Lattice)
#Lattice = np.reshape(LatticeF,(L,L))
j =  0 
i = 0
Neighbours_type_1 = ["Al","Al","Al","Al"]
Neighbours_type_2 = ["Fe","Fe","Fe","Fe"]
#Inters_type_3 = ["Fe","Al","Al","Al"]

while(i<len(Interstitial)):

    Neighbours = Lattice[j:j+4]
    #print Neighbours
    a = 0 
    s1 = 0
    s2 = 0
    for a in range(0,len(Neighbours)-1):
      
      if Neighbours[a] == Neighbours_type_1[a]:
	s1 = s1 + 1
	#print s1
      if Neighbours[a] == Neighbours_type_2[a]:
	s2 = s2 + 1
	#print s2
      
      a = a + 1
                
    if s1 == 3:
	print "Interstitial Site type 1 found"
    elif s2 == 3:
	print "Interstitial Site type 2 found"
    else:
	print "Interstitial Site Not defined"

    
    j = j + 4
    i = i + 1


