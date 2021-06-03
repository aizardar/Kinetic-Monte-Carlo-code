import random 
import multiprocessing
from multiprocessing import Pool
import numpy as np

 
n = 100000000

def monte_carlo_pi_part(n):

	count = 0
	for i in range(n):
		x = random.random()
		y = random.random()

		if x*x + y*y <=1:
			count = count + 1 

	return count

summ = monte_carlo_pi_part(n)

print "Estimated value of Pi:: ", np.sum(summ)/(n*1.0)*4
