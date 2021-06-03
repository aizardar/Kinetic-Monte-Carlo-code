import random 
import multiprocessing
from multiprocessing import Pool
import numpy as np

 

def monte_carlo_pi_part(n):

	count = 0
	for i in range(n):
		x = random.random()
		y = random.random()

		if x*x + y*y <=1:
			count = count + 1 

	return count

if __name__ == "__main__":

	np = multiprocessing.cpu_count()
	print "You have {0:1d} CPUs".format(np)

	n = 100000000

	part_count = [n/np for i in range(np)]
	print part_count

	pool = Pool(processes = np)


	count = pool.map(monte_carlo_pi_part, part_count)
	print count
	print " Estimated value of Pi:: ", sum(count)/(n*1.0)*4
