#!/anaconda/bin/python

import re
import Rate

def umist_rates(species, network):
	r1 = species.name
	with open("RATE12.txt", "r") as f:
		lines = f.readlines()
		f.close()
		for i,line in enumerate(lines):  
			reaction = re.split(":?", line)
			if (reaction[2]==r1 or reaction[3]==r1 or reaction[4]==r1 or reaction[5]==r1): 
				temp = network.T
				Rate.get_rate(reaction, temp)
		


				

           	
