#!/anaconda/bin/python

import re
import Rate

def umist_rates(species):
	r1 = species.name
	with open("RATE12.txt", "r") as f:
		lines = f.readlines()
		f.close()
		for i,line in enumerate(lines):  
			reaction = re.split(":?", line)
			if (reaction[4]==r1 or reaction[5]==r1): 
				print "The reactants chosen are %s and %s. The reaction type is %s" % (reaction[2], reaction[3], reaction[1])
				temp = network.T
				rate = Rate.get_rate(reaction, temp)[0]
				units = Rate.get_rate(reaction, temp)[1]
				print "The reaction rate at %s K is %s %s." %(temp, rate, units)
		


				

           	