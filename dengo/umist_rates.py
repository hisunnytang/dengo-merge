#!/anaconda/bin/python

import re
import Rate

def umist_rates(speciesOne, speciesTwo):
	r1 = speciesOne.name
	r2 = speciesTwo.name
	with open("RATE12.txt", "r") as f:
		lines = f.readlines()
		f.close()
		for i,line in enumerate(lines):  
			reaction = re.split(":?", line)
			if (reaction[2]==r2 or reaction[3]==r2) and (reaction[2]==r1 or reaction[2]==r1): 
				print "The reactants chosen are %s and %s. The reaction type is %s" % (reaction[2], reaction[3], reaction[1])
				temp = network.T
				rate = Rate.get_rate(reaction, temp)[0]
				units = Rate.get_rate(reaction, temp)[1]
				print "The reaction rate at %s K is %s %s." %(temp, rate, units)
		


				

           	