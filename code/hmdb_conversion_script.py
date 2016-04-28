# Loads the individual xml files from hmdb and then extracts the name, mass and chemical formula

import os
import xml.etree.ElementTree
# Get the list of files



prefix = '/Users/simon/Downloads/hmdb_metabolites/'
files = os.listdir(prefix)
files = [f for f in files if f.startswith('HMDB')]

outfile = 'hmdb_test.txt'
with open(outfile,'w') as f:


	for i,file in enumerate(files):
		if i%100 == 0:
			print "{} out of {}".format(i,len(files))

		e = xml.etree.ElementTree.parse(prefix+file).getroot()

		name = e.findall('name')[0].text
		formula = e.findall('chemical_formula')[0].text
		mass = e.findall('monisotopic_moleculate_weight')[0].text

		out_line = "{},{},{}\n".format(formula,mass,name.encode('ascii','ignore'))
		f.write(out_line)


