
class HmdbMol(object):
	def __init__(self,name,mass,formula):
		self.name = name
		self.mass = mass
		self.formula = formula

formulas = []
count = 0
with open('hmdb_test.txt','r') as f:
	for line in f:
		count += 1
		split_line = line.split(',')
		formula = split_line[0]
		mass = split_line[1]
		print mass,count,line
		name = split_line[2]
		if len(split_line)>3:
			for i in range(len(split_line)-3):
				name += ',' + split_line[i+3]
		if mass == 'None':
			mass = 0.0
		else:
			mass = float(mass)
		formulas.append(HmdbMol(name,mass,formula))

formulas = sorted(formulas,key = lambda x: x.mass)

with open('hmdb_sorted.txt','w') as f:
	for form in formulas:
		line = "{},{},{}".format(form.formula,form.mass,form.name)
		f.write(line)

