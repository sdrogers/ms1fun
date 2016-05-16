class HmdbMol(object):
	def __init__(self,name,mass,formula):
		self.name = name
		self.mass = mass
		self.formula = formula

class HMDB(object):
	def __init__(self,dname='../dbs/hmdb.txt'):
		self.dname = dname
		self.load_db()
		# Sort by mass (ascending)
		self.mols = sorted(self.mols,key=lambda x: x.mass)
		self.n_mols = len(self.mols)

	def load_db(self):
		self.mols = []
		count = 0
		with open(self.dname,'r') as f:
			for line in f:
				count += 1
				split_line = line.split(',')
				formula = split_line[0]
				mass = split_line[1]
				name = split_line[2]
				if len(split_line)>3:
					for i in range(len(split_line)-3):
						name += ',' + split_line[i+3]
				if mass == 'None':
					mass = 0.0
				else:
					mass = float(mass)
				self.mols.append(HmdbMol(name,mass,formula))

	def get_hits(self,M,tol=10):
		min_val = M - tol*(M/1e6)
		max_val = M + tol*(M/1e6)
		hits = []
		# Find the first value that is larger than min_val
		pos = next(i for i,mol in enumerate(self.mols) if mol.mass >= min_val)
		while self.mols[pos].mass <= max_val:
			hits.append(self.mols[pos])
			pos += 1
		return hits

