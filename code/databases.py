from formula import Formula

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

	def get_group_hits(self,groups,tol = 10):
		hits = {}
		for g in group:
			temp_hits = self.get_hits(g.M,tol = tol)
			if len(temp_hits)>0:
				hits[g] = temp_hits
		return hits

# Standards

class Mol(object):
	def __init__(self,name,formula,mass,rt):
		self.name = name.strip()
		self.formula = formula.strip()
		self.mass = mass
		self.rt = rt
	def __str__(self):
	    return "{} ({},{})".format(self.name,self.mass,self.rt)
	def __repr__(self):
	    return "{} ({},{})".format(self.name,self.mass,self.rt)

class Standards(object):
	def __init__(self,std_file = '/Users/simon/Dropbox/Bioresearch/Meta_clustering/ms1fundata/Beer/Std1_1_20150422_150810_combined.csv'):
		self.std_file = std_file
		self.mols = []
		self.load()

	def load(self):
		with open(self.std_file,'rU') as f:
		    for i in range(9):
		        f.readline() # remove heads
		    for line in f:
		        split_line = line.split(',')
		        polarity = split_line[4]
		        rts = split_line[9] # observed, not predicted
		        if rts == '-':
		            rt = 0.0
		        else:
		            rt = float(rts)
		        if polarity == '+' and rt > 0.0:
		            name = split_line[2]
		            formula = split_line[3]
		            f = Formula(formula)
		            new_mol = Mol(name,formula,f.compute_exact_mass(),rt*60.0)
		            self.mols.append(new_mol)
			self.mols = sorted(self.mols,key = lambda x: x.mass)

	def get_group_hits(self,groups,use_max_vote=True,mtol = 10,rttol = 10):
	    hits = {}
	    for g in groups:
	        for mol in self.mols:
	            if self.hit(mol.mass,g.M,mol.rt,g.rt,mtol=mtol,rttol=rttol):
	                if not mol in hits:
	                    hits[mol] = g
	                else:
	                    if use_max_vote:
	                        if g.vote > hits[mol]:
	                            hits[mol] = g
	                    else:
	                        current = hits[mol]
	                        if abs(g.rt - mol.rt) < abs(current.rt - mol.rt):
	                            hits[mol] = g
	                
	    return hits

	def hit(self,m1,m2,rt1,rt2,mtol,rttol):
	    if 1e6*abs(m1-m2)/m2 < mtol and abs(rt1-rt2)<rttol:
	        return True
	    else:
	        return False

	def get_hits(self,mass,rt,mtol=10,rttol=10):
		hits = []
		for mol in self.mols:
			if self.hit(mass,mol.mass,rt,mol.rt,mtol = mtol,rttol = rttol):
				hits.append(mol)
		return hits


