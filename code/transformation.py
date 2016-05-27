ELECTRON_MASS = 0.00054857990924

class Transformation(object):
	def __init__(self,name,adduct_mass,charge = 1,multiplicity = 1,fragment_mass = 0,isotope_diff = 0,parent = None,vote = 1.0,adducts = [],fragments = []):
		self.name = name
		self.adduct_mass = adduct_mass
		self.charge = charge
		self.multiplicity = multiplicity
		self.fragment_mass = fragment_mass
		self.isotope_diff = isotope_diff
		self.parent = parent
		self.vote = vote
		self.fragments = fragments
		self.adducts = adducts

	def transform(self,peak):
		M = peak.mass*abs(self.charge) + self.fragment_mass + self.charge*ELECTRON_MASS - self.adduct_mass
		# M /= (1.0*self.multiplicity)
		# M -= self.isotope_diff
		M -= self.isotope_diff
		M /= (1.0*self.multiplicity)
		return M

	def reversetransform(self,M):
		# reverse transform
		p = M
		p += self.isotope_diff
		p *= (1.0*self.multiplicity)
		p = p + self.adduct_mass - self.charge*ELECTRON_MASS - self.fragment_mass
		p = p/self.charge
		return p

	def __str__(self):
		return self.name

	def __repr__(self):
		return self.__str__()

def load_from_file(file_name):
	import yaml,re
	transformations = []
	all_vals = yaml.load(open(file_name,'r'))
	# Loop over adducts and fragments
	for tr in all_vals['transformations']:
		multiplicity = all_vals['transformations'][tr]['n']
		charge = 0
		adduct_mass = 0
		vote = all_vals['transformations'][tr]['v']
		for a in all_vals['transformations'][tr]['g']:
			if a in all_vals['charges']:
				charge += all_vals['charges'][a]
			if a in all_vals['masses']:
				adduct_mass += all_vals['masses'][a]
		adduct_list = all_vals['transformations'][tr]['g']
		transformations.append(Transformation(tr,adduct_mass,charge=charge,multiplicity=multiplicity,vote=vote,adducts = adduct_list))

		# for i in all_vals['isotopes']:
		# 	isotope_diff = all_vals['masses'][i]
		# 	newname = "{} [{}]".format(tr,i)
		# 	transformations.append(Transformation(newname,adduct_mass,charge=charge,multiplicity=multiplicity,isotope_diff=isotope_diff))

		if 'fragments' in all_vals:
			for f in all_vals['fragments']:
				fragment_mass = all_vals['masses'][f]
				if multiplicity == 1:
					newname = "[M-{}]".format(f)
				else:
					newname = "[2M-{}]".format(f)
				for a in all_vals['transformations'][tr]['g']:
					if a[0] == '-':
						newname += a
					else:
						newname += '+' + a

				frag_vote = all_vals['fragments'][f]['v']
				transformations.append(Transformation(newname,adduct_mass,charge=charge,multiplicity=multiplicity,fragment_mass=fragment_mass,vote = vote*frag_vote,adducts = adduct_list,fragments = [f]))	
				# for i in all_vals['isotopes']:
				# 	isotope_diff = all_vals['masses'][i]
				# 	newname = "{} [{}]".format(newname,i)
				# 	transformations.append(Transformation(newname,adduct_mass,charge=charge,
				# 											multiplicity=multiplicity,
				# 											isotope_diff=isotope_diff,
				# 											fragment_mass = fragment_mass))

		
	final_transformations = []
	# if 'isotopes' in all_vals:
	# 	for t in transformations:
	# 		final_transformations.append(t)
	# 		for i in all_vals['isotopes']:
	# 			new_name = t.name + '[{}]'.format(i)
	# 			isotope_diff = all_vals['masses'][i]
	# 			iso_vote = all_vals['isotopes'][i]['v']
	# 			final_transformations.append(Transformation(new_name,t.adduct_mass,charge=t.charge,
	# 														multiplicity = t.multiplicity,
	# 														isotope_diff = isotope_diff,
	# 														fragment_mass = t.fragment_mass,
	# 														parent = t,
	# 														vote = t.vote*iso_vote,
	# 														adducts = t.adducts,
	# 														fragments = t.fragments))
	for t in transformations:
		final_transformations.append(t)
		iso_vote = 0.9 # fairly random choice!
		c13_diff = 1.0033548378
		new_name = t.name + '[C13]'
		c13 = Transformation(new_name,t.adduct_mass,charge=t.charge,
			multiplicity = t.multiplicity,
			isotope_diff = c13_diff,
			fragment_mass = t.fragment_mass,
			parent = t,
			vote = t.vote*iso_vote,
			adducts = t.adducts,
			fragments=t.fragments)
		final_transformations.append(c13)
		new_name = t.name + '[2C13]'
		c13_2 = Transformation(new_name,t.adduct_mass,charge=t.charge,
			multiplicity = t.multiplicity,
			isotope_diff = 2*c13_diff,
			fragment_mass = t.fragment_mass,
			parent = c13,
			vote = c13.vote*iso_vote,
			adducts = t.adducts,
			fragments=t.fragments)

		final_transformations.append(c13_2)

		if t.multiplicity > 1:
			# add two more!
			new_name = t.name + '[3C13]'
			c13_3 = Transformation(new_name,t.adduct_mass,charge=t.charge,
			multiplicity = t.multiplicity,
			isotope_diff = 3*c13_diff,
			fragment_mass = t.fragment_mass,
			parent = c13_2,
			vote = c13.vote*iso_vote,
			adducts = t.adducts,
			fragments=t.fragments)
			final_transformations.append(c13_3)

			new_name = t.name + '[4C13]'
			c13_4 = Transformation(new_name,t.adduct_mass,charge=t.charge,
			multiplicity = t.multiplicity,
			isotope_diff = 3*c13_diff,
			fragment_mass = t.fragment_mass,
			parent = c13_3,
			vote = c13.vote*iso_vote,
			adducts = t.adducts,
			fragments=t.fragments)
			final_transformations.append(c13_4)


	# else:
	# 	final_transformations = transformations

	# print transformations
	return final_transformations


if __name__=='__main__':
	load_from_file('pos_transformations.yml')