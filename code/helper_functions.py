
import pylab as plt
import numpy as np

def collect_group_stats(groups,thresh = 1e5):
	n_peaks = 0
	n_singletons = 0
	n_singletons_under_thresh = 0
	non_singletons_under_thresh = 0
	tiny_groups = 0
	tiny_size = 0
	tiny_vote = 0
	votes = []
	intensities = []
	maxivotes = []
	maxi = []
	for group in groups:
		n_peaks += len(group.members)
		if len(group.members) == 1:
			n_singletons += 1
			if group.members[0][0].intensity < thresh:
				n_singletons_under_thresh += 1
			intensities.append(group.members[0][0].intensity)
			votes.append(group.vote)
			maxi.append(group.members[0][0].intensity)
			maxivotes.append(group.vote)
		else:
			mi = 0
			for p,_,_ in group.members:
				if p.intensity < thresh:
					non_singletons_under_thresh += 1
				intensities.append(p.intensity)
				votes.append(group.vote)
				if p.intensity > mi:
					mi = p.intensity
			if mi < thresh:
				tiny_groups += 1
				tiny_size += len(group.members)
				tiny_vote += group.vote
			maxi.append(mi)
			maxivotes.append(group.vote)
	print "{} groups, consisting of {} peaks".format(len(groups),n_peaks)
	print "{} singleton groups ({:.0f}% of groups, {:.0f}% of peaks)".format(n_singletons,100.0*n_singletons/len(groups),100.0*n_singletons/n_peaks)
	print "{} peaks under threshold ({:.0e}), {:.0f}% of peaks".format(n_singletons_under_thresh + non_singletons_under_thresh,thresh,100.0*(n_singletons_under_thresh+non_singletons_under_thresh)/n_peaks)
	print "\t{} of which are singletons ({:.0f}%)".format(n_singletons_under_thresh,100.0*n_singletons_under_thresh/(n_singletons_under_thresh+non_singletons_under_thresh))
	print "{} peaks below the threshold in groups of size > 1".format(non_singletons_under_thresh)
	print "{} groups where the most intense peak is below the threshold (avg size = {:.2f} avg vote = {:.2f})".format(tiny_groups,1.0*tiny_size/tiny_groups,1.0*tiny_vote/tiny_groups)
	plt.figure()
	plt.plot(np.log(intensities),votes,'k.')
	from scipy.stats.stats import pearsonr
	r,p = pearsonr(np.log(intensities),votes)
	print "Test between intensity and vote for all peaks: corr coef = {}, p-value = {}".format(r,p)
	plt.xlabel('Log intensity')
	plt.ylabel('Vote of enclosing group')


	plt.figure()
	plt.plot(np.log(maxi),maxivotes,'r.')
	from scipy.stats.stats import pearsonr
	r,p = pearsonr(np.log(maxi),maxivotes)
	print "Test between maximum group intensity and vote: corr coef = {}, p-value = {}".format(r,p)
	plt.xlabel('Log maximum group intensity')
	plt.ylabel('Vote of group')

