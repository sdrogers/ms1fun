
import pylab as plt
import numpy as np
import copy
from random import shuffle

PROTON = 1.00727645199076


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


def standard_analysis(groups,mtol=10,rttol=10,mode='pos',filename='blah',v_thresh = np.arange(1.0,7.0,0.5)):
	outpre = 'output/' + mode + '/' + filename
	from databases import Standards
	st = Standards()
	group_hits = 0
	raw_hits = 0
	if mode == 'pos':
		mu = -PROTON
	else:
		mu = PROTON

	vh = []
	ns_group_peak_hits = []

	ns_groups = 0
	ns_group_hits = 0

	n_peaks = 0

	for group in groups:
		if len(group.members) > 1:
			ns_groups += 1

		for p,_,_ in group.members:
			n_peaks += 1
			hits = st.get_hits(p.mass + mu,p.rt,mtol=mtol,rttol=rttol)
			if len(hits) > 0:
				raw_hits += 1
			if len(group.members) > 1:
				if len(hits) > 0:
					ns_group_peak_hits.append(1)
				else:
					ns_group_peak_hits.append(0)

		hits = st.get_hits(group.M,group.rt,mtol=mtol,rttol=rttol)
		if len(hits) > 0:
			h = 1
			group_hits += 1
			if len(group.members) > 1:
				ns_group_hits += 1
		else:
			h = 0
		vh.append((group.vote,h))

	print "Number of raw hits (i.e. comparing all peaks): {} ({:.0f}% of peaks)".format(raw_hits,100.0*raw_hits/n_peaks)
	print "Number of group hits (i.e. hits on group Ms): {} ({:.0f}% of groups)".format(group_hits,100.0*group_hits/len(groups))

	make_thresh_plots(vh,v_thresh,outpre)
	perm_plot(ns_group_peak_hits,ns_groups,ns_group_hits,outpre)


def hmdb_analysis(groups,mtol=10,mode = 'pos',filename = 'blah',v_thresh = np.arange(1.0,7.0,0.5),hmdb_filter = []):

	outpre = 'output/' + mode + '/' + filename

	from databases import HMDB
	hmdb = HMDB(hmdb_filter = hmdb_filter)
	# Find the hits across all peaks at mtol
	group_hits = 0
	raw_hits = 0
	if mode == 'pos':
		mu = -PROTON
	else:
		mu = PROTON

	vh = []
	ns_group_peak_hits = []

	ns_groups = 0
	ns_group_hits = 0

	n_peaks = 0

	from transformation import Counts
	counts = Counts()

	for group in groups:
		if len(group.members) > 1:
			ns_groups += 1

		for p,_,_ in group.members:
			n_peaks += 1
			hits = hmdb.get_hits(p.mass + mu,tol=mtol)
			if len(hits) > 0:
				raw_hits += 1
			if len(group.members) > 1:
				if len(hits) > 0:
					ns_group_peak_hits.append(1)
				else:
					ns_group_peak_hits.append(0)

		hits = hmdb.get_hits(group.M,tol=mtol)
		if len(hits) > 0:
			h = 1
			group_hits += 1
			if len(group.members) > 1:
				ns_group_hits += 1
			for p,t,_ in group.members:
				counts.update(t)
		else:
			h = 0
		vh.append((group.vote,h))


	print "Number of raw hits (i.e. comparing all peaks): {} ({:.0f}% of peaks)".format(raw_hits,100.0*raw_hits/n_peaks)
	print "Number of group hits (i.e. hits on group Ms): {} ({:.0f}% of groups)".format(group_hits,100.0*group_hits/len(groups))


	make_thresh_plots(vh,v_thresh,outpre)
	# Make permuatation plot
	# i.e. of all the peaks in non-singleton groups, if we randomly selected
	# M from them, what would we get
	perm_plot(ns_group_peak_hits,ns_groups,ns_group_hits,outpre)
	return counts


def make_thresh_plots(vh,v_thresh,outpre):
	# make the threshold v hits plot
	overs = []
	unders = []
	sizes = []
	for vth in v_thresh:
	    over = []
	    under = []
	    for v,h in vh:
	        if v > vth:
	            over.append(h)
	        else:
	            under.append(h)

	    over = np.array(over)
	    if len(over) == 0:
	    	overs.append(0)
	    else:
	    	overs.append(over.mean())
	    under = np.array(under)
	    unders.append(under.mean())
	    sizes.append(len(over))

	overs = np.array(overs)
	unders = np.array(unders)
	
	plt.figure(figsize=(10,10))
	plo = []
	overplot, = plt.plot(v_thresh,overs,label='Prop over thresh with hits')
	underplot, = plt.plot(v_thresh,unders,label='Prop under thresh with hits')
	plt.xlabel('Vote threshold')
	plt.ylabel('Proportion of groups')
	plt.grid()
	plt.legend(handles=[overplot,underplot])
	plt.title('Proportion of groups with hits above and below (or equal to) vote threshold')
	plt.savefig(outpre + '_thresh_hits.png',dpi=100)
	plt.figure()
	plt.plot(v_thresh,sizes)
	plt.xlabel('Vote threshold')
	plt.ylabel('Number of groups')
	plt.title('Proportion of groups above vote threshold')
	plt.savefig(outpre + '_thresh_groups.png',dpi=100)


def perm_plot(ns_group_peak_hits,ns_groups,ns_group_hits,outpre):
	counts = []
	ns_copy = copy.deepcopy(ns_group_peak_hits)
	for i in range(10000):
		shuffle(ns_copy)
		counts.append(sum(ns_copy[:ns_groups]))

	plt.figure()
	plt.hist(counts)
	_,ymax = plt.ylim()

	plt.plot([ns_group_hits,ns_group_hits],[0,ymax],'r',linewidth=2)
	plt.ylabel('count')
	plt.xlabel('number of hits')
	plt.title('True number of hits in non-singleton groups, versus randomly chosen')
	plt.savefig(outpre + '_permute.png',dpi=100)

