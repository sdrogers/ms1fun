class PeakGroup(object):
    
    def __init__(self,mim,mam,midm):
        self.minmass = mim
        self.maxmass = mam
        self.midmass = midm
        self.members = []
        self.M = 0.0
        self.vote = 0
        self.rt = 0.0
        
    def add_peak(self,peak,transformation):
        self.members.append((peak,transformation,transformation.transform(peak)))
        self.vote += transformation.vote
        self.M = 0.0
        toti = 0.0
        for a in self.members:
            self.M += a[0].intensity*a[2]
            self.rt += a[0].intensity*a[0].rt
            toti += a[0].intensity
        self.M /= toti
        self.rt /= toti
        
    def __str__(self):
        outstr = ""
        for p in self.members:
            outstr += str(p[2])
        return outstr
    def __repr__(self):
        outstr = ""
        for p in self.members:
            outstr += str(p[2])+","
        return outstr


class Voter(object):

    def __init__(self,transformations):
        self.transformations = transformations

    def make_groups(self,peaks):
        transformed_masses = []
        actual_transforms = []
        original_peaks = []
        for peak in peaks:
            for tr in self.transformations:
                transformed_masses.append(tr.transform(peak))
                actual_transforms.append(tr)
                original_peaks.append(peak)
        temp = zip(transformed_masses,actual_transforms,original_peaks)
        temp = sorted(temp,key = lambda x: x[0])
        transformed_masses,actual_transforms,original_peaks = zip(*temp)
        transformed_masses = list(transformed_masses)
        original_peaks = list(original_peaks)
        actual_transforms = list(actual_transforms)


        # Loop forever
        mtol = 10.0
        groups = []
        while len(transformed_masses) > 0:
            biggest = 0
            best_ind = []
            for i in range(len(transformed_masses)):
                a = transformed_masses[i]
                mid = 1e6*a/(1e6-mtol)
                up = 1e6*mid/(1e6-mtol)
                pos = i
                
                tran_ind = [i]
                temp_tran = [actual_transforms[i]]
                found_peaks = [original_peaks[i]]
                while True:
                    pos += 1
                    if pos >= len(transformed_masses):
                        break
                    if transformed_masses[pos]<=up:
                        if not original_peaks[pos] in found_peaks:
                            tran_ind.append(pos)
                            temp_tran.append(actual_transforms[pos])
                            found_peaks.append(original_peaks[pos])
                        else:
                            already = found_peaks.index(original_peaks[pos])
                            if temp_tran[already].vote < actual_transforms[pos].vote:
                                # replace the original one with this one if it has a higher vote
                                tran_ind[already] = pos
                                temp_tran[already] = actual_transforms[pos]
                    else:
                        break

                # This loop should enable second order dependencies
                finished = False
                while not finished:
                    to_remove = []
                    for k,j in enumerate(tran_ind):
                        if not actual_transforms[j].parent == None:
                            if actual_transforms[j].parent in temp_tran:
                                pass
                            else:
                                to_remove.append(k)
                    if len(to_remove) == 0:
                        finished = True
                    to_remove = sorted(to_remove,reverse=True)
                    for j in to_remove:
                        tran_ind.pop(j)

                total_vote = 0.0
                for k,j in enumerate(tran_ind):
                    total_vote += actual_transforms[j].vote
                
                if total_vote >= biggest:
                    best_ind = [j for j in tran_ind]
                    biggest = total_vote
            

            new_group = PeakGroup(a,up,mid)
            groups.append(new_group)
            best_ind = sorted(best_ind,reverse=True)
            for i in best_ind:
                new_group.add_peak(original_peaks[i],actual_transforms[i])
            peaks_involved = [original_peaks[i] for i in best_ind]
            for p in peaks_involved:
                positions = [i for i,e in enumerate(original_peaks) if e == p]
                positions = sorted(positions,reverse = True)
                for po in positions:
                    original_peaks.pop(po)
                    transformed_masses.pop(po)
                    actual_transforms.pop(po)
        return groups

class ReverseVoter(object):
    def __init__(self,transformations,peaks):
        self.transformations = transformations
        self.peaks = peaks
        self.peaks_by_rt = sorted(self.peaks,key=lambda x:x.rt)

    def find_mol(self,mol,rttol = 10,mtol=10,remove_found_peaks = False,verbose=True):
        
        if verbose:
            print "Searching for {} with mz: {} and rt: {}".format(mol.name,mol.mass,mol.rt)
        # Find peaks within the rt tolerance
        candidate_peaks = []
        predicted_masses = []
        for tr in self.transformations:
            predicted_masses.append(tr.reversetransform(mol.mass))

        for peak in self.peaks_by_rt:
            if peak.rt > mol.rt + rttol:
                break
            if peak.rt > mol.rt - rttol:
                # Check its mass
                for i,m in enumerate(predicted_masses):
                    if self.hit(m,peak.mass,mtol):
                        candidate_peaks.append((peak,self.transformations[i]))
        if verbose:
            print "Found {} explainable peaks within rt range".format(len(candidate_peaks))
            for p,t in candidate_peaks:
                print p.mass,p.rt,p.intensity,t

        if remove_found_peaks:
            for p,t in candidate_peaks:
                if p in self.peaks:
                    self.peaks.remove(p)
                if p in self.peaks_by_rt:
                    self.peaks_by_rt.remove(p)

        return candidate_peaks

    def hit(self,m1,m2,mtol):
        if 1e6*abs(m1-m2)/m1 < mtol:
            return True
        else:
            return False