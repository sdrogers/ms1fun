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
                while True:
                    pos += 1
                    if pos >= len(transformed_masses):
                        break
                    if transformed_masses[pos]<=up:
                        tran_ind.append(pos)
                        temp_tran.append(actual_transforms[pos])
                    else:
                        break

                to_remove = []
                for k,j in enumerate(tran_ind):
                    if not actual_transforms[j].parent == None:
                        if actual_transforms[j].parent in temp_tran:
                            pass
                        else:
                            to_remove.append(k)
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


