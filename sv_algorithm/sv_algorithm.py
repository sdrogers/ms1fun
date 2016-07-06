# Questions:
# 1. Calculation of the mass in PeakGroup - is it needed?
# 2. Let's say we ignore isotopes. But the results don't match. Is it a problem? Or have we just avoided false groups?
# 3. Testing: both for thesis & to make sure algorithm is bug-free. What's the best approach? One part would be
# profiling. Others - unit tests?
# 4. [M-C2H2]+2H[2C13] parent is [M-C2H2]+2H[C13]? What about [M-C2H2]+H[C13]? At the moment code relies on former.

import os, sys
sys.path.insert(0, '../code/')
from corr_cluster import Peak


class PeakGroup(object):
    def __init__(self):
        self.members = []   # Peak that fits the transformation, transformation itself and transformed mass
        self.vote = 0.0     # Total group vote
        self.M = 0.0        # Group mass

    def add_peak(self, peak, transformation):
        # Rather awkward way to check if the peak already exists in the group. If it does, and has a higher vote,
        # the lower one should be deleted. If it has a lower vote, there is no reason to add it to the group.
        for p, t, tm in self.members:
            if peak is p and transformation.vote > t.vote:
                self.members.remove((p, t, tm))
            elif peak is p and transformation.vote < t.vote:
                return

        # Reached either when there is no second peak, or the vote of the same peak's transformation is higher.
        self.members.append((peak, transformation, transformation.transform(peak)))
        self.recalculate()

    # Recalculates votes, intensity, mass. Used when members are added/removed from a group.
    def recalculate(self):
        total_intensity = 0.0   # Total intensity of all peaks in the group.
        self.M = 0.0            # Total mass of peaks in the group.
        self.vote = 0.0         # Total vote of all transformations in the group.

        # In case group does not have any members left (e.g. peak without a parent was removed leaving it empty)
        if len(self.members) <= 0:
            return

        # Calculates mass and vote of a group. Adapted from voter.py.
        for member in self.members:
            self.vote += member[1].vote                 # total vote += transformation.vote
            self.M += member[0].intensity * member[2]   # mass = peak's intensity * transformed peak
            total_intensity += member[0].intensity      # total intensity += peak intensity
        self.M /= total_intensity

    # For printing in command line. No other purpose, can be deleted later on to make the code tidier.
    def __str__(self):
        representation = "\nTotal vote: {} Mass: {}".format(self.vote, self.M)
        for peak, transformation, transformed_mass in sorted(self.members, key=lambda x: x[1].vote, reverse=True):
            representation += "\n   {:.4f}, {:.4f}, {:.2e}, {} ({:.4f}, {:.4f})".format(peak.mass, peak.rt,
                                peak.intensity, transformation, transformed_mass, transformation.vote)
        return representation


class IntensityClustering(object):
    def __init__(self, peak_file, transformations, rt_thresh=3, tolerance=10):
        self.peaks = []                 # Peaks extracted from a file
        self.peak_file = peak_file      # Data file with M/RT/intensity
        self.transformations = transformations  # Transformations we will be working with
        # These 2 would probably be better suited as parameters in functions. But more convenient as attributes for now.
        self.rt_thresh = rt_thresh      # RT threshold according to which peaks will be grouped
        self.tolerance = tolerance      # Tolerance in parts per million

        self.load_peaks(peak_file)      # Initialize peaks from given file

    # Produces a list of peaks from a given file. Numpy, etc. required because of corr_cluster.py code.
    def load_peaks(self, file_path):
        print ("\nCreating peaks from a file at: \n" + file_path + "...")

        if os.path.isfile(file_path):
            with open(file_path) as f:
                f.readline()  # Removes heading
                for line in f:
                    split = line.split(',')
                    peak_id = int(split[0])
                    peak_mass = float(split[1])     # Peak mass (or rather, m/z: mass to charge ratio)
                    peak_rt = float(split[2])       # Peak retention time
                    peak_int = float(split[3])      # Peak intensity
                    self.peaks.append(Peak(peak_id, peak_mass, peak_rt, peak_int))  # Creates peak with above fields
        else:
            print('ERROR! Given path does not lead to a file!')

        print ("\n{} peaks loaded!".format(str(len(self.peaks))))

    # Processes and groups the peaks around the most intense peaks.
    def heavylifting(self):
        groups = []     # Peak groups that are explained by highest voted m (transformed molecular mass)

        # For testing. m/z: 221.0420, RT: 611.6010,Intensity: 9.84e+05, M+K (182.0789,0.6)
        # most_intense = self.peaks[4316]

        # Sorts the peaks into separate lists, so that we don't need to sort them every time loop starts.
        intensity_sorted = sorted(self.peaks, key=lambda x: x.intensity, reverse=True)  # most-intense to least-intense
        rt_sorted = sorted(self.peaks, key=lambda x: x.rt)                              # lowest rt to highest rt

        # Main loop where grouping happens
        while len(intensity_sorted) > 0:
            most_intense = intensity_sorted[0]  # P:: Always the first one, as it's sorted in desc. order

            # Creates a set that will contain peaks that are within |n|s RT from most intense peak
            index = rt_sorted.index(most_intense)   # Finds most intense peak's index in rt_sorted list
            thresh_group = []                       # R:: peaks within rt threshold
            thresh_group = self.forward_pass(thresh_group, most_intense, index, rt_sorted)  # Go right (higher)
            index -= 1                              # Immutable, thus 'index' is not changed by the functions.
            thresh_group = self.backward_pass(thresh_group, most_intense, index, rt_sorted) # Go left (lower)

            # Finds possible molecular masses for the most intense peak by applying all transformations to the most
            # intense peak from the list. As it is the most intense (abundant) peak, we can assume it is mono-isotopic.
            # Thus, we do not need C13 and 2C13, which we can skip when determining what to add to the molecular masses.
            molecular_masses = []   # M:: list of molecular masses
            for transformation in self.transformations:
                if transformation.name[-4:-1] != 'C13':
                    molecular_masses.append(transformation.transform(most_intense))

            highest_voted = PeakGroup()     # Highest voted group of peaks

            # 1. For each molecular mass, calculates all possible peak masses. Records it, and the transformation.
            # 2. For each derived mass, checks if similar mass is in a group created around the most intense peak.
            # 3. If it is, stores it in a group of peaks that are explained by the transformations on M.
            # 4. Removes transformations that are missing parents from a group.
            # 5. Assesses whether the group that was just created has a higher vote than previous highest voted group.

            # At least one of the masses must be in that group, as one of the reverse transformations will inevitably
            # arrive at the same mass most intense peak was before transforming it. Other peaks from rt_thresh group
            # might be in there as well. This means that the initial molecule was transformed in different way, and thus
            # has a different mass when measured. To record these molecules, create an object that holds total vote
            # (from transformations) and mass, in addition to information on the peak and transformation used on it to
            # arrive at the initial mass (or to calculate end mass from initial mass).
            for mass in molecular_masses:
                reverse_masses = []   # N:: reverse transformation on molecular masses. Stores mass + transformation
                # Performs reverse transformations on a mass. Gives possible masses for the peak. Includes C13/2C13.
                for transformation in self.transformations:
                    peak_mass = transformation.reversetransform(mass)
                    reverse_masses.append([peak_mass, transformation])

                # For each peak mass in reverse masses, checks if there is anything similar in rt_thresh group.
                peak_group = PeakGroup()
                for peak_mass, transformation in reverse_masses:
                    for peak in thresh_group:
                        # 1000 000 * |(a-b)/a| <= tolerance
                        ppm = 1000000 * abs((peak.mass - peak_mass) / peak.mass)    # parts per million (ppm) difference
                        # If there is, adds it to the possible peak group, which will later be assessed on total vote.
                        if ppm <= self.tolerance:
                            peak_group.add_peak(peak, transformation)
                peak_group = self.remove_orphans(peak_group)
                # Only keeps the highest voted group out of all groups around various masses
                if peak_group.vote > highest_voted.vote:
                    highest_voted = peak_group

            # Stores the group with the highest vote for the current peak.
            groups.append(highest_voted)

            # Remove peaks that have been grouped (that are explained by m)
            grouped_peaks = [x[0] for x in highest_voted.members]   # List of peaks extracted from a tuple
            for peak in grouped_peaks:
                rt_sorted.remove(peak)
                intensity_sorted.remove(peak)

            # Command line: prints how many peaks are left to process. Delete later.
            # print ("{} peaks left to process...".format(len(rt_sorted)))

        return groups

    """                                 HELPER FUNCTIONS                                    """
    # Removes transformations whose parents are not in the group (orphans). Also applies to those that
    # have a parent, but the parent itself is an orphan. In such a case, both transformations are deleted.
    def remove_orphans(self, group):
        finished = False
        while not finished:
            to_remove = []
            for peak, transformation, transformed_mass in group.members:
                if transformation.parent is not None and transformation.parent not in [x[1] for x in group.members]:
                    to_remove.append((peak, transformation, transformed_mass))

            # Remove orphans
            for member in to_remove:
                group.members.remove(member)

            # If no orphans were removed, it means every transformation had its parent in a group of transformations.
            if len(to_remove) <= 0:
                finished = True
            # Only need to recalculate mass and votes if something was removed. Otherwise nothing changes.
            else:
                group.recalculate()
        return group

    # Finds the peaks that are [0-threshold] s higher than the target peak.
    def forward_pass(self, thresh_group, most_intense, index, rt_sorted):
        while True:
            if index >= len(rt_sorted):
                break
            if abs(most_intense.rt - rt_sorted[index].rt) > self.rt_thresh:
                break
            thresh_group.append(rt_sorted[index])
            index += 1
        return thresh_group

    # Finds the peaks that are [0-threshold] s lower than the target peak
    def backward_pass(self, thresh_group, most_intense, index, rt_sorted):
        while True:
            if index < 0:
                break
            if abs(most_intense.rt - rt_sorted[index].rt > self.rt_thresh):
                break
            thresh_group.append(rt_sorted[index])
            index -= 1
        return thresh_group

    """                         TEMPORARY: for testing                         """
    def print_peaks(self, peak_group=None):
        # Default will print peaks in an order they were retrieved from a file.
        # However, if some parameter is given, this will be skipped and other group of peaks printed.
        if peak_group is None:
            peak_group = self.peaks

        # In case list is empty
        if not peak_group:
            print ("List of peaks is empty!")
        else:
            for peak in peak_group:
                print ("Peak ID: {}  Mass: {:.4f}  RT: {:.4f}  Int: {:.4f}").format(peak.pid, peak.mass, peak.rt,
                                                                                    peak.intensity)
