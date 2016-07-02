# Questions:
# 1. Total vote is very low (0.32967032967 compared to 6.765 in ..._by_vote.txt). Changes to transformations' files?
# 2. Calculation of the mass in PeakGroup - is it needed?

import os, sys
sys.path.insert(0, os.path.normpath(os.getcwd() + os.sep + os.pardir) + '/code/')
from corr_cluster import Peak


class PeakGroup(object):
    def __init__(self):
        self.members = []   # Peak that fits the transformation, transformation itself and transformed mass
        self.vote = 0.0     # Total group vote
        self.mass = 0.0     # Group mass

    def add_peak(self, peak, transformation):
        # Rather awkward way to check if the peak already exists in the group. If it does, and has a higher vote,
        # the lower one should be deleted. If it has a lower vote, there is no reason to add it to the group.
        for p, t, tm in self.members:
            if peak is p and transformation.vote > t.vote:
                self.members.remove((p, t, tm))
            elif peak is p and transformation.vote < t.vote:
                return

        # This is reached either when the vote of the same peak's transformation is higher, or there is no second peak.
        self.members.append((peak, transformation, transformation.transform(peak)))
        self.vote += transformation.vote

        # Calculates mass of a group. Adapted from voter.py.
        total_intensity = 0.0       # Total intensity of all peaks in the group.
        self.mass = 0.0             # To reset mass each time the peak is added.
        for member in self.members:
            self.mass += member[0].intensity * member[2]    # Peak's intensity * transformed peak
            total_intensity += member[0].intensity          # total intensity += peak intensity
        self.mass /= total_intensity

    # For console printing. No other purpose, can be deleted later on to make the code tidier.
    def __str__(self):
        representation = "\nTotal vote: {} Mass: {}".format(self.vote, self.mass)
        for peak, transformation, transformed_mass in sorted(self.members, key=lambda x: x[1].vote, reverse=True):
            representation += "\n   {:.4f}, {:.4f}, {:.2e}, {} ({:.4f}, {:.4f})".format(peak.mass, peak.rt,
                                peak.intensity, transformation, transformed_mass, transformation.vote)
        return representation


class Algorithm(object):
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
                heading = f.readline()  # Removes heading
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

    # Processes and groups the peaks. Currently takes ~614-616s to process 7733 peaks.
    def heavylifting(self):
        groups = []     # Peak groups that are explained by highest voted m (transformed molecular mass)

        """            1st step of the algorithm           """
        # For testing. m/z: 221.0420, RT: 611.6010,Intensity: 9.84e+05, M+K (182.0789,0.6)
        # most_intense = self.peaks[4316]
        intensity_sorted = sorted(self.peaks, key=lambda x: x.intensity, reverse=True)  # most-intense to least-intense
        rt_sorted = sorted(self.peaks, key=lambda x: x.rt)                              # least rt to most rt

        # Main loop where grouping happens
        print "Processing..."
        while len(intensity_sorted) > 0:
            """            1st step of the algorithm           """
            most_intense = intensity_sorted[0]  # P:: Always the first one, as it's sorted in desc. order

            """            2nd step of the algorithm           """
            # Creates a set that will contain peaks that are within |n|s RT from most intense peak
            index = rt_sorted.index(most_intense)   # Finds most intense peak's index in rt_sorted list
            thresh_group = []                       # R:: peaks within rt threshold
            thresh_group = self.forward_pass(thresh_group, most_intense, index, rt_sorted)  # Go right (higher)
            index -= 1
            thresh_group = self.backward_pass(thresh_group, most_intense, index, rt_sorted) # Go left (lower)

            """            3rd step of the algorithm           """
            # Finds the possible molecular masses for the most intense peak. Found by applying all transformation to
            # the most intense peak from the list.
            molecular_masses = []   # M:: list of molecular masses
            for transformation in self.transformations:
                molecular_masses.append(transformation.transform(most_intense))

            """            4th step of the algorithm: a           """
            highest_voted = PeakGroup()     # Highest voted group of peaks

            # 1. For each molecular mass, calculates all possible peak masses.
            # 2. Records the derived mass and the transformation that was used to arrive at that mass.
            # 3. For each derived mass, checks if such a mass is in a group created around the most intense peak.
            # 4. If it is, stores it in a group of peaks that are explained by the transformations on M.
            # At least one of the masses must be in that group, as one of the reverse transformations will inevitably
            # arrive at the same mass most intense peak was before transforming it. Other peaks from rt_thresh group
            # might be in there as well. This means that the initial molecule was transformed in different way, and thus
            # has a different mass when measured. We want to record all such molecules. For this, create an object that
            # holds total vote (from transformations) and mass, in addition to information on the peak and transformas
            # that can be used to arrive at the initial mass (or to calculate end mass from initial mass).
            for mass in molecular_masses:
                # Performs reverse transformations on a given mass. Gives possible masses for the peak.
                reverse_masses = []   # N:: reverse transformation on molecular masses. Stores mass + transformation
                for transformation in self.transformations:
                    peak_mass = transformation.reversetransform(mass)
                    reverse_masses.append([peak_mass, transformation])

                """            4th step of the algorithm: c           """
                # For each peak mass in reverse masses...
                peak_group = PeakGroup()
                for peak_mass, transformation in reverse_masses:
                    for peak in thresh_group:
                        # 1000 000 * |(a-b)/a| <= tolerance
                        ppm = 1000000 * abs((peak.mass - peak_mass) / peak.mass)    # ppm difference
                        if ppm <= self.tolerance:
                            peak_group.add_peak(peak, transformation)

                """            5th step of the algorithm           """
                # Only keeps the highest voted group
                if peak_group.vote > highest_voted.vote:
                    highest_voted = peak_group

            """            6th step of the algorithm           """
            # Stores the group with the highest vote in a list of groups
            groups.append(highest_voted)

            """            7th step of the algorithm           """
            # Remove peaks that have been grouped (that are explained by m)
            grouped_peaks = [x[0] for x in highest_voted.members]   # List of peaks extracted from a tuple
            for peak in grouped_peaks:
                rt_sorted.remove(peak)
                intensity_sorted.remove(peak)
            # print ("\n{} peaks left to process...".format(len(rt_sorted)))

        """         The End         """
        return groups



    """                                 HELPER FUNCTIONS                                    """
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

    """                         FUNCTIONS FOR QUICK & DIRTY TESTING                         """
    # Sort on intensity - least to most
    def intensity_sort(self):
        int_sort = sorted(self.peaks, key=lambda x: x.intensity, reverse=False)
        return int_sort

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
