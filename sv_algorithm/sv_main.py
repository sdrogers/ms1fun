"""
Code that launches the algorithm.
"""

import sys, os, time
from sv_algorithm import IntensityClustering
# import cProfile, pstats, StringIO
sys.path.insert(0, '../code/')
import transformation
from sv_tests import singleton_test
from helper_functions import hmdb_analysis

# Profiling stuff
# pr = cProfile.Profile()
# pr.enable()

# The path where data files live: currently, it is in grandparent_folder/data/<file>.csv
grandparent_dir = os.path.normpath(os.getcwd() + os.sep + os.pardir + os.sep + os.pardir)
data_dir = grandparent_dir + '/data/'
data_file = 'Urine_37_fullscan1_POS.csv'
file_path = data_dir + data_file

# Path where transformation files live
transformation_folder = os.path.normpath(os.getcwd() + os.sep + os.pardir) + '/dbs/'
transformation_file = 'pos_transformations_camra.yml'
transformation_path = transformation_folder + transformation_file

print ("\nLoading transformations from:\n" + transformation_path + "...\n")
transformations = transformation.load_from_file(transformation_path)
print ("\nLoaded " + str(len(transformations)) + " transformations.")

# Initialises an object that loads peaks from a file and gives the RT threshold in seconds
cluster = IntensityClustering(file_path, transformations, rt_thresh=2, tolerance=10)

# Taking most intense peak, figuring out the groups around it (within given RT window)
print ("Processing...")
start_time = time.clock()
groups = cluster.heavylifting()
end_time = time.clock()
print ("Data processing time: " + str(end_time-start_time))

# Write out the results to a file. Top voted group presented at the top.
groups = sorted(groups, key=lambda x: x.vote, reverse=True)
output_dir = grandparent_dir + '/output/'
output_file = 'Urine_37_fullscan1_POS_by_vote.txt'
output_path = output_dir + output_file

with open(output_path, 'w') as f:
    for group in groups:
        line = "vote: {}, M: {}\n".format(group.vote, group.M)
        f.write(line)
        line = '\tPeak m/z,Peak rt,Peak intensity,transformation (transformed mass,vote)\n'
        f.write(line)
        # members = (peak, transformation, transformed_mass); x[1] = list(transformation)
        for (peak, transformation, transformed_mass) in sorted(group.members, key=lambda x: x[1].vote, reverse=True):
            line = "\t{:.4f},{:.4f},{:.2e},{} ({:.4f},{:.4f})\n".format(peak.mass, peak.rt, peak.intensity,
                    transformation, transformed_mass, transformation.vote)
            f.write(line)
        f.write('\n')

print ("\nOutput file has been prepared:\n{}".format(output_path))
print ("Total groups found: " + str(len(groups)))

# Analyse the output, whether it makes sense.
# counts = hmdb_analysis(groups, filename=output_file[:-4])

# Command line printing
for group in groups:
    print group

# Test specific parts of an algorithm.
print("Singleton groups test. Every singleton group contains M+H: {}"
        .format(singleton_test(groups)))

# Profiling continued
# pr.disable()
# s = StringIO.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print s.getvalue()

