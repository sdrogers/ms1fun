import sys, os, time
from sv_algorithm import Algorithm
sys.path.insert(0, os.path.normpath(os.getcwd() + os.sep + os.pardir) + '/code/')
import transformation

# The path where data files live: currently, it is in grandparent_folder/data/<file>.csv
# Or in other words, to navigate from the folder where the code is in: >>> cd ../../data.
# In theory, should work on Unix/Windows, as long as files are in the right folders.
data_dir = os.path.normpath(os.getcwd() + os.sep + os.pardir + os.sep + os.pardir) + '/data/'
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
cluster = Algorithm(file_path, transformations, rt_thresh=2, tolerance=10)

# Taking most intense peak, figuring out the groups around it (within given RT window)
print ("Processing...")
start_time = time.clock()
groups = cluster.heavylifting()
end_time = time.clock()
print ("\nData processing time: " + str(end_time-start_time))

# Write out the results to a file
groups = sorted(groups, key=lambda x: x.vote, reverse=True)  # top voted at the top
output_dir = os.path.normpath(os.getcwd() + os.sep + os.pardir + os.sep + os.pardir) + '/output/'
output_file = 'Urine_37_fullscan1_POS_by_vote.txt'
output_path = output_dir + output_file

with open(output_path, 'w') as f:
    for group in groups:
        line = "vote: {}, M: {}\n".format(group.vote, group.mass)
        f.write(line)
        line = '\tPeak m/z,Peak rt,Peak intensity,transformation (transformed mass,vote)\n'
        f.write(line)
        # members = (peak, transformation, transformed_mass); x[1] = transformation
        for (peak, transformation, transformed_mass) in sorted(group.members, key=lambda x: x[1].vote, reverse=True):
            line = "\t{:.4f},{:.4f},{:.2e},{} ({:.4f},{:.4f})\n".format(peak.mass, peak.rt, peak.intensity,
                    transformation, transformed_mass, transformation.vote)
            f.write(line)
        f.write('\n')

print_time = time.clock()
print ("Output file has been prepared! Total execution time: " + str(print_time-start_time))
print ("Total groups found: " + str(len(groups)))

# for group in groups:
#     print group
