import hash
import screed
import sys
import zlib
import cPickle
from screed import dna
from optparse import OptionParser

# TODO: encode these values into the bfs stored to disk
N_HT=4
HT_SIZES=100e6
k = 16
sizes = hash.get_n_primes_above_x(N_HT, int(HT_SIZES))

parser = OptionParser()

parser.add_option("-b", 
				  "--bloom_filters",
				  dest="bf_filenames",
				  help="Comma seperated list of archived bloom filters")

parser.add_option("-q", 
				  "--query_file",
				  dest="sample_filename",
				  help="Query FASTA or FASTQ file name")

parser.add_option("-m", 
				  "--max_miss_size",
				  dest="max_miss_size",
				  default=0,
				  help="Maximum size (bp) for a miss within a read")



(options, args) = parser.parse_args()


#bf_filenames = sys.argv[1]
#sample_filename = sys.argv[2]

bfs = []

for filename in options.bf_filenames.split(","):
	bf = hash.BloomFilter(sizes, k)
	# clear the Bloom filter table
	bf.tables = None

	# now load the stored Bloom filter tables into the bf.table attribute
	# for subsequent querying
	with open(filename, 'rb') as fp:
		data = zlib.decompress(fp.read())
		bf.tables = cPickle.loads(data)
	bfs.append(bf)

counts = {}

for read in screed.open(options.sample_filename):
	# error_state tracks the state of the read
	# 0 : no error yet
	# 1 : in the first error
	# 2 : past the first error
	# 3 : second error found
	error_state = [0] * len(bfs)
	# error_len tracks the length of the first error
	error_len = [0] * len(bfs)

	for i in range(len(read.sequence) - k + 1):
		kmer = read.sequence[i:i+k]
		for i in range(0,len(bfs)):
			# if the state is 3, then no need to check
			if error_state[i] != 3:
				bf = bfs[i]
				if (kmer in bf) or (dna.reverse_complement(kmer) in bf):
					# move out of the first error
					if error_state[i] == 1:
						error_state[i] = 2
				else:
					# first error, switch states and up the length
					if error_state[i] == 0:
						error_state[i] = 1
						error_len[i] += 1
					# continuation of the first error, up the length
					elif error_state[i] == 1:
						error_len[i] += 1
					#second error found
					elif error_state[i] == 2:
						error_state[i] = 3

	# a read will be put into was_in if there finishes in
	# state 0, or error_len is less thank k * options.max_miss_size
	# and is in either state 1 or 2
	was_in = []
	for i in range(0,len(bfs)):
		if error_state[i] == 0:
			was_in.append(i)
		if ( (error_state[i] == 1) or (error_state[i] == 2) ) and \
				(error_len[i] <= k * options.max_miss_size) : 
			was_in.append(i)


	if str(was_in) in counts:
		counts[str(was_in)] += 1
	else:
		counts[str(was_in)] = 1

sorted_keys=sorted(counts, key=counts.get, reverse=True)
for key in sorted_keys:
	print key,counts[key]
