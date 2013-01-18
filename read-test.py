import hash
import screed
import sys
import zlib
import cPickle
from screed import dna

N_HT=4
HT_SIZES=100e6

k = 16
sizes = hash.get_n_primes_above_x(N_HT, int(HT_SIZES))

bf_filenames = sys.argv[1]

sample_filename = sys.argv[2]

bfs = []

for filename in bf_filenames.split(","):
	bf = hash.BloomFilter(sizes, k)
	# clear the Bloom filter table
	bf.tables = None

	# now load the stored Bloom filter tables into the bf.table attribute
	# for subsequent querying
	with open(filename, 'rb') as fp:
		data = zlib.decompress(fp.read())
		bf.tables = cPickle.loads(data)
	bfs.append(bf)

#counts = {}

for read in screed.open(sample_filename):
	misses = 0
	miss_str = ''
	for i in range(len(read.sequence) - k + 1):
		kmer = read.sequence[i:i+k]
		#was_in = []
		for i in range(0,len(bfs)):
			bf = bfs[i]
			if kmer in bf:
				misses += 0
				miss_str += '1'
				#was_in.append(i)
			elif dna.reverse_complement(kmer) in bf:
				misses += 0
				miss_str += '1'
				#was_in.append(i)
			else:
				misses += 1
				miss_str += '0'
	if misses > 0 :
		print misses
		print miss_str
		
#		if str(was_in) in counts:
#			counts[str(was_in)] += 1
#		else:
#			counts[str(was_in)] = 1

#sorted_keys=sorted(counts, key=counts.get, reverse=True)
#for key in sorted_keys:
#	print key,counts[key]
