import hash
import screed
import sys
import zlib
import cPickle

N_HT=4
HT_SIZES=100e6

print 'creating data structure'

k = 16
sizes = hash.get_n_primes_above_x(N_HT, int(HT_SIZES))

print 'loading fastq into Bloom filter'

filenames = sys.argv[1]

bfs = []

for filename in filenames.split(","):
	bf = hash.BloomFilter(sizes, k)
	# clear the Bloom filter table
	bf.tables = None

	# now load the stored Bloom filter tables into the bf.table attribute
	# for subsequent querying
	with open(filename, 'rb') as fp:
		data = zlib.decompress(fp.read())
		bf.tables = cPickle.loads(data)
	bfs.append(bf)
