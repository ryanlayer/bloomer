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

counts = [0] * len(bfs)

for read in screed.open(sample_filename):
	for i in range(len(read.sequence) - k + 1):
		kmer = read.sequence[i:i+k]
		for i in range(0,len(bfs)):
			bf = bfs[i]
			if kmer in bf:
				counts[i] += 1

i = 0
for filename in bf_filenames.split(","):
	print filename,counts[i]
	i+=1
