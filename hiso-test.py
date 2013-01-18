import hash
import screed
import sys
import zlib
import cPickle
from screed import dna
from optparse import OptionParser

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

parser.add_option("-n", 
				  "--number_of_buckets",
				  type="int",
				  dest="num_buckets",
				  default=10,
				  help="Number of buckets in histograph")


(options, args) = parser.parse_args()

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

buckets = [[0]*options.num_buckets for x in range(0,len(bfs))]
step = float(1)/float(options.num_buckets);

for read in screed.open(options.sample_filename):
	kmers_found = [0] * len(bfs)
	num_kmers = len(read.sequence) - k + 1
	for i in range(num_kmers):
		kmer = read.sequence[i:i+k]
		for i in range(0,len(bfs)):
			bf = bfs[i]
			if kmer in bf or dna.reverse_complement(kmer) in bf:
				kmers_found[i] += 1
	print kmers_found,
	for i in range(len(bfs)):
		bucket = \
			min(9,int((float(kmers_found[i])/float(num_kmers))/float(step)))
		print bucket,
		buckets[i][bucket] += 1	
	print

for bucket in buckets:
	print "\t".join(map(str, bucket))

		#if str(was_in) in counts:
			#counts[str(was_in)] += 1
		#else:
			#counts[str(was_in)] = 1

#sorted_keys=sorted(counts, key=counts.get, reverse=True)
#for key in sorted_keys:
#	print key,counts[key]
