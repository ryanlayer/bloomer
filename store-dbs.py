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

for filename in filenames.split(","):
	bf = hash.BloomFilter(sizes, k)
	for read in screed.open(filename):
		# adds each kmer in the sequence to the database
		bf.insert_text(read.sequence)


	# compress and store the Bloom filter table as a file.
	archive_filename = filename + ".archive.gz" 
	with open(archive_filename, 'wb') as fp:
		fp.write(zlib.compress(cPickle.dumps(bf.tables,\
							   cPickle.HIGHEST_PROTOCOL),9))
	fp.close()
