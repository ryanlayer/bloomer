import string
import numpy
from bitarray import bitarray

MAX_K=16
DEFAULT_K=8

def is_prime(n):
   '''
   checks if a number is prime
   '''
   if n < 2:
      return False
   if n == 2:
      return True
   for x in range(2, int(n**0.5)+1, 2):
      if n % x == 0:
         return False
   return True

def get_n_primes_above_x(n, x):
   '''
   steps forward until n primes (other than 2) have been
   found that are smaller than x.
   '''
   primes = []
   i = x+1
   if i % 2 == 0:
      i += 1
   while len(primes) != n and i > 0:
      if is_prime(i):
         primes.append(i)
      i += 2
   return primes

def hash(word):
    assert len(word) <= MAX_K

    value = 0
    for n, ch in enumerate(word):
        value += ord(ch) * 128**n

    return value

class BloomFilter(object):
    allchars = "".join([ chr(i) for i in range(128) ])
    
    def __init__(self, tablesizes, k=DEFAULT_K):
        #self.tables = [ (size, numpy.zeros(size, dtype=numpy.bool)) \
        #                 for size in tablesizes ]

        self.tables = [ (size, bitarray(size)) \
                         for size in tablesizes ]
        for table in self.tables:
            table[1].setall(False)

        self.k = k

    def add(self, word):
        val = hash(word)
        for size, ht in self.tables:
            ht[val % size] = 1

    def __contains__(self, word):
        val = hash(word)
        return all( ht[val % size] for (size, ht) in self.tables )

    def insert_text(self, text):
        for i in range(len(text) - self.k + 1):
            self.add(text[i:i+self.k])

    def occupancy(self):
        return [ sum(t)/float(len(t)) for _, t in self.tables ]

def first_next_word(bf, word):
    prefix = word[1:]
    for ch in bf.allchars:
        word = prefix + ch
        if word in bf:
            return ch, word
        
    return None, None

def retrieve_first_sentence(bf, start):
    word = start[-bf.k:]

    while 1:
        ch, word = first_next_word(bf, word)
        if ch is None:
            break

        start += ch

    return start

def next_words(bf, word):
    prefix = word[1:]
    for ch in bf.allchars:
        word = prefix + ch
        if word in bf:
            yield ch
            
def previous_words(bf, word):
    suffix = word[:-1]
    for ch in bf.allchars:
        word = ch + suffix
        if word in bf:
            yield ch
        

def retrieve_all_sentences(bf, start, max_depth=50):
    if max_depth == 0:
        yield start + '...'
        return
        
    word = start[-bf.k:]

    n = -1
    for n, ch in enumerate(next_words(bf, word)):
        for sentence in retrieve_all_sentences(bf, start + ch, max_depth-1):
            yield sentence

    if n < 0:
        yield start

def count_connected_graph(bf, word, cutoff=10, keeper=None):
    assert len(word) == bf.k
    if keeper is None:
        keeper = set()
    
    if word in keeper:
        return len(keeper)

    keeper.add(word)
    if len(keeper) >= cutoff:
        return len(keeper)
    
    for n, ch in enumerate(next_words(bf, word)):
        count_connected_graph(bf, word[1:] + ch, cutoff=cutoff, keeper=keeper)

    for n, ch in enumerate(previous_words(bf, word)):
        count_connected_graph(bf, ch + word[:-1], cutoff=cutoff, keeper=keeper)
        
    return len(keeper)
