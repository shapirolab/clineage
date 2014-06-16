import random
from unittest import TestCase
from Bio import SeqIO
from random import randint, sample
from linapp.models import Sequence
import re
import hashlib
__author__ = 'Yair'


class SeqsUtil:

    def IsSimpleSeq(self,seq):
        return re.match('^[ACTGactg]+$',seq)

    def IsDegenerateSeq(self,seq):
        return re.match('^[ACGTRYKMSWBDHVN]+$',seq, re.IGNORECASE)

    def RandomDNASeq(self,len):
        '''
        Genrates a simple random DNA sequence of the requested length
        '''
        return ''.join([random.choice('ACTG') for i in range (len)])

    def MergeSeqs(self,s1,s2):
        minOlap = 2
        if not s1 or not s2:
            return ''.join(s1,s2)
        for i in range(minOlap,min(len(s1),len(s2))):
            if s1[-i:] != s2[0:i]:
                if i == minOlap:
                    raise Exception ( "No minimal match")
                return s1 + s2[i-1:]

    def SeqDotPlot(self,s1,s2,window = 7):
        dict_one = {}
        dict_two = {}
        for (seq, section_dict) in [(s1.upper(), dict_one),
            (s2.upper(), dict_two)]:
            for i in range(len(seq)-window):
                section = seq[i:i+window]
                try:
                    section_dict[section].append(i)
                except KeyError:
                    section_dict[section] = [i]
        #Now find any sub-sequences found in both sequences
        #(Python 2.3 would require slightly different code here)
        matches = set(dict_one).intersection(dict_two)
        print "%i unique matches" % len(matches)
        #Create lists of x and y co-ordinates for scatter plot
        x = []
        y = []
        for section in matches:
            for i in dict_one[section]:
                for j in dict_two[section]:
                    x.append(i)
                    y.append(j)
            # Draw the dot plot as scattered plot
        import pylab
        pylab.cla() #clear any prior graph
        pylab.gray()
        pylab.scatter(x,y)
        pylab.xlim(0, len(s1)-window)
        pylab.ylim(0, len(s2)-window)
        pylab.xlabel("%s (length %i bp)" % ("Sequence1", len(s1)))
        pylab.ylabel("%s (length %i bp)" % ("Sequecne2", len(s2)))
        pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
        pylab.show()

class SeqsUtilTest (TestCase):

    def setUp(self):
        self.su = SeqsUtil()

    def test_IsSimpleSeq(self):
        s1  = 'ACTGCTGC'
        s2 = ''
        s3 = 'CGTCD'
        self.su.IsSimpleSeq(s1)
        self.assert_(self.su.IsSimpleSeq(s1))
        self.assertFalse(self.su.IsSimpleSeq(s2))
        self.assertFalse(self.su.IsSimpleSeq(s3))


    def test_MergeSeqs(self):
        """
        Test whether 2 sequences are merged to one correctly
        """
        s1 = 'ACCTG'
        s2 = 'TGAAA'
        s3 = 'TAGCT'
        s4 =''
        su = SeqsUtil()
        merged = self.su.MergeSeqs(s1,s2)
        self.assertEqual(merged,'ACCTGAAA')
        self.assertRaises(Exception,self.su.MergeSeqs,s1,s3)
        self.assertRaises(Exception,self.su.MergeSeqs,s1,s4)

    def test_SeqDotPlot(self):
        s1 = self.su.RandomDNASeq(100)
        print s1
        s2 = s1[0:30] + self.su.RandomDNASeq(70)
        self.su.SeqDotPlot(s1,s2,6)


num2base = {0:'A', 1:'C', 2:'G', 3:'T'}
basecomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'t', 't':'a', 'c':'g', 'g':'c'}
bases = set('ACTG')
purines = set('AG')
basestring = 'ACTG'

#import string
#_translation_table = string.maketrans('TAGCtagc', 'ATCGATCG')
##This will be the point we should change when we decide to preserve sequences casing. TODO: decide global casing policy.
#def complement(strand):
#    """fast, python translation table based, implementation of complement"""
#    return strand.translate(_translation_table)

def randseq(length):
    seq = ''
    for i in range(length):
        seq += num2base[randint(0, 3)]
    return seq

def complement(seq):
    assert re.match('^[ACTGactg]+$',seq)
    compseq = ''
    for base in seq:
        compseq += basecomp[base]
    return compseq

def mutate(seq,i,base): #index is zero based!
    assert re.match('^[ACTGactg]+$',seq)
    return seq[:i] + base + seq[i+1:]

def regionalmutation(seq, i, j, percentage): #randomize region i to j according percentage (given as fraction)
    assert re.match('^[ACTGactg]+$',seq)
    if i>j: #correction for bad input where indices are switched
        i, j = j, i
    mutations = sample(range(i,j), int(round((j-i)*percentage)))
    for m in mutations:
        seq = mutate(seq, m, sample(list(bases - set(seq[m])), 1)[0])
    return seq

def save_sequence(sequence):
    seq = sequence.strip().upper()
    try:
        sequence = Sequence.objects.get(hash=hashlib.md5(seq).hexdigest())
    except Sequence.DoesNotExist:
        sequence = Sequence(length=len(seq), sequence=seq, hash=hashlib.md5(seq).hexdigest())
        sequence.save()
    return sequence