from copy import deepcopy

import numpy as np


def entropy(p):

    nz = p > 0
    return -np.sum(p[nz] * np.log2(p[nz]))


def hamming(s1, s2):
    return np.sum([c1 != c2 for c1,c2 in zip(s1, s2)])


def generate_kmers(k, pos=0, s='AAAA', alpha=['A', 'G', 'C', 'T']):

    assert len(s) == k

    kmers = list()
    kmers.append(s)
    for a in alpha:
        snew = s[:pos] + a + s[(pos+1):]
        kmers.append(snew)
        if pos < k-1:
            kmers.extend(generate_kmers(k, pos=pos+1, s=snew))
    return np.unique(kmers)


def min_hamming(pattern, str):
    k = len(pattern)
    min_dist = np.inf
    min_str = None
    for n in range(len(str)-k):
        s = str[n:(n+k)]
        d = hamming(pattern, s)
        if d < min_dist:
            min_dist = d
            min_str = s

    return min_dist,min_str


def d(pattern, strs):
    return np.sum([min_hamming(pattern, s)[0] for s in strs])


def median_string(strs, k, kmers=None):
    min_dist = np.inf
    min_str = None

    display_score = True
    if kmers is None:
        display_score = False
        kmers = generate_kmers(k, s='A'*k)

    for kmer in kmers:
        dist = d(kmer, strs)
        if display_score:
            print '%s: %d' % (kmer, dist)
        if dist < min_dist:
            min_dist = dist
            min_str = kmer

    return min_str


def profile_prob(P, s, alpha=['A', 'C', 'G', 'T']):

    p = 1.
    for k,c in enumerate(s):
        i = alpha.index(c)
        p *= P[i, k]
    return p


if __name__ == '__main__':

    # print entropy(np.array([0.5, 0, 0, 0.5]))
    # print entropy(np.array([0.25, 0.25, 0.25, 0.25]))
    # print entropy(np.array([0., 0, 0, 1]))
    # print entropy(np.array([0.25, 0, 0.5, 0.25]))

    strs = ['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
            'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
            'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG']

    # kmers = generate_kmers(3, s='AAA')
    # print 'len(kmers)=%d, expected=%d' % (len(kmers), 4**3)
    # med_str = median_string(strs, 7)
    # print 'Median String: %s' % med_str

    med_str = median_string(strs, 7, kmers=['AATCCTA', 'TCTGAAG', 'CGTGTAA', 'GAACCAC', 'GTCAGCG', 'GATGAGT'])

    P = np.array([ [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
                   [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
                   [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
                   [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
                 ])
    print 'Profile p=%0.9f' % profile_prob(P, 'TCGGTA')
    print 'Profile p=%0.9f' % profile_prob(P, 'AAGCGA')
