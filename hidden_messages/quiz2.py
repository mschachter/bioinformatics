import numpy as np
import matplotlib.pyplot as plt

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def hamming(s1, s2):
    return np.sum([c1 != c2 for c1,c2 in zip(s1, s2)])

def skew(s, i):
    s = s.upper()
    num_c = np.sum([x == 'C' for x in s[i:]])
    num_g = np.sum([x == 'G' for x in s[i:]])
    return num_c - num_g

def all_skew(s):
    the_skew = list()
    for i in range(len(s)):
        the_skew.append(skew(s, i))
    return np.array(the_skew)

def Count(s, p, d):

    cnt = 0
    for k in range(len(s)-len(p)+1):
        s1 = s[k:(k+len(p))]
        hd = hamming(s1, p)
        if hd <= d:
            print 'Found approx match: s1=%s, p=%s, hd=%d' % (s1, p, hd)
            cnt += 1
    return cnt

def dNeighborhood(s, d):

    K = len(s)

    n = 0
    for k in range(K-1):
        n += choose(K, k+1) * 3

    return n


if __name__ == '__main__':

    """
    d = hamming('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 'CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG')
    print 'distance=%d' % d

    sk = all_skew('CATTCCAGTACTTCGATGATGGCGTGAAGA')
    print 'min skew=%d' % sk.argmin()
    plt.plot(sk)
    plt.show()
    """

    # ac = Count('CATGCCATTCGCATTGTCCCAGTGA', 'CCC', 2)
    # print 'Approx count=%d' % ac

    print '# of elements in k-neighborhood: %d' % dNeighborhood('ACGT', 3)
