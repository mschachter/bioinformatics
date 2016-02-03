import operator

def Count(s, pat):

    cnt = 0
    N = len(pat)
    for k in range(len(s)-N):
        if s[k:(k+N)] == pat:
            cnt += 1

    return cnt

def Kmer(s, n):

    kmers = dict()
    for k in range(len(s)-n):
        w = s[k:(k+n)]
        if w not in kmers:
            kmers[w] = 0
        kmers[w] += 1

    lst = [(k,c) for k,c in kmers.items()]
    lst.sort(key=operator.itemgetter(1), reverse=True)

    for k,c in lst:
        print '%s: %d' % (k, c)


if __name__ == '__main__':
    pass


