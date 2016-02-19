import numpy as np
import operator


def bwt(s):

    sd = list(s + "$")

    # form the first and last column of the cyclic permutation matrix
    strs = list()
    for k in range(len(sd)):
        sdp = np.roll(sd, k)
        strs.append((''.join(sdp), k))

    strs.sort(key=operator.itemgetter(0))
    fc = [s[0] for s,k in strs]
    lc = [s[-1] for s,k in strs]
    sarr = [k for s,k in strs]

    return fc,lc,sarr


def invert_bwt(lc):

    fc = sorted(lc)

    fc_cnt = dict()
    fc_cnt2index = dict()
    fc_index2cnt = list()
    for k,c in enumerate(fc):
        if c not in fc_cnt:
            fc_cnt[c] = 0
        fc_cnt[c] += 1
        fc_cnt2index[(c, fc_cnt[c])] = k
        fc_index2cnt.append(fc_cnt[c])
        
    lc_cnt = dict()
    lc_cnt2index = dict()
    lc_index2cnt = list()
    for k,c in enumerate(lc):
        if c not in lc_cnt:
            lc_cnt[c] = 0
        lc_cnt[c] += 1
        lc_cnt2index[(c, lc_cnt[c])] = k
        lc_index2cnt.append(lc_cnt[c])

    i = 0
    s = '$'
    while lc[i] != '$':
        l = lc[i]
        s = l + s
        cnt = lc_index2cnt[i]
        i = fc_cnt2index[(l, cnt)]

    return s


if __name__ == '__main__':

    # fc,lc,sarr = bwt('panamabananas')
    # istr = invert_bwt(lc)
    # print 'inverted: %s' % istr

    # fc,lc,sarr = bwt('banana')
    # print sarr

    # fc,lc,sarr = bwt('TCAGGGCTTG')
    # print lc

    # istr = invert_bwt('AT$AAACTTCG')
    # print 'inverted: %s' % istr

    print 363. / 6.
