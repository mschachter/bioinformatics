import numpy as np


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


def hamming(s1, s2):
    return np.sum([c1 != c2 for c1,c2 in zip(s1, s2)])


def entropy(p):

    nz = p > 0
    return -np.sum(p[nz] * np.log2(p[nz]))


def profile_prob(P, s, alpha=['A', 'C', 'G', 'T']):

    p = 1.
    for k,c in enumerate(s):
        i = alpha.index(c)
        p *= P[i, k]
    return p


def make_profile(motifs, alpha=['A', 'C', 'G', 'T'], normalize=True, laplace=True):

    N = len(motifs)
    M = len(motifs[0])
    Pcnt = np.zeros([len(alpha), M])

    for m in range(M):
        for n in range(N):
            i = alpha.index(motifs[n][m])
            Pcnt[i, m] += 1

    if laplace:
        Pcnt += 1

    if normalize:
        Pcnt /= Pcnt.sum(axis=0)

    return Pcnt


def profile_prob_kmer(P, s, k):

    best_p = 0.
    best_kmer = None
    for i in range(len(s)-k):
        kmer = s[i:(i+k)]
        p = profile_prob(P, kmer)
        if p > best_p:
            best_p = p
            best_kmer = kmer

    return best_kmer


def motif_score(motifs, alpha=['A', 'C', 'T', 'G']):

    M = len(motifs[0])
    score = list()
    P = make_profile(motifs, normalize=False)

    for m in range(M):
        imax = P[:, m].argmax()
        scorem = np.sum([P[i, m] for i in range(len(alpha)) if i != imax])
        score.append(scorem)

    return np.sum(score)


def greedy_motif_search(dna, k):

    N = len(dna)

    # get initial guess for motif matrix from first kmer of each DNA string
    best_motifs = list()
    for n in range(N):
        best_motifs.append(dna[n][:k])

    for m in range(len(dna[0])-k):
        kmer = dna[0][m:(m+k)]
        motifs = [kmer,]
        for n2 in range(1, N):
            P = make_profile(motifs)
            motif = profile_prob_kmer(P, dna[n2], k)
            motifs.append(motif)

        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs

    return best_motifs


def randomized_motif_search(dna, k, maxiters=10):

    N = len(dna)
    M = len(dna[0])
    best_motifs = list()
    for n in range(N):
        i = np.random.randint(M-k)
        best_motifs.append(dna[n][i:(i+k)])

    motifs = best_motifs
    Pold = np.zeros([4, k])
    for iter in range(maxiters):
        P = make_profile(motifs)
        motifs = [profile_prob_kmer(P, s, k) for s in dna]
        Pdiff = np.sum((P - Pold)**2) / (4*k)
        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs
        Pold = P
        if Pdiff < 1e-6:
            break

    return best_motifs


def random_motif(P, k, alpha=['A', 'C', 'T', 'G']):

    s = ''
    Pcsum = np.cumsum(P, axis=0)
    for i in range(k):
        F = Pcsum[:, i]
        r = np.random.rand()
        j = np.where(F >= r)[0]
        ii = np.min(j)
        s += alpha[ii]
    return s


def gibbs_motif_search(dna, k, maxiters=10):

    N = len(dna)
    M = len(dna[0])
    best_motifs = list()
    for n in range(N):
        i = np.random.randint(M-k)
        best_motifs.append(dna[n][i:(i+k)])

    motifs = best_motifs
    Pold = np.zeros([4, k])
    for iter in range(maxiters):
        P = make_profile(motifs)
        motifs = [profile_prob_kmer(P, s, k) for s in dna]

        # randomly generate one of the motifs
        ri = np.random.randint(N)
        motifs[ri] = random_motif(P, k)

        Pdiff = np.sum((P - Pold)**2) / (4*k)
        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs
        Pold = P
        if Pdiff < 1e-6:
            break

    return best_motifs


def generate_random_dna(k, alpha=['A','C', 'T', 'G']):
    alpha = np.array(alpha)
    return list(alpha[np.random.randint(len(alpha), size=k)])


if __name__ == '__main__':

    motifs = [ ['A', 'A', 'T', 'G', 'T', 'G', 'G'],
               ['T', 'G', 'T', 'G', 'C', 'C', 'G'],
               ['A', 'G', 'T', 'G', 'G', 'C', 'A'],
               ['A', 'A', 'T', 'G', 'T', 'C', 'A']
             ]

    dna = ['CAAGCGTGCAGAAAACGGCC',
           'CCCCGAAGCCGCAATACCCC',
           'TAAAGAAACCAAGCAGCGTT',
           'ACATCAAGCTCAAAACCCAG',
           'TGCAAGCCAATTCCCAAGCA']

    P = np.array([[0.2, 0.5, 1.0, 0.1],
                  [0.3, 0.2, 0.,  0.3],
                  [0.4, 0.2, 0.,  0.4],
                  [0.1, 0.1, 0.,  0.2]
                 ])

    # best_motifs = greedy_motif_search(dna, 4)
    # best_motifs = randomized_motif_search(dna, 4)
    best_motifs = gibbs_motif_search(dna, 4)

    # print 'Best Motifs:'
    # print best_motifs

    dna = ['TGACGTTC',
           'TAAGAGTT',
           'GGACGAAA',
           'CTGTTCGC']

    motifs = ['TGA', 'GTT', 'GAA', 'TGT']

    P = make_profile(motifs, laplace=True)
    new_motifs = [profile_prob_kmer(P, s, 3) for s in dna]

    print P
    print new_motifs









