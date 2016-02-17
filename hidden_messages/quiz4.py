

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


def make_profile(motifs, alpha=['A', 'C', 'G', 'T']):

    N = len(motifs)
    M = len(motifs[0])
    Pcnt = np.zeros([len(alpha), M])

    for m in range(M):
        for n in range(N):
            i = alpha.index(motifs[n][m])
            Pcnt[i, m] += 1

    P = Pcnt / float(len(alpha))

    return P


def profile_prob_kmer(P, s, k):

    best_p = 0.
    best_kmer = None
    for i in range(len(s)-k):
        kmer = s[i:(i+k)]
        p = profile_prob(P, s)
        if p > best_p:
            best_p = p
            best_kmer = kmer

    return best_kmer


def motif_score(motifs):

    N = len(motifs)
    M = len(motifs[0])

    for m in range(M):
        for n in range(N):





def greedy_motif_search(dna, k):

    N = len(dna)

    # get initial guess for motif matrix from first kmer of each DNA string
    best_motifs = list()
    for n in range(N):
        best_motifs.append(dna[n][:k])

    for m in range(len(dna[0])-k):
        kmer = dna[0][m:(m+k)]
        motifs = [kmer]
        for n2 in range(1, N):
            P = make_profile(motifs)
            motif = profile_prob_kmer(P, dna[n2])
            motifs.append(motif)

        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs

    return best_motifs









