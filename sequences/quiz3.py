import numpy as np


aa_mass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115,
           'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}


aa_short2long = {'G':'Gly', 'V':'Val', 'Y':'Tyr', 'S':'Ser', 'C':'Cys', 'W':'Trp', 'F':'Phe',
                 'L':'Leu', 'N':'Asn', 'K':'Lys', 'T':'Thr', 'R':'Arg', 'I':'Ile',
                 'M':'Met', 'H':'His', 'Q':'Gln', 'P':'Pro', 'A':'Ala', 'E':'Glu'}


aa2codons = {'M':['AUG'],
             'I':['AUA', 'AUC', 'AUU'],
             'R':['AGG', 'AGA', 'CGG', 'CGA', 'CGC', 'CGU'],
             'S':['AGC', 'AGU', 'UCG', 'UCA', 'UCC', 'UCU'],
             'T':['ACG', 'ACA', 'ACC', 'ACU'],
             'K':['AAG', 'AAA'],
             'N':['AAC', 'AAU'],
             'L':['UUG', 'UUA', 'CUG', 'CUA', 'CUC', 'CUU'],
             'F':['UUC', 'UUU'],
             'W':['UGG'],
             '*':['UGA', 'UAG', 'UAA'],
             'C':['UGC', 'UGU'],
             'Y':['UAC', 'UAU'],
             'V':['GUG', 'GUA', 'GUC', 'GUU'],
             'G':['GGG', 'GGA', 'GGC', 'GGU'],
             'A':['GCG', 'GCA', 'GCC', 'GCU'],
             'E':['GAG', 'GAA'],
             'D':['GAC', 'GAU'],
             'P':['CCG', 'CCA', 'CCC', 'CCU'],
             'Q':['CAG', 'CAA'],
             'H':['CAC', 'CAU']
           }


codons2aa = dict()
for aa,cdns in aa2codons.items():
    for cd in cdns:
        codons2aa[cd] = aa


def decode_rna(s):
    assert len(s) % 3 == 0
    naa = len(s) / 3

    aa_str = ''
    for k in range(naa):
        si = k*3
        ei = si + 3
        r = s[si:ei]
        aa = codons2aa[r]
        aa_str += aa

    return aa_str


def rna_multiplicity(s):

    nstrs = 1
    for aa in s:
        codons = aa2codons[aa]
        nstrs *= len(codons)

    return nstrs


def peptide_mass(s):
    m = 0
    for aa in s:
        m += aa_mass[aa]
    return m


def linear_peptide_spectrum(s):

    mass = list()

    maxk = len(s)
    for k in range(maxk):
        for pos in range(len(s)-k):
            kmer = s[pos:(pos+k+1)]
            # print 'kmer=%s' % kmer
            mass.append(peptide_mass(kmer))

    return sorted(mass)


def cyclopeptide_spectrum(s):

    all_masses = dict()

    maxk = len(s)

    nstr = np.array([c for c in s])

    for shift_amt in range(maxk):

        aa_str = np.roll(nstr, shift_amt)

        for k in range(maxk):
            for pos in range(maxk-k):
                kmer = ''.join(aa_str[pos:(pos+k+1)])
                if kmer not in all_masses:
                    # print 'kmer=%s' % kmer
                    all_masses[kmer] = peptide_mass(kmer)

    return sorted(np.unique(all_masses.values()))


if __name__ == '__main__':

    # print decode_rna('CCACGUACUGAAAUUAAC')
    # print decode_rna('CCCAGUACCGAGAUGAAU')
    # print decode_rna('CCGAGGACCGAAAUCAAC')
    # print decode_rna('CCCAGUACCGAAAUUAAC')

    # print rna_multiplicity('LEADER')

    # print cyclopeptide_spectrum('MAIT')
    # print cyclopeptide_spectrum('TMIA')
    # print cyclopeptide_spectrum('TLAM')
    # print cyclopeptide_spectrum('MLAT')
    # print cyclopeptide_spectrum('TALM')
    # print cyclopeptide_spectrum('TAIM')

    print linear_peptide_spectrum('TVQ')
    print linear_peptide_spectrum('AQV')
    print linear_peptide_spectrum('CTQ')
    print linear_peptide_spectrum('CET')
    print linear_peptide_spectrum('AVQ')
    print linear_peptide_spectrum('CTV')
