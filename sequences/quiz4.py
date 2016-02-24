import operator
from copy import deepcopy

import numpy as np


aa2mass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115,
           'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}

masses = dict()
for aa,m in aa2mass.items():
    if m not in masses:
        masses[m] = list()
    masses[m].append(aa)
mass2aa = {m:'/'.join(aa) for m,aa in masses.items()}
del masses

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


def spectral_convolution(spec):
    diffs = dict()
    for k in range(len(spec)):
        for j in range(len(spec)-k-1):
            d = spec[len(spec)-k-1] - spec[j]
            if d not in diffs:
                diffs[d] = 0
            diffs[d] += 1

    aa_diffs = list()
    num_unknown = 0
    for mass,ndiffs in diffs.items():
        if mass in mass2aa:
            aa = mass2aa[mass]
        else:
            aa = 'Unk%d' % num_unknown
            num_unknown += 1
        aa_diffs.append( (aa, mass, ndiffs) )
    aa_diffs.sort(key=operator.itemgetter(-1), reverse=True)

    return aa_diffs


def peptide_mass(s, aa2mass_dict=aa2mass):
    m = 0
    for aa in s:
        m += aa2mass_dict[aa]
    return m


def cyclopeptide_spectrum(s, aa2mass_dict=aa2mass):

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
                    all_masses[kmer] = peptide_mass(kmer, aa2mass_dict)
    all_masses['None'] = 0

    return sorted(all_masses.values())


def prune_spectra(exp_spec, specs, topn):

    parent_mass = np.sum(exp_spec)

    spec_list = [(aa_str,spec,spec_score(exp_spec, spec)) for aa_str,spec in specs]
    spec_list.sort(key=operator.itemgetter(-1), reverse=True)
    sorted_strs = [x[0] for x in spec_list if np.sum(x) <= parent_mass]
    maxi = min(topn, len(sorted_strs))
    return sorted_strs[:maxi]


def linear_peptide_spectrum(s):

    mass = [0]

    maxk = len(s)
    for k in range(maxk):
        for pos in range(len(s)-k):
            kmer = s[pos:(pos+k+1)]
            # print 'kmer=%s' % kmer
            mass.append(peptide_mass(kmer))

    return sorted(mass)


def leaderboard_spectrum(exp_spec, topn=5, nd_thresh=2):

    allowed_masses = deepcopy(aa2mass)

    aa_diffs = spectral_convolution(exp_spec)
    for aa,m,nd in aa_diffs:
        if aa.startswith('Unk') and nd > nd_thresh:
            allowed_masses[aa] = m

    parent_mass_met = False

    contenders = allowed_masses.keys()

    while not parent_mass_met:

        good_spectra = [(aa_str, cyclopeptide_spectrum(aa_str, allowed_masses)) for aa_str in contenders]
        pruned_contenders = prune_spectra(exp_spec, good_spectra, topn=topn)
        if len(pruned_contenders) == 0:
            parent_mass_met = True
            break

        # create new spectra by adding amino acids
        new_contenders = list()
        for aa_str in pruned_contenders:
            for aa in allowed_masses.keys():
                new_contenders.append(aa_str + aa)

        contenders = new_contenders

def spec_score(spec1, spec2):

    shorter_spec = spec1
    longer_spec = spec2
    if len(spec2) < spec1:
        shorter_spec = spec2
        longer_spec = spec1

    num_matches = 0
    imatch = list()
    for k,m in enumerate(shorter_spec):
        for j,m2 in enumerate(longer_spec):
            if m == m2 and j not in imatch:
                imatch.append(j)
                num_matches += 1
                break

    return num_matches


if __name__ == '__main__':

    # real_seq = 'MAIT'
    # spec = [71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]

    # aa_diffs = spectral_convolution(spec)
    # for aa,m,nd in aa_diffs:
    #     print '%s: mass=%d, # of diffs=%d' % (aa, m, nd)

    """
    pep = 'NQEL'
    spec = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]

    pep_spec = linear_peptide_spectrum(pep)
    print 'pep=%s' % pep
    print 'spec=',spec
    print 'pep_spec=',pep_spec
    print 'Score=%d' % spec_score(spec, pep_spec)
    """

    """
    pep = 'MAMA'
    spec = [0, 71, 178, 202, 202, 202, 333, 333, 333, 404, 507, 507]

    pep_spec = cyclopeptide_spectrum(pep)
    print 'pep=%s' % pep
    print 'spec=',spec
    print 'pep_spec=',pep_spec
    print 'Score=%d' % spec_score(spec, pep_spec)
    """

    """
    pep = 'PEEP'
    spec = [0, 97, 129, 129, 129, 194, 226, 323, 323, 355, 452]

    pep_spec = linear_peptide_spectrum(pep)
    print 'pep=%s' % pep
    print 'spec=',spec
    print 'pep_spec=',pep_spec
    print 'Score=%d' % spec_score(spec, pep_spec)
    """

    spec = [0, 86, 160, 234, 308, 320, 382]
    aa_diffs = spectral_convolution(spec)
    for aa_name,m,ndiffs in aa_diffs:
        print '%s: m=%d, ndiffs=%d' % (aa_name, m, ndiffs)
