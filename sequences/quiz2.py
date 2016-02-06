import numpy as np
import matplotlib.pyplot as plt

import networkx as nx

def composition(mers):

    g = nx.DiGraph()

    for kmer in mers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        print '%s %s' % (prefix, suffix)

        if prefix not in g:
            g.add_node(prefix)
        if suffix not in g:
            g.add_node(suffix)

        g.add_edge(prefix, suffix, label=kmer)

    return g

def graph_from_adjacency(d):

    g = nx.DiGraph()
    for n in d.keys():
        g.add_node(n)

    for n,conns in d.items():
        for n2 in conns:
            g.add_edge(n, n2)

    return g


def pair_composition(pair_mers):

    g = nx.DiGraph()

    for kmer in pair_mers:

        k1,k2 = kmer.split('|')

        pre1 = k1[:-1]
        pre2 = k2[:-1]
        suf1 = k1[1:]
        suf2 = k2[1:]

        n1 = '%s,%s' % (pre1,pre2)
        n2 = '%s,%s' % (suf1,suf2)

        if n1 not in g:
            g.add_node(n1)
        if n2 not in g:
            g.add_node(n2)

        g.add_edge(n1, n2, label=kmer)

    return g


if __name__ == '__main__':

    """
    kmer_str = 'AAAT AATG ACCC ACGC ATAC ATCA ATGC CAAA CACC CATA CATC CCAG CCCA CGCT CTCA GCAT GCTC TACG TCAC TCAT TGCA'
    g = composition(kmer_str.split(' '))

    pos = nx.pygraphviz_layout(g, prog='dot', root=0)
    nx.draw_networkx(g, pos=pos, node_color='w', node_size=600,
                     labels={n:n for n in g.nodes()})
    nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=nx.get_edge_attributes(g, 'label'))
    plt.show()
    """

    """
    g = graph_from_adjacency({1:[2, 3, 5], 2:[1, 4], 3:[2, 5], 4:[1, 2, 5], 5:[3, 4]})

    pos = nx.pygraphviz_layout(g, prog='dot', root=0)
    nx.draw_networkx(g, pos=pos, node_color='w', node_size=600,
                     labels={n:n for n in g.nodes()})
    plt.show()
    """
    
    pairs = ['ACC|ATA',
            'ACT|ATT',
            'ATA|TGA',
            'ATT|TGA',
            'CAC|GAT',
            'CCG|TAC',
            'CGA|ACT',
            'CTG|AGC',
            'CTG|TTC',
            'GAA|CTT',
            'GAT|CTG',
            'GAT|CTG',
            'TAC|GAT',
            'TCT|AAG',
            'TGA|GCT',
            'TGA|TCT',
            'TTC|GAA']

    g = pair_composition(pairs)
    # pos = nx.pygraphviz_layout(g, prog='dot', root=0)
    # pos = nx.spectral_layout(g)
    pos = nx.spring_layout(g)
    nx.draw_networkx(g, pos=pos, node_color='w', node_size=1700,
                     labels={n:n for n in g.nodes()})
    nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=nx.get_edge_attributes(g, 'label'))
    plt.show()




