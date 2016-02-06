
import numpy as np
import matplotlib.pyplot as plt
import pygraphviz
import networkx as nx


class SuffixTree(object):

    def __init__(self, s):
        self.s = s
        self.g = nx.DiGraph()

        self.g.add_node(0, label='Root')

        # build the tree
        for k in range(len(s)):
            current_node = 0
            print ''
            print 'k=%d, s[k:]=%s' % (k, s[k:])
            for c in s[k:]:
                match_node = None
                for n in self.g.successors(current_node):
                    print '%d[%s] is successor of %d[%s]' % (n, self.g.node[n]['label'], current_node, self.g.node[current_node]['label'])
                    if self.g.node[n]['label'] == c:
                        print 'Found match at node %d[%s]' % (n, self.g.node[n]['label'])
                        if match_node is None:
                            match_node = n
                        # break

                if match_node is None:
                    nid = self.addconn(current_node, c)
                    match_node = nid

                current_node = match_node

            self.addconn(current_node, '$%d' % k)

        # merge unbranched jawns
        leaves = [n for n in self.g.nodes() if self.g.node[n]['label'][0] == '$']
        for n in leaves:
            self.merge_up(n)

        # count the number of leaves
        leaves = [n for n in self.g.nodes() if self.g.node[n]['label'].find('$') > -1]
        print '# of leaves: %d' % len(leaves)

    def merge_up(self, n):
        if n == 0:
            return

        assert len(self.g.predecessors(n)) == 1, "Node %d[%s] has %d predecessors" % (n, self.g.node[n]['label'], len(self.g.predecessors(n)))
        parent = self.g.predecessors(n)[0]
        if parent == 0:
            return

        nchildren = len(self.g.successors(parent))
        if nchildren == 1:
            mlabel = self.g.node[parent]['label'] + self.g.node[n]['label']
            self.g.node[parent]['label'] = mlabel
            self.g.remove_node(n)
            self.merge_up(parent)

    def addconn(self, current_node, label):
        nid = self.g.number_of_nodes()
        self.g.add_node(nid, label=label)
        self.g.add_edge(current_node, nid)
        return nid

    def show(self):
        pos = nx.pygraphviz_layout(self.g, prog='dot', root=0)
        nx.draw_networkx(self.g, pos=pos, node_color='w', node_size=500,
                         labels={n:self.g.node[n]['label'] for n in self.g.nodes()})
        plt.show()


if __name__ == '__main__':

    st = SuffixTree('TCTGAGCCCTACTGTCGAGAAATATGTATCTCGCCCCCGCAGCTT')
    st.show()




