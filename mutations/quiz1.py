
import numpy as np
import matplotlib.pyplot as plt

import networkx as nx


class SuffixTree(object):

    def __init__(self, s):
        self.s = s
        self.g = nx.DiGraph()

        self.g.add_node(0, label='Root')

        for k in range(len(s)):
            current_node = 0
            for c in s[k:]:
                match_node = None
                for n in self.g.successors(current_node):
                    if self.g.node[n]['label'] == c:
                        match_node = n
                        break

                if match_node is None:
                    nid = self.addconn(current_node, c)
                    match_node = nid

                current_node = match_node

            self.addconn(current_node, '$')

    def addconn(self, current_node, label):
        nid = self.g.number_of_nodes()
        self.g.add_node(nid, label=label)
        self.g.add_edge(current_node, nid)
        return nid

    def show(self):
        sl = nx.spectral_layout(self.g)

        nx.draw(sl)

        nx.draw_networkx_labels(self.g, pos=sl)
        plt.show()


if __name__ == '__main__':

    # st = SuffixTree('panamabananas')
    st = SuffixTree('abcc')

    st.show()




