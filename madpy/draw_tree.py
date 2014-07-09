import networkx as nx
# from matplotlib.pyplot import *
# For use in ipynb use %pylab inline and comment the above line.

def draw_compressed(f, labels=False):
    '''
    Use networkx to draw a directed graph of the madness function f.
    Currently, this does not guarantee the ordering of leaves within a
    refinement level, n. If you really need to know what you're looking
    at, use labels=True.
    '''
    if not f.compressed: f.compress()
    G = nx.DiGraph()
    n = 0
    while f.d[n] != {}:
        for l in f.d[n]:
            G.add_node((n,l))
            if n != 0:
                G.add_edge((n,l),(n-1,l//2))
        n += 1

    G.reverse(copy=False)
    pos = nx.graphviz_layout(G, prog='dot')

    nx.draw(G, pos, node_size=25, font_size=10, with_labels=labels, arrows=False)
