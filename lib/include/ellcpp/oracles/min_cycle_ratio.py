# -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import networkx as nx
from parametric import max_parametric


def set_default(G, weight, value):
    for (u, v) in G.edges:
        if G[u][v].get(weight, None) is None:
            G[u][v][weight] = value


def calc_weight(G, r, e):
    u, v = e
    return G[u][v]['cost'] - r * G[u][v]['time']


def calc_ratio(G, C):
    """Calculate the ratio of the cycle

    Arguments:
        G {Networkx Graph} -- [description]
        C {list} -- [description]

    Returns:
        float -- cycle ratio
    """
    total_cost = sum(G[u][v]['cost'] for (u, v) in C)
    total_time = sum(G[u][v]['time'] for (u, v) in C)
    return total_cost/total_time


def min_cycle_ratio(G):
    mu = 'cost'
    sigma = 'time'
    set_default(G, mu, 1)
    set_default(G, sigma, 1)
    max_cost = max(cost for _, _, cost in G.edges.data(mu))
    min_time = min(time for _, _, time in G.edges.data(sigma))
    r0 = max_cost * G.number_of_edges() / min_time
    return max_parametric(G, r0, calc_weight, calc_ratio)


# if __name__ == "__main__":
#     import networkx as nx
#     from neg_cycle import *
#     from networkx.utils import generate_unique_node

#     G = create_test_case1()
#     G[1][2]['cost'] = 5
#     r, c, dist = min_cycle_ratio(G)
#     assert c != None
#     print(r)
#     print(c)
#     print(dist.items())

#     G = nx.cycle_graph(5, create_using=nx.DiGraph())
#     G[1][2]['cost'] = -6.
#     newnode = generate_unique_node()
#     G.add_edges_from([(newnode, n) for n in G])
#     r, c, dist = min_cycle_ratio(G)
#     assert c != None
#     print(r)
#     print(c)
#     print(dist.items())
