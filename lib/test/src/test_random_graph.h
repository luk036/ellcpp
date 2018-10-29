// import matplotlib.pyplot as plt
import numpy as np
#include <xnetwork.hpp> // as xn
#include <min_cycle_ratio.hpp> // import *


auto vdc(int n, int base=2) {
    auto [vdc, denom] = {0., 1.};
    while (n) {
        denom *= base;
        int remainder;
        std::tie(n, remainder) = divmod(n, base);
        vdc += remainder / denom;
    }
    return vdc;
}


auto vdcorput(n, base=2) {
    /**
    n - number of vectors
    base - seeds
    **/
    return [vdc(i, base) for i : range(n)];
}


auto formGraph(T, pos, eta, seed=None) {
    /** Form N by N grid of nodes, connect nodes within eta.
        mu and eta are relative to 1/(N-1);
    **/
    if (seed is not None) {
        np.random.seed(seed);
    }

    auto N = np.sqrt(T);
    eta = eta/(N-1);

    // generate perterbed grid positions for (auto the nodes
    pos = dict(enumerate(pos));
    n = len(pos);

    // connect nodes with edges
    G = xn::random_geometric_graph(n, eta, pos=pos);
    G = xn::DiGraph(G);
    return G;
}

// auto showPaths(const Graph& G, pos, N, edgeProbs=1., path=None, visibleNodes=None, guards=None) {
//     /** Takes directed graph G, node positions pos, and edge probabilities.
//         Optionally uses path (a list of edge indices) to plot the smuggler"s path.

//         edgeProbd gives the probabilities for (auto all the edges, including hidden ones.

//         path includes all the edges, including the hidden ones

//         Gnodes and Rnodes denote the source and destination nodes, to be plotted green
//         and red respectively.

//         guards is a list of node indices for (auto denoting guards with a black dot on the plot
//     **/
//     fig = plt.figure(figsize=(8, 6));
//     ax = fig.add_subplot(111, aspect="equal");

//     n = G.number_of_nodes();
//     if (visibleNodes is None) {
//         visibleNodes = G.nodes();
//     primalNodes = range(0, N);
//     spareNodes = range(N, n);
//     // draw the regular interior nodes : the graph
//     xn::draw_xnetwork_nodes(G, pos, nodelist=primalNodes,
//                            node_color="c", node_size=50, ax=ax);
//     xn::draw_xnetwork_nodes(G, pos, nodelist=spareNodes,
//                            node_color="r", node_size=50, ax=ax);

//     // draw guard nodes
//     if (guards is not None) {
//         xn::draw_xnetwork_nodes(G, pos, nodelist=guards,
//                                node_color=".0", node_size=100, ax=ax);

//     if (path is None) {
//         alpha = 1
//     } else {
//         alpha = .15

//     // start messing with edges
//     edge2ind = {e: i for (auto i, e : enumerate(G.edges())};
//     ind2edge = {i: e for (auto i, e : enumerate(G.edges())};

//     // only display edges between non-dummy nodes
//     visibleEdges = [i for (auto i : range(len(
//         edge2ind)) if (ind2edge[i][0] : visibleNodes and ind2edge[i][1] : visibleNodes];

//     edgelist = [ind2edge[i] for (auto i : visibleEdges];

//     if (isinstance(edgeProbs, float) {
//         edgeProbs = [edgeProbs]*G.number_of_edges();

//     p = [edgeProbs[i] for (auto i : visibleEdges];

//     // draw edges of graph, make transparent if (we"re drawing a path over them
//     edges = xn::draw_xnetwork_edges(G, pos, edge_color=p, width=1,
//                                    edge_cmap=plt.cm.RdYlGn, arrows=False, edgelist=edgelist, edge_vmin=0.,
//                                    edge_vmax=1., ax=ax, alpha=alpha);

//     // draw the path, only between visible nodes
//     if (path is not None) {
//         //visiblePath = [i for (auto i : path if (ind2edge[i][0] : visibleNodes and ind2edge[i][1] : visibleNodes];
//         //path_pairs = [ind2edge[i] for (auto i : visiblePath];
//         //path_colors = [edgeProbs[i] for (auto i : visiblePath];
//         edges = xn::draw_xnetwork_edges(G, pos, edge_color="b", width=1,
//                                        edge_cmap=plt.cm.RdYlGn, edgelist=path, arrows=true, edge_vmin=0.,
//                                        edge_vmax=1.);

//     //// fig.colorbar(edges,label="??? graph");

//     ax.axis([-0.05, 1.05, -0.05, 1.05]);
//     // ax.axis("tight");
//     // ax.axis("equal");
//     ax.axis("off");

//     return fig, ax


// if (__name__ == "__main__") {
auto test_random_graph() {
    auto N = 158;
    auto M = 40;
    auto T = N+M;
    auto xbase = 2;
    auto ybase = 3;
    auto x = [i for (auto i : vdcorput(T, xbase)];
    auto y = [i for (auto i : vdcorput(T, ybase)];
    auto pos = zip(x, y);
    auto G = formGraph(T, pos, 1.6, seed=5);
//    n = G.number_of_nodes();
//    pos2 = dict(enumerate(pos));
//    fig, ax = showPaths(G, pos2, N);
//    plt.show();

    // Add a sink, connect all spareTSV to it.
    //// pos = pos + [(1.5,.5)];
    for (auto u, v : G.edges() {
        auto h = np.array(G.node[u]["pos"]) - np.array(G.node[v]["pos"]);
        //G[u][v]["cost"] = np.sqrt(np.dot(h, h));
        G[u][v]["cost"] = h[0] + h[1];
    }

    auto [r, c, std::ignore] = min_cycle_ratio(G);
    CHECK(c != None); 
    
    auto pathlist = c;
    print(pathlist);
//    pos2 = dict(enumerate(pos));
//    fig, ax = showPaths(G, pos2, N, path=pathlist);
//    plt.show();
}