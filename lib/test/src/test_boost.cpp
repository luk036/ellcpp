#include <catch.hpp>
#include <py2cpp/nx2bgl.hpp>

#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// using namespace boost; 
  template <class grAdaptor> struct exercise_vertex {
    //...
    typedef typename boost::graph_traits<grAdaptor>
      ::vertex_descriptor Vertex;

    exercise_vertex(grAdaptor& g_) : g(g_) {}
    //...
    grAdaptor& g;

    void operator()(const Vertex& v) const
    {
      // typedef boost::graph_traits<Graph> GraphTraits;
      // typename boost::property_map<Graph, boost::vertex_index_t>::type 
      auto index = boost::get(boost::vertex_index, g);

      std::cout << "out-edges: ";
      // typename GraphTraits::out_edge_iterator out_i, out_end;
      // typename GraphTraits::edge_descriptor e;
      for (auto e : g.neighbors(v)) {
        auto src = g.source(e);
        auto targ = g.target(e);
        std::cout << "(" << index[src] << "," 
                  << index[targ] << ") ";
      }
      std::cout << std::endl;
    }
    //...
  };


TEST_CASE("Test Boost", "[test_boost]") {
    // create a typedef for the Graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;

    // Make convenient labels for the vertices
    enum { A, B, C, D, E, N };
    const int num_vertices = N;
    const char* name = "ABCDE";

    // writing out the edges in the graph
    typedef std::pair<int, int> Edge;
    Edge edge_array[] = 
    { Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
      Edge(C,E), Edge(B,D), Edge(D,E) };
    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

    // declare a graph object
    Graph g(num_vertices);
    xn::grAdaptor<Graph> G(g);

    // add the edges to the graph object
    for (int i = 0; i < num_edges; ++i)
        G.add_edge(edge_array[i].first, edge_array[i].second);

    // typedef graph_traits<Graph>::vertex_descriptor Vertex;

    // get the property map for vertex indices
    // typedef property_map<Graph, vertex_index_t>::type IndexMap;
    auto index = boost::get(boost::vertex_index, G);

    std::cout << "vertices(g) = ";
    // typedef graph_traits<Graph>::vertex_iterator vertex_iter;
    // std::pair<vertex_iter, vertex_iter> vp;
    for (auto v : G) {
      std::cout << index[v] <<  " ";
    }
    std::cout << std::endl;


    std::cout << "edges(g) = ";
    // graph_traits<Graph>::edge_iterator ei, ei_end;
    for (auto e : G.edges())
        std::cout << "(" << index[boost::source(e, G)] 
                  << "," << index[boost::target(e, G)] << ") ";
    std::cout << std::endl;

    std::for_each(boost::vertices(G).first, boost::vertices(G).second,
                  exercise_vertex<xn::grAdaptor<Graph>>(G));
    
}