//    Copyright (C) 2004-2018 by
//    Wai-Shing Luk <luk036@gmail.com>
//
//
//    All rights reserved.
//    BSD license.
//
// Author:  Wai-Shing Luk (luk036@gmail.com),
//          Pieter Swart (swart@lanl.gov),
//          Dan Schult(dschult@colgate.edu);
/** Base class for undirected graphs.

The Graph class allows any hashable object as a node
and can associate key/value attribute pairs with each undirected edge.

Self-loops are allowed but multiple edges are not (see MultiGraph).

For directed graphs see DiGraph && MultiDiGraph.
*/
// from __future__ import division
// import warnings
// from copy import deepcopy
// from collections import Mapping

#include <xnetwork.hpp> // as xn
#include <xnetwork/classes/coreviews.hpp> // import AtlasView, AdjacencyView
#include <xnetwork/classes/reportviews.hpp> // import NodeView, EdgeView, DegreeView
#include <xnetwork/exception.hpp> // import XNetworkError
// #include <xnetwork/convert.hpp> // as convert
// #include <xnetwork/utils.hpp> // import pairwise
#include <py2cpp/py2cpp.hpp>
#include <any>
#include <memory> // for unique_ptr

struct Node : py::dict<const char*, std::any> {
};

struct Edge {
    Node* source;
    Node* target;
};

struct object : py::dict<const char*, std::any> {
};


class Graph: public object {
    /**
    Base class for undirected graphs.

    A Graph stores nodes && edges with optional data, || attributes.

    Graphs hold undirected edges.  Self loops are allowed but multiple
    (parallel) edges are not.

    Nodes can be arbitrary (hashable) Python objects with optional
    key/value attributes. By convention `None` is not used as a node.

    Edges are represented as links between nodes with optional
    key/value attributes.

    Parameters
    ----------
    incoming_graph_data : input graph (optional, default: None);
        Data to initialize graph. If None (default) an empty
        graph is created.  The data can be any format that is supported
        by the to_xnetwork_graph() function, currently including edge list,
        dict of dicts, dict of lists, XNetwork graph, NumPy matrix
        || 2d ndarray, SciPy sparse matrix, || PyGraphviz graph.

    attr : keyword arguments, optional (default= no attributes);
        Attributes to add to graph as key=value pairs.

    See Also
    --------
    DiGraph
    MultiGraph
    MultiDiGraph
    OrderedGraph

    Examples
    --------
    Create an empty graph structure (a "null graph") with no nodes &&
    no edges.

    >>> G = xn::Graph();

    G can be grown : several ways.

    **Nodes:**

    Add one node at a time) {

    >>> G.add_node(1);

    Add the nodes from any container (a list, dict, set or
    even the lines from a file || the nodes from another graph).

    >>> G.add_nodes_from([2, 3]);
    >>> G.add_nodes_from(range(100, 110));
    >>> H = xn::path_graph(10);
    >>> G.add_nodes_from(H);

    In addition to strings && integers any hashable Python object
    (exception (None) can represent a node, e.g. a customized node object,
    || even another Graph.

    >>> G.add_node(H);

    **Edges:**

    G can also be grown by adding edges.

    Add one edge,

    >>> G.add_edge(1, 2);

    a list of edges,

    >>> G.add_edges_from([(1, 2), (1, 3)]);

    || a collection of edges,

    >>> G.add_edges_from(H.edges);

    If some edges connect nodes not yet : the graph, the nodes
    are added automatically.  There are no errors when adding
    nodes || edges that already exist.

    **Attributes:**

    Each graph, node, && edge can hold key/value attribute pairs
    : an associated attribute dictionary (the keys must be hashable).
    By default these are empty, but can be added || changed using
    add_edge, add_node || direct manipulation of the attribute
    dictionaries named graph, node && edge respectively.

    >>> G = xn::Graph(day="Friday");
    >>> G.graph
    {"day": "Friday"}

    Add node attributes using add_node(), add_nodes_from() || G.nodes

    >>> G.add_node(1, time="5pm");
    >>> G.add_nodes_from([3], time="2pm");
    >>> G.nodes[1];
    {"time": "5pm"}
    >>> G.nodes[1]["room"] = 714  // node must exist already to use G.nodes
    >>> del G.nodes[1]["room"];  // remove attribute
    >>> list(G.nodes(data=true));
    [(1, {"time": "5pm"}), (3, {"time": "2pm"})];

    Add edge attributes using add_edge(), add_edges_from(), subscript
    notation, || G.edges.

    >>> G.add_edge(1, 2, weight=4.7 );
    >>> G.add_edges_from([(3, 4), (4, 5)], color="red");
    >>> G.add_edges_from([(1, 2, {"color": "blue"}), (2, 3, {"weight": 8})]);
    >>> G[1][2]["weight"] = 4.7
    >>> G.edges[1, 2]["weight"] = 4

    Warning: we protect the graph data structure by making `G.edges` a
    read-only dict-like structure. However, you can assign to attributes
    : e.g. `G.edges[1, 2]`. Thus, use 2 sets of brackets to add/change
    data attributes: `G.edges[1, 2]["weight"] = 4`
    (For multigraphs: `MG.edges[u, v, key][name] = value`).

    **Shortcuts:**

    Many common graph features allow python syntax to speed reporting.

    >>> 1 : G     // check if (node : graph
    true
    >>> [n for n : G if (n < 3];  // iterate through nodes
    [1, 2];
    >>> len(G);  // number of nodes : graph
    5

    Often the best way to traverse all edges of a graph is via the neighbors.
    The neighbors are reported as an adjacency-dict `G.adj` || as `G.adjacency()`

    >>> for n, nbrsdict : G.adjacency() {
    ...     for (auto nbr, eattr : nbrsdict.items() {
    ...        if ("weight" : eattr) {
    ...            // Do something useful with the edges
    ...            // pass;

    But the edges() method is often more convenient) {

    >>> for u, v, weight : G.edges.data("weight") {
    ...     if (weight is not None) {
    ...         // Do something useful with the edges
    ...         // pass;

    **Reporting:**

    Simple graph information is obtained using object-attributes && methods.
    Reporting typically provides views instead of containers to reduce memory
    usage. The views update as the graph is updated similarly to dict-views.
    The objects `nodes, `edges` && `adj` provide access to data attributes
    via lookup (e.g. `nodes[n], `edges[u, v]`, `adj[u][v]`) && iteration
    (e.g. `nodes.items()`, `nodes.data("color")`,
    `nodes.data("color", default="blue")` && similarly for `edges`);
    Views exist for `nodes`, `edges`, `neighbors()`/`adj` && `degree`.

    For details on these && other miscellaneous methods, see below.

    **Subclasses (Advanced) {**

    The Graph class uses a dict-of-dict-of-dict data structure.
    The outer dict (node_dict) holds adjacency information keyed by node.
    The next dict (adjlist_dict) represents the adjacency information && holds
    edge data keyed by neighbor.  The inner dict (edge_attr_dict) represents
    the edge data && holds edge attribute values keyed by attribute names.

    Each of these three dicts can be replaced : a subclass by a user defined
    dict-like object. In general, the dict-like features should be
    maintained but extra features can be added. To replace one of the
    dicts create a new graph class by changing the class(!) variable
    holding the factory for that dict-like structure. The variable names are
    node_dict_factory, adjlist_inner_dict_factory, adjlist_outer_dict_factory,
    && edge_attr_dict_factory.

    node_dict_factory : function, (default: dict);
        Factory function to be used to create the dict containing node
        attributes, keyed by node id.
        It should require no arguments && return a dict-like object

    adjlist_outer_dict_factory : function, (default: dict);
        Factory function to be used to create the outer-most dict
        : the data structure that holds adjacency info keyed by node.
        It should require no arguments && return a dict-like object.

    adjlist_inner_dict_factory : function, (default: dict);
        Factory function to be used to create the adjacency list
        dict which holds edge data keyed by neighbor.
        It should require no arguments && return a dict-like object

    edge_attr_dict_factory : function, (default: dict);
        Factory function to be used to create the edge attribute
        dict which holds attrbute values keyed by attribute name.
        It should require no arguments && return a dict-like object.

    Examples
    --------

    Create a low memory graph class that effectively disallows edge
    attributes by using a single attribute dict for all edges.
    This reduces the memory used, but you lose edge attributes.

    >>> class ThinGraph(xn::Graph) {
    ...     all_edge_dict = {"weight": 1}
    ...     auto single_edge_dict( ) {
    ...         return this->all_edge_dict
    ...     edge_attr_dict_factory = single_edge_dict
    >>> G = ThinGraph();
    >>> G.add_edge(2, 1);
    >>> G[2][1];
    {"weight": 1}
    >>> G.add_edge(2, 2);
    >>> G[2][1] is G[2][2];
    true

    Please see :mod:`~xnetwork.classes.ordered` for more examples of
    creating graph subclasses by overwriting the base class `dict` with
    a dictionary-like object.
     */
  private:
    using edge_attr_dict_factory = py::dict<const char*, std::any>;
    using adjlist_outer_dict_factory = py::dict<Node*, edge_attr_dict_factory>;
    using adjlist_inner_dict_factory = py::dict<Node*, edge_attr_dict_factory>;
    using node_dict_factory = py::dict<Node*, adjlist_outer_dict_factory>;
    using ndf = node_dict_factory;

  private:
    std::vector<std::unique_ptr<Node> > _Nodes{};

    Graph&  root_graph{*this};
    ndf _node{};  // empty node attribute dict
    edge_attr_dict_factory graph{};   // dictionary for graph attributes
    adjlist_outer_dict_factory _adj{};  // empty adjacency dict

    // auto __getstate__( ) {
    //     attr = this->__dict__.copy();
    //     // remove lazy property attributes
    //     if ("nodes" : attr) {
    //         del attr["nodes"];
    //     }
    //     if ("edges" : attr) {
    //         del attr["edges"];
    //     }
    //     if ("degree" : attr) {
    //         del attr["degree"];
    //     }
    //     return attr;
    // }

    // explicit _Self( incoming_graph_data=None, **attr) {
    explicit Graph( ) {
        /** Initialize a graph with edges, name, || graph attributes.

        Parameters
        ----------
        incoming_graph_data : input graph (optional, default: None);
            Data to initialize graph. If None (default) an empty
            graph is created.  The data can be an edge list, || any
            XNetwork graph object.  If the corresponding optional Python
            packages are installed the data can also be a NumPy matrix
            || 2d ndarray, a SciPy sparse matrix, || a PyGraphviz graph.

        attr : keyword arguments, optional (default= no attributes);
            Attributes to add to graph as key=value pairs.

        See Also
        --------
        convert

        Examples
        --------
        >>> G = xn::Graph();  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G = xn::Graph(name="my graph");
        >>> e = [(1, 2), (2, 3), (3, 4)];  // list of edges
        >>> G = xn::Graph(e);

        Arbitrary graph attribute pairs (key=value) may be assigned

        >>> G = xn::Graph(e, day="Friday");
        >>> G.graph
        {"day": "Friday"}

         */
        // attempt to load graph with data
        // if (incoming_graph_data is not None) {
        //     convert.to_xnetwork_graph(incoming_graph_data, create_using=*this);
        // }
        // load graph attributes (must be after convert);
        // this->graph.update(attr);
    }

    /// @property
    auto get_adj( ) {
        /** Graph adjacency object holding the neighbors of each node.

        This object is a read-only dict-like structure with node keys
        && neighbor-dict values.  The neighbor-dict is keyed by neighbor
        to the edge-data-dict.  So `G.adj[3][2]["color"] = "blue"` sets
        the color of the edge `(3, 2)` to `"blue"`.

        Iterating over G.adj behaves like a dict. Useful idioms include
        `for nbr, datadict : G.adj[n].items() {`.

        The neighbor information is also provided by subscripting the graph.
        So `for nbr, foovalue : G[node].data("foo", default=1) {` works.

        For directed graphs, `G.adj` holds outgoing (successor) info.
         */
        return AdjacencyView(this->_adj);
    }

    /// @property
    auto get_name( ) {
        /** String identifier of the graph.

        This graph attribute appears : the attribute dict G.graph
        keyed by the string `"name"`. as well as an attribute (technically
        a property) `G.name`. This is entirely user controlled.
         */
        if (!this->graph.contains("name"))
            return "";
        return std::any_cast<const char*>(this->graph["name"]);
    }

    // @name.setter
    auto set_name(const char* s) {
        this->graph["name"] = std::any(s);
    }

    friend const char* str(Graph&);

    auto __iter__( ) {
        /** Iterate over the nodes. Use: "for n : G".

        Returns
        -------
        niter : iterator
            An iterator over all nodes : the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> [n for n : G];
        [0, 1, 2, 3];
        >>> list(G);
        [0, 1, 2, 3];
         */
        return iter(this->_node);
    }

    bool contains(const Node *n) {
        /** Return true if (n is a node, false otherwise. Use: "n : G".

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> 1 : G
        true
         */
        return this->_node.contains(n);
    }

    friend std::size_t len(const Graph&);

    auto operator[](const Node* n) {
        /** Return a dict of neighbors of node n.  Use: "G[n]".

        Parameters
        ----------
        n : node
           A node : the graph.

        Returns
        -------
        adj_dict : dictionary
           The adjacency dictionary for nodes connected to n.

        Notes
        -----
        G[n] is the same as G.adj[n] && similar to G.neighbors(n);
        (which is an iterator over G.adj[n]);

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G[0];
        AtlasView({1: {}});
         */
        return this->adj[n];
    }

    // auto add_node(std::any node_for_adding, **attr) {
    auto add_node(std::any node_for_adding,
                  std::initializer_list<std::pair<const char*, std::any>> attr) {
        /** Add a single node `node_for_adding` && update node attributes.

        Parameters
        ----------
        node_for_adding : node
            A node can be any hashable Python object } catch (None.
        attr : keyword arguments, optional
            Set || change node attributes using key=value.

        See Also
        --------
        add_nodes_from

        Examples
        --------
        >>> G = xn::Graph();  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_node(1);
        >>> G.add_node("Hello");
        >>> K3 = xn::Graph([(0, 1), (1, 2), (2, 0)]);
        >>> G.add_node(K3);
        >>> G.number_of_nodes();
        3

        Use keywords set/change node attributes) {

        >>> G.add_node(1, size=10);
        >>> G.add_node(3, weight=0.4, UTM=("13S", 382871, 3972649));

        Notes
        -----
        A hashable object is one that can be used as a key : a Python
        dictionary. This includes strings, numbers, tuples of strings
        && numbers, etc.

        On many platforms hashable items also include mutables such as
        XNetwork Graphs, though one should be careful that the hash
        doesn"t change on mutables.
         */
        if (!this->_node.contains(node_for_adding) ) {
            auto n = _Nodes.emplace_back(std::unique_ptr{Node});
            this->_adj[node_for_adding] = this->adjlist_inner_dict_factory{};
            this->_node[node_for_adding] = attr;
        } else { //update attr even if (node already exists
            this->_node[node_for_adding].update(attr);
        }
    }

    auto add_nodes_from( nodes_for_adding, **attr) {
        /** Add multiple nodes.

        Parameters
        ----------
        nodes_for_adding : iterable container
            A container of nodes (list, dict, set, etc.).
            OR
            A container of (node, attribute dict) tuples.
            Node attributes are updated using the attribute dict.
        attr : keyword arguments, optional (default= no attributes);
            Update attributes for all nodes : nodes.
            Node attributes specified : nodes as a tuple take
            precedence over attributes specified via keyword arguments.

        See Also
        --------
        add_node

        Examples
        --------
        >>> G = xn::Graph();  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_nodes_from("Hello");
        >>> K3 = xn::Graph([(0, 1), (1, 2), (2, 0)]);
        >>> G.add_nodes_from(K3);
        >>> sorted(G.nodes(), key=str);
        [0, 1, 2, "H", "e", "l", "o"];

        Use keywords to update specific node attributes for every node.

        >>> G.add_nodes_from([1, 2], size=10);
        >>> G.add_nodes_from([3, 4], weight=0.4);

        Use (node, attrdict) tuples to update attributes for specific nodes.

        >>> G.add_nodes_from([(1, dict(size=11)), (2, {"color":"blue"})]);
        >>> G.nodes[1]["size"];
        11
        >>> H = xn::Graph();
        >>> H.add_nodes_from(G.nodes(data=true));
        >>> H.nodes[1]["size"];
        11

         */
        for (auto n : nodes_for_adding) {
            // keep all this inside try/} catch (because
            // CPython throws TypeError on n not : this->_node,
            // while (pre-2.7.5 ironpython throws on this->_adj[n];
            try {
                if (n not : this->_node) {
                    this->_adj[n] = this->adjlist_inner_dict_factory();
                    this->_node[n] = attr.copy();
                } else {
                    this->_node[n].update(attr);
                }
            } catch (TypeError) {
                nn, ndict = n;
                if (nn not : this->_node) {
                    this->_adj[nn] = this->adjlist_inner_dict_factory();
                    newdict = attr.copy();
                    newdict.update(ndict);
                    this->_node[nn] = newdict;
                } else {
                    olddict = this->_node[nn];
                    olddict.update(attr);
                    olddict.update(ndict);
                }
            }
        }
    }

    auto remove_node( n) {
        /** Remove node n.

        Removes the node n && all adjacent edges.
        Attempting to remove a non-existent node will throw an exception.

        Parameters
        ----------
        n : node
           A node : the graph

        Raises
        -------
        XNetworkError
           If n is not : the graph.

        See Also
        --------
        remove_nodes_from

        Examples
        --------
        >>> G = xn::path_graph(3);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> list(G.edges);
        [(0, 1), (1, 2)];
        >>> G.remove_node(1);
        >>> list(G.edges);
        [];

         */
        adj = this->_adj;
        try {
            nbrs = list(adj[n]);  // list handles self-loops (allows mutation);
            del this->_node[n];
        } catch (KeyError) { //XNetworkError if (n not : self
            throw XNetworkError("The node %s is not in the graph." % (n,));
        }
        for (auto u : nbrs) {
            del adj[u][n];   // remove all edges n-u : graph
        }
        del adj[n];          // now remove node
    }

    auto remove_nodes_from( nodes) {
        /** Remove multiple nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes (list, dict, set, etc.).  If a node
            : the container is not : the graph it is silently
            ignored.

        See Also
        --------
        remove_node

        Examples
        --------
        >>> G = xn::path_graph(3);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> e = list(G.nodes);
        >>> e
        [0, 1, 2];
        >>> G.remove_nodes_from(e);
        >>> list(G.nodes);
        [];

         */
        adj = this->_adj;
        for (auto n : nodes) {
            try {
                del this->_node[n];
                for (auto u : list(adj[n]) {   // list handles self-loops
                    del adj[u][n];  // (allows mutation of dict : loop);
                }
                del adj[n];
            } catch (KeyError) {
                // pass;
            }
        }
    }

    /// @property
    auto nodes( ) {
        /** A NodeView of the Graph as G.nodes || G.nodes().

        Can be used as `G.nodes` for data lookup && for set-like operations.
        Can also be used as `G.nodes(data="color", default=None)` to return a
        NodeDataView which reports specific node data but no set operations.
        It presents a dict-like interface as well with `G.nodes.items()`
        iterating over `(node, nodedata)` 2-tuples && `G.nodes[3]["foo"]`
        providing the value of the `foo` attribute for node `3`. In addition,
        a view `G.nodes.data("foo")` provides a dict-like interface to the
        `foo` attribute of each node. `G.nodes.data("foo", default=1)`
        provides a default for nodes that do not have attribute `foo`.

        Parameters
        ----------
        data : string || bool, optional (default=false);
            The node attribute returned : 2-tuple (n, ddict[data]).
            If true, return entire node attribute dict as (n, ddict).
            If false, return just the nodes n.

        default : value, optional (default=None);
            Value used for nodes that don"t have the requested attribute.
            Only relevant if (data is not true || false.

        Returns
        -------
        NodeView
            Allows set-like operations over the nodes as well as node
            attribute dict lookup && calling to get a NodeDataView.
            A NodeDataView iterates over `(n, data)` && has no set operations.
            A NodeView iterates over `n` && includes set operations.

            When called, if (data == false, an iterator over nodes.
            Otherwise an iterator of 2-tuples (node, attribute value);
            where the attribute is specified : `data`.
            If data is true then the attribute becomes the
            entire data dictionary.

        Notes
        -----
        If your node data is not needed, it is simpler && equivalent
        to use the expression ``for n : G``, || ``list(G)``.

        Examples
        --------
        There are two simple ways of getting a list of all nodes : the graph) {

        >>> G = xn::path_graph(3);
        >>> list(G.nodes);
        [0, 1, 2];
        >>> list(G);
        [0, 1, 2];

        To get the node data along with the nodes) {

        >>> G.add_node(1, time="5pm");
        >>> G.nodes[0]["foo"] = "bar";
        >>> list(G.nodes(data=true));
        [(0, {"foo": "bar"}), (1, {"time": "5pm"}), (2, {})];
        >>> list(G.nodes.data());
        [(0, {"foo": "bar"}), (1, {"time": "5pm"}), (2, {})];

        >>> list(G.nodes(data="foo"));
        [(0, "bar"), (1, None), (2, None)];
        >>> list(G.nodes.data("foo"));
        [(0, "bar"), (1, None), (2, None)];

        >>> list(G.nodes(data="time"));
        [(0, None), (1, "5pm"), (2, None)];
        >>> list(G.nodes.data("time"));
        [(0, None), (1, "5pm"), (2, None)];

        >>> list(G.nodes(data="time", default="Not Available"));
        [(0, "Not Available"), (1, "5pm"), (2, "Not Available")];
        >>> list(G.nodes.data("time", default="Not Available"));
        [(0, "Not Available"), (1, "5pm"), (2, "Not Available")];

        If some of your nodes have an attribute && the rest are assumed
        to have a default attribute value you can create a dictionary
        from node/attribute pairs using the `default` keyword argument
        to guarantee the value is never None:) {

            >>> G = xn::Graph();
            >>> G.add_node(0);
            >>> G.add_node(1, weight=2);
            >>> G.add_node(2, weight=3);
            >>> dict(G.nodes(data="weight", default=1));
            {0: 1, 1: 2, 2: 3}

         */
        nodes = NodeView( );
        // Lazy View creation: overload the (class) property on the instance
        // Then future G.nodes use the existing View
        // setattr doesn"t work because attribute already exists
        this->__dict__["nodes"] = nodes;
        return nodes;
    }

    // for backwards compatibility with 1.x, will be removed for 3.x
    node = nodes;

    auto add_path( nodes, **attr) {
        const auto msg = "add_path is deprecated. Use xn::add_path instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::add_path( nodes, **attr);
    }

    auto add_cycle( nodes, **attr) {
        const auto msg = "add_cycle is deprecated. Use xn::add_cycle instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::add_cycle( nodes, **attr);
    }

    auto add_star( nodes, **attr) {
        const auto msg = "add_star is deprecated. Use xn::add_star instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::add_star( nodes, **attr);
    }

    auto nodes_with_selfloops( ) {
        const auto msg = "nodes_with_selfloops is deprecated." \
              "Use xn::nodes_with_selfloops instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::nodes_with_selfloops( );
    }

    auto number_of_selfloops( ) {
        const auto msg = "number_of_selfloops is deprecated." \
              "Use xn::number_of_selfloops instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::number_of_selfloops( );
    }

    auto selfloop_edges( data=false, keys=false, default=None) {
        const auto msg = "selfloop_edges is deprecated. Use xn::selfloop_edges instead."
        warnings.warn(msg, DeprecationWarning);
        return xn::selfloop_edges( data=false, keys=false, default=None);
    }
    // Done with backward compatibility methods for 1.x

    auto number_of_nodes( ) {
        /** Return the number of nodes : the graph.

        Returns
        -------
        nnodes : int
            The number of nodes : the graph.

        See Also
        --------
        order, __len__  which are identical

        Examples
        --------
        >>> G = xn::path_graph(3);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> len(G);
        3
         */
        return len(this->_node);
    }

    auto order( ) {
        /** Return the number of nodes : the graph.

        Returns
        -------
        nnodes : int
            The number of nodes : the graph.

        See Also
        --------
        number_of_nodes, __len__  which are identical

         */
        return len(this->_node);
    }

    auto has_node( n) {
        /** Return true if (the graph contains the node n.

        Identical to `n : G`

        Parameters
        ----------
        n : node

        Examples
        --------
        >>> G = xn::path_graph(3);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.has_node(0);
        true

        It is more readable && simpler to use

        >>> 0 : G
        true

         */
        try {
            return n : this->_node;
        } catch (TypeError) {
            return false;
        }
    }

    auto add_edge( u_of_edge, v_of_edge, **attr) {
        /** Add an edge between u && v.

        The nodes u && v will be automatically added if (they are
        not already : the graph.

        Edge attributes can be specified with keywords || by directly
        accessing the edge"s attribute dictionary. See examples below.

        Parameters
        ----------
        u, v : nodes
            Nodes can be, for example, strings || numbers.
            Nodes must be hashable (and not None) Python objects.
        attr : keyword arguments, optional
            Edge data (or labels || objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edges_from : add a collection of edges

        Notes
        -----
        Adding an edge that already exists updates the edge data.

        Many XNetwork algorithms designed for weighted graphs use
        an edge attribute (by default `weight`) to hold a numerical value.

        Examples
        --------
        The following all add the edge e=(1, 2) to graph G) {

        >>> G = xn::Graph()   // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> e = (1, 2);
        >>> G.add_edge(1, 2)           // explicit two-node form
        >>> G.add_edge(*e)             // single edge as tuple of two nodes
        >>> G.add_edges_from([(1, 2)]);  // add edges from iterable container

        Associate data to edges using keywords) {

        >>> G.add_edge(1, 2, weight=3);
        >>> G.add_edge(1, 3, weight=7, capacity=15, length=342.7);

        For non-string attribute keys, use subscript notation.

        >>> G.add_edge(1, 2);
        >>> G[1][2].update({0: 5});
        >>> G.edges[1, 2].update({0: 5});
         */
        auto [u, v] = u_of_edge, v_of_edge
        // add nodes
        if (u not : this->_node) {
            this->_adj[u] = this->adjlist_inner_dict_factory();
            this->_node[u] = {};
        }
        if (v not : this->_node) {
            this->_adj[v] = this->adjlist_inner_dict_factory();
            this->_node[v] = {};
        }
        // add the edge
        datadict = this->_adj[u].get(v, this->edge_attr_dict_factory());
        datadict.update(attr);
        this->_adj[u][v] = datadict;
        this->_adj[v][u] = datadict;
    }

    auto add_edges_from( ebunch_to_add, **attr) {
        /** Add all the edges : ebunch_to_add.

        Parameters
        ----------
        ebunch_to_add : container of edges
            Each edge given : the container will be added to the
            graph. The edges must be given as as 2-tuples (u, v) or
            3-tuples (u, v, d) where d is a dictionary containing edge data.
        attr : keyword arguments, optional
            Edge data (or labels || objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edge : add a single edge
        add_weighted_edges_from : convenient way to add weighted edges

        Notes
        -----
        Adding the same edge twice has no effect but any edge data
        will be updated when each duplicate edge is added.

        Edge attributes specified : an ebunch take precedence over
        attributes specified via keyword arguments.

        Examples
        --------
        >>> G = xn::Graph()   // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edges_from([(0, 1), (1, 2)]) // using a list of edge tuples
        >>> e = zip(range(0, 3), range(1, 4));
        >>> G.add_edges_from(e) // Add the path graph 0-1-2-3

        Associate data to edges

        >>> G.add_edges_from([(1, 2), (2, 3)], weight=3);
        >>> G.add_edges_from([(3, 4), (1, 4)], label="WN2898");
         */
        for (auto e : ebunch_to_add) {
            ne = len(e);
            if (ne == 3) {
                u, v, dd = e;
            } else if (ne == 2) {
                auto [u, v] = e;
                dd = {};  // doesn"t need edge_attr_dict_factory;
            } else {
                throw XNetworkError(
                    "Edge tuple %s must be a 2-tuple || 3-tuple." % (e,));
            }
            if (u not : this->_node) {
                this->_adj[u] = this->adjlist_inner_dict_factory();
                this->_node[u] = {};
            }
            if (v not : this->_node) {
                this->_adj[v] = this->adjlist_inner_dict_factory();
                this->_node[v] = {};
            }
            datadict = this->_adj[u].get(v, this->edge_attr_dict_factory());
            datadict.update(attr);
            datadict.update(dd);
            this->_adj[u][v] = datadict;
            this->_adj[v][u] = datadict;
        }
    }

    auto add_weighted_edges_from( ebunch_to_add, weight="weight", **attr) {
        /** Add weighted edges : `ebunch_to_add` with specified weight attr

        Parameters
        ----------
        ebunch_to_add : container of edges
            Each edge given : the list || container will be added
            to the graph. The edges must be given as 3-tuples (u, v, w);
            where w is a number.
        weight : string, optional (default= "weight");
            The attribute name for the edge weights to be added.
        attr : keyword arguments, optional (default= no attributes);
            Edge attributes to add/update for all edges.

        See Also
        --------
        add_edge : add a single edge
        add_edges_from : add multiple edges

        Notes
        -----
        Adding the same edge twice for Graph/DiGraph simply updates
        the edge data. For MultiGraph/MultiDiGraph, duplicate edges
        are stored.

        Examples
        --------
        >>> G = xn::Graph()   // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_weighted_edges_from([(0, 1, 3.0), (1, 2, 7.5)]);
         */
        this->add_edges_from(((u, v, {weight: d}) for u, v, d : ebunch_to_add),
                            **attr);
    }

    auto remove_edge( u, v) {
        /** Remove the edge between u && v.

        Parameters
        ----------
        u, v : nodes
            Remove the edge between nodes u && v.

        Raises
        ------
        XNetworkError
            If there is not an edge between u && v.

        See Also
        --------
        remove_edges_from : remove a collection of edges

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, etc
        >>> G.remove_edge(0, 1);
        >>> e = (1, 2);
        >>> G.remove_edge(*e) // unpacks e from an edge tuple
        >>> e = (2, 3, {"weight":7}) // an edge with attribute data
        >>> G.remove_edge(*e[:2]) // select first part of edge tuple
         */
        try {
            del this->_adj[u][v];
            if (u != v) { //self-loop needs only one entry removed
                del this->_adj[v][u];
            }
        } catch (KeyError) {
            throw XNetworkError("The edge %s-%s is not : the graph" % (u, v));
        }
    }

    auto remove_edges_from( ebunch) {
        /** Remove all edges specified : ebunch.

        Parameters
        ----------
        ebunch: list || container of edge tuples
            Each edge given : the list || container will be removed
            from the graph. The edges can be) {

                - 2-tuples (u, v) edge between u && v.
                - 3-tuples (u, v, k) where k is ignored.

        See Also
        --------
        remove_edge : remove a single edge

        Notes
        -----
        Will fail silently if (an edge : ebunch is not : the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> ebunch=[(1, 2), (2, 3)];
        >>> G.remove_edges_from(ebunch);
         */
        adj = this->_adj;
        for (auto e : ebunch) {
            auto [u, v] = e[:2];  // ignore edge data if (present
            if (u : adj && v : adj[u]) {
                del adj[u][v];
                if (u != v) { //self loop needs only one entry removed
                    del adj[v][u];
                }
            }
        }
    }

    auto has_edge( u, v) {
        /** Return true if (the edge (u, v) is : the graph.

        This is the same as `v : G[u]` without KeyError exceptions.

        Parameters
        ----------
        u, v : nodes
            Nodes can be, for example, strings || numbers.
            Nodes must be hashable (and not None) Python objects.

        Returns
        -------
        edge_ind : bool
            true if (edge is : the graph, false otherwise.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.has_edge(0, 1);  // using two nodes
        true
        >>> e = (0, 1);
        >>> G.has_edge(*e);  //  e is a 2-tuple (u, v);
        true
        >>> e = (0, 1, {"weight":7});
        >>> G.has_edge(*e[:2]);  // e is a 3-tuple (u, v, data_dictionary);
        true

        The following syntax are equivalent) {

        >>> G.has_edge(0, 1);
        true
        >>> 1 : G[0];  // though this gives KeyError if (0 not : G
        true

         */
        try {
            return v : this->_adj[u];
        } catch (KeyError) {
            return false;
        }
    }

    auto neighbors( n) {
        /** Return an iterator over all neighbors of node n.

        This is identical to `iter(G[n])`

        Parameters
        ----------
        n : node
           A node : the graph

        Returns
        -------
        neighbors : iterator
            An iterator over all neighbors of node n

        Raises
        ------
        XNetworkError
            If the node n is not : the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> [n for n : G.neighbors(0)];
        [1];

        Notes
        -----
        It is usually more convenient (and faster) to access the
        adjacency dictionary as ``G[n]``) {

        >>> G = xn::Graph()   // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge("a", "b", weight=7);
        >>> G["a"];
        AtlasView({"b": {"weight": 7}});
        >>> G = xn::path_graph(4);
        >>> [n for n : G[0]];
        [1];
         */
        try {
            return iter(this->_adj[n]);
        } catch (KeyError) {
            throw XNetworkError("The node %s is not : the graph." % (n,));
        }
    }

    /// @property
    auto edges( ) {
        /** An EdgeView of the Graph as G.edges || G.edges().

        edges( nbunch=None, data=false, default=None);

        The EdgeView provides set-like operations on the edge-tuples
        as well as edge attribute lookup. When called, it also provides
        an EdgeDataView object which allows control of access to edge
        attributes (but does not provide set-like operations).
        Hence, `G.edges[u, v]["color"]` provides the value of the color
        attribute for edge `(u, v)` while
        `for (auto u, v, c] : G.edges.data("color", default="red") {`
        iterates through all the edges yielding the color attribute
        with default `"red"` if (no color attribute exists.

        Parameters
        ----------
        nbunch : single node, container, || all nodes (default= all nodes);
            The view will only report edges incident to these nodes.
        data : string || bool, optional (default=false);
            The edge attribute returned : 3-tuple (u, v, ddict[data]).
            If true, return edge attribute dict : 3-tuple (u, v, ddict).
            If false, return 2-tuple (u, v).
        default : value, optional (default=None);
            Value used for edges that don"t have the requested attribute.
            Only relevant if (data is not true || false.

        Returns
        -------
        edges : EdgeView
            A view of edge attributes, usually it iterates over (u, v);
            || (u, v, d) tuples of edges, but can also be used for
            attribute lookup as `edges[u, v]["foo"]`.

        Notes
        -----
        Nodes : nbunch that are not : the graph will be (quietly) ignored.
        For directed graphs this returns the out-edges.

        Examples
        --------
        >>> G = xn::path_graph(3)   // or MultiGraph, etc
        >>> G.add_edge(2, 3, weight=5);
        >>> [e for e : G.edges];
        [(0, 1), (1, 2), (2, 3)];
        >>> G.edges.data();  // default data is {} (empty dict);
        EdgeDataView([(0, 1, {}), (1, 2, {}), (2, 3, {"weight": 5})]);
        >>> G.edges.data("weight", default=1);
        EdgeDataView([(0, 1, 1), (1, 2, 1), (2, 3, 5)]);
        >>> G.edges([0, 3]);  // only edges incident to these nodes
        EdgeDataView([(0, 1), (3, 2)]);
        >>> G.edges(0);  // only edges incident to a single node (use G.adj[0]?);
        EdgeDataView([(0, 1)]);
         */
        this->__dict__["edges"] = edges = EdgeView( );
        return edges;
    }

    auto get_edge_data( u, v, default=None) {
        /** Return the attribute dictionary associated with edge (u, v).

        This is identical to `G[u][v]` } catch (the default is returned
        instead of an exception is the edge doesn"t exist.

        Parameters
        ----------
        u, v : nodes
        default:  any Python object (default=None);
            Value to return if (the edge (u, v) is not found.

        Returns
        -------
        edge_dict : dictionary
            The edge attribute dictionary.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G[0][1];
        {}

        Warning: Assigning to `G[u][v]` is not permitted.
        But it is safe to assign attributes `G[u][v]["foo"]`

        >>> G[0][1]["weight"] = 7
        >>> G[0][1]["weight"];
        7
        >>> G[1][0]["weight"];
        7

        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.get_edge_data(0, 1);  // default edge data is {}
        {}
        >>> e = (0, 1);
        >>> G.get_edge_data(*e);  // tuple form
        {}
        >>> G.get_edge_data("a", "b", default=0);  // edge not : graph, return 0
        0
         */
        try {
            return this->_adj[u][v];
        } catch (KeyError) {
            return default;
        }
    }

    auto adjacency( ) {
        /** Return an iterator over (node, adjacency dict) tuples for all nodes.

        For directed graphs, only outgoing neighbors/adjacencies are included.

        Returns
        -------
        adj_iter : iterator
           An iterator over (node, adjacency dictionary) for all nodes in
           the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> [(n, nbrdict) for n, nbrdict : G.adjacency()];
        [(0, {1: {}}), (1, {0: {}, 2: {}}), (2, {1: {}, 3: {}}), (3, {2: {}})];

         */
        return iter(this->_adj.items());
    }

    /// @property
    auto degree( ) {
        /** A DegreeView for the Graph as G.degree || G.degree().

        The node degree is the number of edges adjacent to the node.
        The weighted node degree is the sum of the edge weights for
        edges incident to that node.

        This object provides an iterator for (auto node, degree) as well as
        lookup for the degree for a single node.

        Parameters
        ----------
        nbunch : single node, container, || all nodes (default= all nodes);
            The view will only report edges incident to these nodes.

        weight : string || None, optional (default=None);
           The name of an edge attribute that holds the numerical value used
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        If a single node is requested
        deg : int
            Degree of the node

        OR if (multiple nodes are requested
        nd_view : A DegreeView object capable of iterating (node, degree) pairs

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.degree[0];  // node 0 has degree 1
        1
        >>> list(G.degree([0, 1, 2]));
        [(0, 1), (1, 2), (2, 2)];
         */
        this->__dict__["degree"] = degree = DegreeView( );
        return degree;
    }

    auto clear( ) {
        /** Remove all nodes && edges from the graph.

        This also removes the name, && all graph, node, && edge attributes.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.clear();
        >>> list(G.nodes);
        [];
        >>> list(G.edges);
        [];

         */
        this->_adj.clear();
        this->_node.clear();
        this->graph.clear();
    }

    auto is_multigraph( ) {
        /** Return true if (graph is a multigraph, false otherwise. */
        return false;
    }

    auto is_directed( ) {
        /** Return true if (graph is directed, false otherwise. */
        return false;
    }

    auto fresh_copy( ) {
        /** Return a fresh copy graph with the same data structure.

        A fresh copy has no nodes, edges || graph attributes. It is
        the same data structure as the current graph. This method is
        typically used to create an empty version of the graph.

        Notes
        -----
        If you subclass the base class you should overwrite this method
        to return your class of graph.
         */
        return Graph();
    }

    auto copy( as_view=false) {
        /** Return a copy of the graph.

        The copy method by default returns a shallow copy of the graph
        && attributes. That is, if (an attribute is a container, that
        container is shared by the original an the copy.
        Use Python"s `copy.deepcopy` for new containers.

        If `as_view` is true then a view is returned instead of a copy.

        Notes
        -----
        All copies reproduce the graph structure, but data attributes
        may be handled : different ways. There are four types of copies
        of a graph that people might want.

        Deepcopy -- The default behavior is a "deepcopy" where the graph
        structure as well as all data attributes && any objects they might
        contain are copied. The entire graph object is new so that changes
        : the copy do not affect the original object. (see Python"s
        copy.deepcopy);

        Data Reference (Shallow) -- For a shallow copy the graph structure
        is copied but the edge, node && graph attribute dicts are
        references to those : the original graph. This saves
        time && memory but could cause confusion if (you change an attribute
        : one graph && it changes the attribute : the other.
        XNetwork does not provide this level of shallow copy.

        Independent Shallow -- This copy creates new independent attribute
        dicts && then does a shallow copy of the attributes. That is, any
        attributes that are containers are shared between the new graph
        && the original. This is exactly what `dict.copy()` provides.
        You can obtain this style copy using) {

            >>> G = xn::path_graph(5);
            >>> H = G.copy();
            >>> H = G.copy(as_view=false);
            >>> H = xn::Graph(G);
            >>> H = G.fresh_copy().__class__(G);

        Fresh Data -- For fresh data, the graph structure is copied while
        new empty data attribute dicts are created. The resulting graph
        is independent of the original && it has no edge, node || graph
        attributes. Fresh copies are not enabled. Instead use) {

            >>> H = G.fresh_copy();
            >>> H.add_nodes_from(G);
            >>> H.add_edges_from(G.edges);

        View -- Inspired by dict-views, graph-views act like read-only
        versions of the original graph, providing a copy of the original
        structure without requiring any memory for copying the information.

        See the Python copy module for more information on shallow
        && deep copies, https://docs.python.org/2/library/copy.html.

        Parameters
        ----------
        as_view : bool, optional (default=false);
            If true, the returned graph-view provides a read-only view
            of the original graph without actually copying any data.

        Returns
        -------
        G : Graph
            A copy of the graph.

        See Also
        --------
        to_directed: return a directed copy of the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> H = G.copy();

         */
        if (as_view == true) {
            return xn::graphviews.GraphView( );
        }
        G = this->fresh_copy();
        G.graph.update(this->graph);
        G.add_nodes_from((n, d.copy()) for n, d : this->_node.items());
        G.add_edges_from((u, v, datadict.copy());
                         for (auto u, nbrs : this->_adj.items();
                         for (auto v, datadict : nbrs.items());
        return G;
    }

    auto to_directed( as_view=false) {
        /** Return a directed representation of the graph.

        Returns
        -------
        G : DiGraph
            A directed graph with the same name, same nodes, && with
            each edge (u, v, data) replaced by two directed edges
            (u, v, data) && (v, u, data).

        Notes
        -----
        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data && references.

        This is : contrast to the similar D=DiGraph(G) which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        && deep copies, https://docs.python.org/2/library/copy.html.

        Warning: If you have subclassed Graph to use dict-like objects
        : the data structure, those changes do not transfer to the
        DiGraph created by this method.

        Examples
        --------
        >>> G = xn::Graph();  // or MultiGraph, etc
        >>> G.add_edge(0, 1);
        >>> H = G.to_directed();
        >>> list(H.edges);
        [(0, 1), (1, 0)];

        If already directed, return a (deep) copy

        >>> G = xn::DiGraph();  // or MultiDiGraph, etc
        >>> G.add_edge(0, 1);
        >>> H = G.to_directed();
        >>> list(H.edges);
        [(0, 1)];
         */
        if (as_view == true) {
            return xn::graphviews.DiGraphView( );
        }
        // deepcopy when not a view
        #include <xnetwork.hpp> // import DiGraph
        G = DiGraph();
        G.graph.update(deepcopy(this->graph));
        G.add_nodes_from((n, deepcopy(d)) for n, d : this->_node.items());
        G.add_edges_from((u, v, deepcopy(data));
                         for (auto u, nbrs : this->_adj.items();
                         for (auto v, data : nbrs.items());
        return G;
    }

    auto to_undirected( as_view=false) {
        /** Return an undirected copy of the graph.

        Parameters
        ----------
        as_view : bool (optional, default=false);
          If true return a view of the original undirected graph.

        Returns
        -------
        G : Graph/MultiGraph
            A deepcopy of the graph.

        See Also
        --------
        Graph, copy, add_edge, add_edges_from

        Notes
        -----
        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data && references.

        This is : contrast to the similar `G = xn::DiGraph(D)` which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        && deep copies, https://docs.python.org/2/library/copy.html.

        Warning: If you have subclassed DiGraph to use dict-like objects
        : the data structure, those changes do not transfer to the
        Graph created by this method.

        Examples
        --------
        >>> G = xn::path_graph(2)   // or MultiGraph, etc
        >>> H = G.to_directed();
        >>> list(H.edges);
        [(0, 1), (1, 0)];
        >>> G2 = H.to_undirected();
        >>> list(G2.edges);
        [(0, 1)];
         */
        if (as_view == true) {
            return xn::graphviews.GraphView( );
        }
        // deepcopy when not a view
        G = Graph();
        G.graph.update(deepcopy(this->graph));
        G.add_nodes_from((n, deepcopy(d)) for n, d : this->_node.items());
        G.add_edges_from((u, v, deepcopy(d));
                         for (auto u, nbrs : this->_adj.items();
                         for (auto v, d : nbrs.items());
        return G;
    }

    auto subgraph( nodes) {
        /** Return a SubGraph view of the subgraph induced on `nodes`.

        The induced subgraph of the graph contains the nodes : `nodes`
        && the edges between those nodes.

        Parameters
        ----------
        nodes : list, iterable
            A container of nodes which will be iterated through once.

        Returns
        -------
        G : SubGraph View
            A subgraph view of the graph. The graph structure cannot be
            changed but node/edge attributes can && are shared with the
            original graph.

        Notes
        -----
        The graph, edge && node attributes are shared with the original graph.
        Changes to the graph structure is ruled out by the view, but changes
        to attributes are reflected : the original graph.

        To create a subgraph with its own copy of the edge/node attributes use) {
        G.subgraph(nodes).copy();

        For an inplace reduction of a graph to a subgraph you can remove nodes) {
        G.remove_nodes_from([n for n : G if (n not : set(nodes)]);

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> H = G.subgraph([0, 1, 2]);
        >>> list(H.edges);
        [(0, 1), (1, 2)];
         */
        induced_nodes = xn::filters.show_nodes(this->nbunch_iter(nodes));
        SubGraph = xn::graphviews.SubGraph
        // if (already a subgraph, don"t make a chain
        if (hasattr( "_NODE_OK") {
            return SubGraph(this->_graph, induced_nodes, this->_EDGE_OK);
        }
        return SubGraph( induced_nodes);
    }

    auto edge_subgraph( edges) {
        /** Return the subgraph induced by the specified edges.

        The induced subgraph contains each edge : `edges` && each
        node incident to any one of those edges.

        Parameters
        ----------
        edges : iterable
            An iterable of edges : this graph.

        Returns
        -------
        G : Graph
            An edge-induced subgraph of this graph with the same edge
            attributes.

        Notes
        -----
        The graph, edge, && node attributes : the returned subgraph
        view are references to the corresponding attributes : the original
        graph. The view is read-only.

        To create a full graph version of the subgraph with its own copy
        of the edge || node attributes, use:) {

            >>> G.edge_subgraph(edges).copy();  // doctest: +SKIP

        Examples
        --------
        >>> G = xn::path_graph(5);
        >>> H = G.edge_subgraph([(0, 1), (3, 4)]);
        >>> list(H.nodes);
        [0, 1, 3, 4];
        >>> list(H.edges);
        [(0, 1), (3, 4)];

         */
        return xn::edge_subgraph( edges);
    }

    auto size( weight=None) {
        /** Return the number of edges || total of all edge weights.

        Parameters
        ----------
        weight : string || None, optional (default=None);
            The edge attribute that holds the numerical value used
            as a weight. If None, then each edge has weight 1.

        Returns
        -------
        size : numeric
            The number of edges or
            (if (weight keyword is provided) the total weight sum.

            If weight.empty(), returns an int. Otherwise a double
            (or more general numeric if (the weights are more general).

        See Also
        --------
        number_of_edges

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.size();
        3

        >>> G = xn::Graph()   // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge("a", "b", weight=2);
        >>> G.add_edge("b", "c", weight=4);
        >>> G.size();
        2
        >>> G.size(weight="weight");
        6.0
         */
        s = sum(d for v, d : this->degree(weight=weight));
        // If `weight`.empty(), the sum of the degrees is guaranteed to be
        // even, so we can perform integer division && hence return an
        // integer. Otherwise, the sum of the weighted degrees is not
        // guaranteed to be an integer, so we perform "real" division.
        return s; // 2 if (weight.empty() else s / 2
    }

    auto number_of_edges( u=None, v=None) {
        /** Return the number of edges between two nodes.

        Parameters
        ----------
        u, v : nodes, optional (default=all edges);
            If u && v are specified, return the number of edges between
            u && v. Otherwise return the total number of all edges.

        Returns
        -------
        nedges : int
            The number of edges : the graph.  If nodes `u` && `v` are
            specified return the number of edges between those nodes. If
            the graph is directed, this only returns the number of edges
            from `u` to `v`.

        See Also
        --------
        size

        Examples
        --------
        For undirected graphs, this method counts the total number of
        edges : the graph) {

        >>> G = xn::path_graph(4);
        >>> G.number_of_edges();
        3

        If you specify two nodes, this counts the total number of edges
        joining the two nodes) {

        >>> G.number_of_edges(0, 1);
        1

        For directed graphs, this method can count the total number of
        directed edges from `u` to `v`) {

        >>> G = xn::DiGraph();
        >>> G.add_edge(0, 1);
        >>> G.add_edge(1, 0);
        >>> G.number_of_edges(0, 1);
        1

         */
        if (u.empty()) {
            return int(this->size());
        }
        if (v : this->_adj[u]) {
            return 1;
        }
        return 0;
    }

    auto nbunch_iter( nbunch=None) {
        /** Return an iterator over nodes contained : nbunch that are
        also : the graph.

        The nodes : nbunch are checked for membership : the graph
        && if (!are silently ignored.

        Parameters
        ----------
        nbunch : single node, container, || all nodes (default= all nodes);
            The view will only report edges incident to these nodes.

        Returns
        -------
        niter : iterator
            An iterator over nodes : nbunch that are also : the graph.
            If nbunch.empty(), iterate over all nodes : the graph.

        Raises
        ------
        XNetworkError
            If nbunch is not a node || or sequence of nodes.
            If a node : nbunch is not hashable.

        See Also
        --------
        Graph.__iter__

        Notes
        -----
        When nbunch is an iterator, the returned iterator yields values
        directly from nbunch, becoming exhausted when nbunch is exhausted.

        To test whether nbunch is a single node, one can use
        "if (nbunch : self:", even after processing with this routine.

        If nbunch is not a node || a (possibly empty) sequence/iterator
        || None, a :exc:`XNetworkError` is raised.  Also, if (any object in
        nbunch is not hashable, a :exc:`XNetworkError` is raised.
         */
        if (nbunch.empty()) {   // include all nodes via iterator
            bunch = iter(this->_adj);
        } else if (nbunch : *this) { //if (nbunch is a single node
            bunch = iter([nbunch]);
        } else {                // if (nbunch is a sequence of nodes
            auto bunch_iter(nlist, adj) {
                try {
                    for (auto n : nlist) {
                        if (n : adj) {
                            yield n;
                        }
                    }
                } catch (TypeError as e) {
                    message = e.args[0];
                    // capture error for non-sequence/iterator nbunch.
                    if ("iter" : message) {
                        const auto msg = "nbunch is not a node || a sequence of nodes."
                        throw XNetworkError(msg);
                    // capture error for unhashable node.
                    } else if ("hashable" : message) {
                        const auto msg = "Node {} : sequence nbunch is not a valid node."
                        throw XNetworkError(msg.format(n));
                    } else {
                        throw;
                    }
                }
            }
            bunch = bunch_iter(nbunch, this->_adj);
        }
        return bunch;
    }
};

// Non-member functions   
    std::string str(const Graph& G) {
        /** Return the graph name.

        Returns
        -------
        name : string
            The name of the graph.

        Examples
        --------
        >>> G = xn::Graph(name="foo");
        >>> str(G);
        "foo";
         */
        return G.get_name();
    }

    std::size_t len(const Graph& G) {
        /** Return the number of nodes. Use: "len(G)".

        Returns
        -------
        nnodes : int
            The number of nodes : the graph.

        Examples
        --------
        >>> G = xn::path_graph(4);  // or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> len(G);
        4

         */
        return G._node.size();
    }
