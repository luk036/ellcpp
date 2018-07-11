#ifndef _HOME_UBUNTU_GITHUB_PY2CPP_NX2BGL_HPP
#define _HOME_UBUNTU_GITHUB_PY2CPP_NX2BGL_HPP 1

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

namespace xn {

template <typename Graph>
class VertexView : public Graph {
  public:
    explicit VertexView() = default;
    explicit VertexView(int num_vertices) : Graph(num_vertices) {}

    auto begin() { 
        auto [v_iter, v_end] = boost::vertices(*this);
        return v_iter;
    }
    auto end() { 
        auto [v_iter, v_end] = boost::vertices(*this);
        return v_end;
    }
    auto cbegin() const { 
        auto [v_iter, v_end] = boost::vertices(*this);
        return v_iter;
    }
    auto cend() const { 
        auto [v_iter, v_end] = boost::vertices(*this);
        return v_end;
    }
};


template <typename Graph>
class EdgeView {
  private:
    const Graph& _G;

  public:
    explicit EdgeView(const Graph& G) : _G{G} {}
    auto begin() { 
        auto [e_iter, e_end] = boost::edges(_G);
        return e_iter;
    }
    auto end() { 
        auto [e_iter, e_end] = boost::edges(_G);
        return e_end;
    }
    auto cbegin() const { 
        auto [e_iter, e_end] = boost::edges(_G);
        return e_iter;
    }
    auto cend() const { 
        auto [e_iter, e_end] = boost::edges(_G);
        return e_end;
    }
};


template <typename Vertex, typename Graph>
class AtlasView {
  private:
    const Vertex& _v;
    const Graph& _G;

  public:
    explicit AtlasView(const Vertex& v, const Graph& G) : _v{v}, _G{G} {}
    auto begin() { 
        auto [e_iter, e_end] = boost::out_edges(_v, _G);
        return e_iter;
    }
    auto end() { 
        auto [e_iter, e_end] = boost::out_edges(_v, _G);
        return e_end;
    }
    auto cbegin() const { 
        auto [e_iter, e_end] = boost::out_edges(_v, _G);
        return e_iter;
    }
    auto cend() const { 
        auto [e_iter, e_end] = boost::out_edges(_v, _G);
        return e_end;
    }
};

template <typename Graph>
class grAdaptor : public VertexView<Graph> {
  public:
    explicit grAdaptor() = default;
    explicit grAdaptor(int num_vertices) : VertexView<Graph>(num_vertices) {}

    auto add_edge(int u, int v) {
        return boost::add_edge(u, v, *this);
    }

    EdgeView<Graph> edges() const { return EdgeView<Graph>(*this); }

    template <typename Vertex>
    AtlasView<Vertex, Graph> neighbors(const Vertex& v) const { 
        return AtlasView<Vertex, Graph>(v, *this); 
    }

};
 


}
#endif