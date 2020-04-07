// -*- coding: utf-8 -*-
#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <type_traits>

namespace xn
{

/*!
 * @brief
 *
 * @tparam Graph
 */
template <typename Graph>
class VertexView : public Graph
{
  public:
    /*!
     * @brief Construct a new Vertex View object
     *
     * @param[in,out] G
     */
    explicit VertexView(Graph& G)
        : Graph(G)
    {
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto begin() const
    {
        // auto [v_iter, v_end] = boost::vertices(*this);
        // return v_iter;
        return boost::vertices(*this).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto end() const
    {
        // auto [v_iter, v_end] = boost::vertices(*this);
        // return v_end;
        return boost::vertices(*this).second;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cbegin() const
    {
        // auto [v_iter, v_end] = boost::vertices(*this);
        // return v_iter;
        return boost::vertices(*this).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cend() const
    {
        // auto [v_iter, v_end] = boost::vertices(*this);
        // return v_end;
        return boost::vertices(*this).second;
    }
};

/*!
 * @brief
 *
 * @tparam Graph
 */
template <typename Graph>
class EdgeView
{
  private:
    const Graph& _G;

  public:
    /*!
     * @brief Construct a new Edge View object
     *
     * @param[in] G
     */
    explicit EdgeView(const Graph& G)
        : _G {G}
    {
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto begin() const
    {
        // auto [e_iter, e_end] = boost::edges(_G);
        // return e_iter;
        return boost::edges(_G).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto end() const
    {
        // auto [e_iter, e_end] = boost::edges(_G);
        // return e_end;
        return boost::edges(_G).second;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cbegin() const
    {
        // auto [e_iter, e_end] = boost::edges(_G);
        // return e_iter;
        return boost::edges(_G).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cend() const
    {
        // auto [e_iter, e_end] = boost::edges(_G);
        // return e_end;
        return boost::edges(_G).second;
    }
};

/*!
 * @brief
 *
 * @tparam Vertex
 * @tparam Graph
 */
template <typename Vertex, typename Graph>
class AtlasView
{
  private:
    Vertex _v;
    const Graph& _G;

  public:
    /*!
     * @brief Construct a new Atlas View object
     *
     * @param[in] v
     * @param[in] G
     */
    AtlasView(Vertex v, const Graph& G)
        : _v {v}
        , _G {G}
    {
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto begin() const
    {
        // auto [e_iter, e_end] = boost::out_edges(_v, _G);
        // return e_iter;
        return boost::out_edges(_v, _G).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto end() const
    {
        // auto [e_iter, e_end] = boost::out_edges(_v, _G);
        // return e_end;
        return boost::out_edges(_v, _G).second;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cbegin() const
    {
        // auto [e_iter, e_end] = boost::out_edges(_v, _G);
        // return e_iter;
        return boost::out_edges(_v, _G).first;
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto cend() const
    {
        // auto [e_iter, e_end] = boost::out_edges(_v, _G);
        // return e_end;
        return boost::out_edges(_v, _G).second;
    }
};

/*!
 * @brief
 *
 * @tparam Graph
 */
template <typename Graph>
class grAdaptor : public VertexView<Graph>
{
  public:
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using node_t = Vertex;
    using edge_t = typename boost::graph_traits<Graph>::edge_descriptor;

    // using edge_wt_t = decltype( boost::get(boost::edge_weight,
    // std::declval<Graph>()) );

    /*!
     * @brief Construct a new gr Adaptor object
     *
     */
    grAdaptor() = delete;

    /*!
     * @brief Construct a new gr Adaptor object
     *
     * @param[in,out] G
     */
    explicit grAdaptor(Graph& G)
        : VertexView<Graph>(G)
    {
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto number_of_nodes() const
    {
        return boost::num_vertices(*this);
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto number_of_edges() const
    {
        return boost::num_edges(*this);
    }

    /*!
     * @brief
     *
     * @return EdgeView<Graph>
     */
    EdgeView<Graph> edges() const
    {
        return EdgeView<Graph>(*this);
    }

    /*!
     * @brief
     *
     * @param[in] v
     * @return AtlasView<Vertex, Graph>
     */
    AtlasView<Vertex, Graph> neighbors(Vertex v) const
    {
        return AtlasView<Vertex, Graph>(v, *this);
    }

    /*!
     * @brief
     *
     * @param[in] u
     * @param[in] v
     * @return auto
     */
    auto add_edge(int u, int v)
    {
        return boost::add_edge(u, v, *this);
    }

    /*!
     * @brief
     *
     * @return Vertex
     */
    static Vertex null_vertex()
    {
        return boost::graph_traits<Graph>::null_vertex();
    }

    /*!
     * @brief
     *
     * @tparam Edge
     * @param[in] e
     * @return Vertex
     */
    template <typename Edge>
    Vertex source(const Edge& e) const
    {
        return boost::source(e, *this);
    }

    /*!
     * @brief
     *
     * @tparam Edge
     * @param[in] e
     * @return Vertex
     */
    template <typename Edge>
    Vertex target(const Edge& e) const
    {
        return boost::target(e, *this);
    }

    /*!
     * @brief
     *
     * @tparam Edge
     * @param[in] e
     * @return auto
     */
    template <typename Edge>
    auto end_points(const Edge& e) const
    {
        auto s = boost::source(e, *this);
        auto t = boost::target(e, *this);
        return std::pair {s, t};
    }
};

} // namespace xn
