//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once

#include <optional>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <edge_list.hpp>
#include <util/intersection_size.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include "algorithms/slinegraph_map.hpp"
#include "algorithms/slinegraph_efficient.hpp"
#include "algorithms/toplexes.hpp"
#include <algorithms/connected_components.hpp>
#include <algorithms/delta_stepping.hpp>
#include <algorithms/betweenness_centrality.hpp>
using namespace nw::graph;
namespace py = pybind11;

namespace nw {
namespace hypergraph {

//forward-declaration
template<class Index_t, typename... Attributes> class Slinegraph;

template<class Index_t, typename... Attributes>
class NWHypergraph {
private:
    //hypergraph in bi-adjacency
    nw::graph::adjacency<0, Attributes...> edges_;
    Index_t max_edge_;
    nw::graph::adjacency<1, Attributes...> nodes_;
    Index_t max_node_;

    //for each hyperedge, map stores its hyperedge neighbor list
    //which the key is the adjacent hyperedge id
    //the value is the number of overlap hypernodes
    std::vector<std::map<size_t, size_t>> edge_neighbor_count_;
    //for each hypernode, map stores its hypernode neighbor list
    //which the key is the adjacent hypernode
    //the value is the number of overlap hyperedges
    std::vector<std::map<size_t, size_t>> node_neighbor_count_;   
public:
    py::array_t<Index_t, py::array::c_style> row_;
    py::array_t<Index_t, py::array::c_style> col_;
    py::array_t<Attributes..., py::array::c_style> data_;
    friend class Slinegraph<Index_t, Attributes...>;
public:
    //constructor
    NWHypergraph() {}
    NWHypergraph(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) : row_(x), col_(y) {
        auto&& el = populate_weighted_edge_list(x, y);
        
        edges_ = nw::graph::adjacency<0, Attributes...>(el);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(el);
        max_node_ = nodes_.size();
    }
    NWHypergraph(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Attributes..., py::array::c_style | py::array::forcecast> &data,
    bool collapse,
    bool remove_duplicates = false) : row_(x), col_(y), data_(data) {
        if (collapse) {
            auto&& [el, dict] = collapse_array_x();
            edges_ = nw::graph::adjacency<0, Attributes...>(el);
            max_edge_ = edges_.size();
            nodes_ = nw::graph::adjacency<1, Attributes...>(el);
            max_node_ = nodes_.size();
        }
        else {
            auto&& el = populate_weighted_edge_list(x, y, data);

            if (remove_duplicates) {
                //Remove duplicate edges
                std::cout << "before collapse: " << el.size() << " edges" << std::endl;
                el.template lexical_sort_by<0>();
                el.uniq();
                std::cout << "after collapse: " << el.size() << " edges" << std::endl;
            }
            edges_ = nw::graph::adjacency<0, Attributes...>(el);
            max_edge_ = edges_.size();
            nodes_ = nw::graph::adjacency<1, Attributes...>(el);
            max_node_ = nodes_.size();
        }
    }
    /*
    NWHypergraph(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Attributes..., py::array::c_style | py::array::forcecast> &data,
    bool collapse,
    bool remove_duplicates = false) : row_(x), col_(y), data_(data) {
        auto&& el = populate_weighted_edge_list(x, y, data);
        if(remove_duplicates) {
            //Remove duplicate edges
            std::cout << "before collapse: " << el.size() << " edges" << std::endl;
            el.template lexical_sort_by<0>();
            el.uniq();
            std::cout << "after collapse: " << el.size() << " edges" << std::endl;
        }
        edges_ = nw::graph::adjacency<0, Attributes...>(el);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(el);
        max_node_ = nodes_.size();
        std::vector<std::set<Index_t>> neighborhoods(max_edge_);
        for (Index_t i = 0; i < max_edge_; ++i) {
            std::set<Index_t> set_i;
            for (auto& i_ngh : edges_[i]) {
                set_i.add(std::get<0>(*i_ngh));
            }
            neighborhoods[i] = set_i;
        }
        std::unordered_map<std::set<Index_t>, std::set<Index_t>> equal_class;
        for (Index_t i = 0; i < max_edge_; ++i) {
            for (Index_t j = i + 1; j < max_edge_; ++j) {   
                if (neighborhoods[i].size() == neighborhoods[j].size()) {
                    if (neighborhoods[i] == neighborhoods[j]) {
                        //then edge i and j are equal class
                        if (equal_class.find(neighborhoods[i]) != equal_class.end()) {
                            equal_class[neighborhoods[i]].add(i);
                            equal_class[neighborhoods[i]].add(j);
                        }
                        else {
                            std::set<Index_t> s;
                            s.add(i);
                            s.add(j);
                            equal_class[neighborhoods[i]] = s;
                        }
                    }
                }          
            }
        }
    }
*/
    decltype(auto) collapse_edges(bool return_equivalence_class = false) {
        auto &&[el, equivalence_class_dict] = collapse_array_x();
        edges_ = nw::graph::adjacency<0, Attributes...>(el);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(el);
        max_node_ = nodes_.size();

        if (return_equivalence_class)
            return std::tuple{*this, equivalence_class_dict};
        else
            return std::tuple{*this, equivalence_class_dict};
    }

    decltype(auto) collapse_nodes(bool return_equivalence_class = false) {
        bool transpose = true;
        auto &&[el, equivalence_class_dict] = collapse_array_x(transpose);
        edges_ = nw::graph::adjacency<0, Attributes...>(el);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(el);
        max_node_ = nodes_.size();

        if (return_equivalence_class)
            return std::tuple{*this, equivalence_class_dict};
        else
            return std::tuple{*this, equivalence_class_dict};   
    }
    
    decltype(auto) collapse_nodes_and_edges(bool return_equivalence_class = false) {
        auto &&[el, equivalence_class_dict] = collapse_array_x();
        edges_ = nw::graph::adjacency<0, Attributes...>(el);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(el);
        max_node_ = nodes_.size();

        if (return_equivalence_class)
            return std::tuple{*this, equivalence_class_dict};
        else
            return std::tuple{*this, equivalence_class_dict};   
    }

    std::vector<std::map<size_t, size_t>> get_edge_neighbor_counts() const { return edge_neighbor_count_; }
    std::vector<std::map<size_t, size_t>> get_node_neighbor_counts() const { return node_neighbor_count_; }

    Slinegraph<Index_t, Attributes...> s_linegraph(int s = 1, bool edges = true) {
        Slinegraph<Index_t, Attributes...> slineg(*this, s, edges);
        return slineg;
    }
    /*
    * return a list of slinegraph based on the s value list
    * */
    py::list s_linegraphs(py::list l, bool edges = true) {
        auto tmp = std::min_element(l.begin(), l.end());
        int min_s = (*tmp).cast<int>();
        if(edges) {
            if (edge_neighbor_count_.empty()) 
                edge_neighbor_count_ = to_two_graph_count_neighbors_cyclic(edges_, nodes_, min_s, 32);
        }
        else {
            if (node_neighbor_count_.empty()) 
                node_neighbor_count_ = to_two_graph_count_neighbors_cyclic(nodes_, edges_, min_s, 32);
        }
        py::list res;
        for (auto obj : l) {
            int s = obj.cast<int>();
            Slinegraph<Index_t, Attributes...> slineg(*this, s, edges);
            res.append(slineg);
        }
        return res;
        
        //return slineg;
    }
    /*
    * Find the connected components for a slinegraph
    */
    py::list s_connected_components(Slinegraph<Index_t, Attributes...>& linegraph) {
        return linegraph.s_connected_components();
    }
    /*
    * Find the distance from the src to dest in its slinegraph
    */
    Index_t distance(Slinegraph<Index_t, Attributes...>& linegraph, Index_t src, Index_t dest) {
        return linegraph.s_distance(src, dest);
    }
    /*
    * Find neighbors of a vertex in its slinegraph
    */
    py::list neighbors(Slinegraph<Index_t, Attributes...>& linegraph, Index_t v) {
        return linegraph.s_neighbors(v);
    } 
    py::ssize_t s_degree(Slinegraph<Index_t, Attributes...>& linegraph, Index_t v) {
        return linegraph.s_degree(v);
    }
    /*
    * Get the edge size distributation of the hypergraph
    */
    py::list edge_size_dist() const {
        auto degs = edges_.degrees();
        py::list l = py::list(max_edge_);
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, max_edge_), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                l[i] = degs[i];
            }
        }, tbb::auto_partitioner());
        return l;
    }
    /*
    * Get the node size distributation of the hypergraph
    */
    py::list node_size_dist() const {
        auto degs = nodes_.degrees();
        py::list l = py::list(max_node_);
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, max_node_), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                l[i] = degs[i];
            }
        }, tbb::auto_partitioner());
        return l;
    }
    /*
    * Get the incident nodes of an edge in the hypergraph
    */
    py::list edge_incidence(Index_t edge) {
        py::list l;
        /*
        std::for_each(edges_[e].begin(), edges_[e].end(), [&](auto &&x) {
            auto hyperN = std::get<0>(x);
            l.append(hyperN);
        });
        */
        for (auto &&[hyperN, w] : edges_[edge]) {
            l.append(hyperN);
        }
        return l;
    }
    /*
    * Get the incident edges of a node in the hypergraph
    */
    py::list node_incidence(Index_t node) {
        py::list l;
        std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
            auto hyperE = std::get<0>(x);
            l.append(hyperE);
        });
        return l;
    }
    /*
    * Get the degree of a node in the hypergraph.
    * If min_size and/or max_size is specified, then either/both are used to filter.
    * */
    Index_t degree(Index_t node, std::size_t min_size = 1, std::optional<std::size_t> max_size = {}) {
        Index_t deg = -1;
        // validate the input
        if (node >= max_node_)
            return deg;
        if (!max_size.has_value()) {
            //if no max_size is specified
            if (1 == min_size)
                return nodes_[node].size();
            deg = 0;
            std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_size
                if (min_size <= edges_[hyperE].size())
                    ++deg;
            });
        }
        else {
            //if max_size is specified
            //and max_size is also 1
            if (1 == min_size && 1 == max_size)
                return nodes_[node].size();
            deg = 0;
            std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_size
                //and is no greater than max_size
                if (min_size <= edges_[hyperE].size() && max_size >= edges_[hyperE].size())
                    ++deg;
            });
        }
        return deg;
    }
    py::ssize_t number_of_nodes() const { return max_node_; }
    py::ssize_t order() const { return max_node_; }
    py::ssize_t size(Index_t edge, std::size_t min_degree = 1, std::optional<std::size_t> max_degree = {}) {
        Index_t size = -1;
        // validate the input
        if (edge >= max_edge_)
            return size;
        if (!max_degree.has_value()) {
            //if no max_size is specified
            if (1 == min_degree)
                return edges_[edge].size();
            size = 0;
            std::for_each(edges_[edge].begin(), edges_[edge].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_degree
                if (min_degree <= edges_[hyperE].size())
                    ++size;
            });
        }
        else {
            //if max_degree is specified
            //and max_degree is also 1
            if (1 == min_degree && 1 == max_degree)
                return edges_[edge].size();
            size = 0;
            std::for_each(edges_[edge].begin(), edges_[edge].end(), [&](auto &&x) {
                auto hyperN = std::get<0>(x);
                //filter the edges whose size is no less than min_degree
                //and is no greater than max_degree
                if (min_degree <= nodes_[hyperN].size() && max_degree >= nodes_[hyperN].size())
                    ++size;
            });
        }
        return size;
    }
    py::ssize_t dim(Index_t edge) {
        if (edge >= max_edge_)
            return -1;
        else
            return edges_[edge].size() - 1;
    }
    py::ssize_t number_of_edges() const { return max_edge_; }
    // a singleton is an edge of size 1
    py::list singletons() {
        py::list l;
        for(Index_t e = 0; e < max_edge_; ++e) {
            if (1 == edges_[e].size()) {
                std::for_each(edges_[e].begin(), edges_[e].end(), [&](auto &&x) {
                    auto hyperN = std::get<0>(x);
                    if (1 == nodes_[hyperN].size())
                        l.append(e);
                });
            }
        }
        return l;
    }
    /*
    * Verified. using intersection_size instead set_intersection
    * */
    py::list toplexes_v2() {
        std::vector<Index_t> tops;
        //std::map<Index_t, bool> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
            bool flag = true;
            std::vector<Index_t> old_tops(tops);
            for (size_t i = 0, end = old_tops.size(); i < end; ++i) {
                Index_t top = old_tops[i];
                //TODO is this necessary?
                if (e == top)
                    continue;
                auto lhs = edges_[e].size();
                auto rhs = edges_[top].size();
                auto s = nw::graph::intersection_size(edges_[e], edges_[top]);
                if (s == rhs)
                    tops.erase(tops.begin() + i);
                else if (s == lhs) {
                    flag = false;
                    break;
                }
            }
            if (flag)
                tops.push_back(e);
        }

        auto n = tops.size();
        py::list l = py::list(n);
        tbb::parallel_for(
            tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t> &r) {
                for (auto i = r.begin(), e = r.end(); i != e; ++i)
                    l[i] = tops[i];
            },
            tbb::auto_partitioner());
        return l;
    }
    /*
    * Verified. using unordered_set
    * */
    py::list toplexes_v1() {
        //check rhs is a subset of lhs
        //create a freq table for all the elements of lhs
        //traverse rhs and search for each element of rhs in the freq table
        //if element is found , then decrease the frequency, otherwise, return false
        //if all elements are found, return true
        //O(m+n), where m is the size of lhs, n is the size of rhs
        auto issubset = []<class A>(A&& lhs, A&& rhs) {
            std::map<Index_t, size_t> frequency;
            std::for_each (lhs.begin(), lhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                ++frequency[v];
            });
            bool res = true;
            std::for_each (rhs.begin(), rhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                if (0 < frequency[v])
                    --frequency[v];
                else {
                    res = false;
                    return;
                }
            }); 
            return res;
        };

        //create an empty toplex set and an empty old toplex set
        std::unordered_set<Index_t> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
            //for each edge, assume it is a toplex by default
            bool flag = true;
            //make a copy of the toplex set as old_tops
            auto old_tops(tops);
            //tops.clear();
            //clear the toplex set
            //TODO Could be parallized:
            //1) if the flag can be set atomically
            //2) once the flag is set, then every thread will exist
            //3) tops.erase has to be within a critical section
            for (auto& top : old_tops) {
                if (e == top)
                    continue; 
               if (issubset(edges_[top], edges_[e])) {
                    //if e is a subset of top, then e is not a toplex
                    flag = false;
                    break;
                } else if (issubset(edges_[e], edges_[top]))
                    tops.erase(top);
            }//for old_tops
            if (flag)
                tops.insert(e);
        }//for each e

        py::list l;
        for (auto& e : tops) {
            l.append(e);
        }
        return l;
    }
    py::list toplexes() {
        auto tops = toplexes_serial_v0(edges_);

        auto n = tops.size();
        py::list l = py::list(n);
        tbb::parallel_for(
            tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t> &r) {
                for (auto i = r.begin(), e = r.end(); i != e; ++i)
                    l[i] = tops[i];
            }, tbb::auto_partitioner());
        return l;
    }
protected:
    auto populate_edge_linegraph (size_t s) {
        //if neighbors have not been counted
        if (edge_neighbor_count_.empty()) {
            edge_neighbor_count_ = to_two_graph_count_neighbors_cyclic(edges_, nodes_);
        }
        return populate_linegraph_from_neighbor_map<nw::graph::undirected, decltype(edges_), Attributes...>(edges_, edge_neighbor_count_, s, false);
    }
    auto populate_node_linegraph (size_t s) {
        //if neighbors have not been counted
        if (0 == node_neighbor_count_.size()) {
            node_neighbor_count_ = to_two_graph_count_neighbors_cyclic(nodes_, edges_);
        }
        return populate_linegraph_from_neighbor_map<nw::graph::undirected, decltype(nodes_), Attributes...>(nodes_, node_neighbor_count_, s, false);
    }
private:
    /*
    * Populate unweighted edge list with numpy arrays.
    * */
    auto populate_unweighted_edge_list(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) {
        nw::graph::edge_list<nw::graph::directed> g(0);
        g.open_for_push_back();
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        for (size_t i = 0; i < n_x; ++i) 
            g.push_back({rx(i), ry(i)});

        g.close_for_push_back(false);
        return g;
    }
    /*
    * Populate weighted edge list with numpy arrays and fake weight 0.
    * */
    auto populate_weighted_edge_list(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) {
        nw::graph::edge_list<nw::graph::directed, Attributes...> g(0);
        g.open_for_push_back();
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        for (size_t i = 0; i < n_x; ++i) 
            g.push_back({rx(i), ry(i), 0});

        g.close_for_push_back(false);
        return g;
    }

    /*
    * Populate weighted edge list with numpy arrays.
    * */
    auto populate_weighted_edge_list(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Attributes..., py::array::c_style | py::array::forcecast> &data) {
        nw::graph::edge_list<nw::graph::directed, Attributes...> g(0);
        g.open_for_push_back();
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        auto rdata = data.template mutable_unchecked<1>();
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        size_t n_data = data.shape(0);
        if (0 == n_data) {
            //when there is no weight passed in, but request a weighted hypergraph
            //we fake a weight
            for (size_t i = 0; i < n_x; ++i) 
                g.push_back({rx(i), ry(i), 0});
        }
        else {
            for (size_t i = 0; i < n_x; ++i) 
                g.push_back({rx(i), ry(i), rdata(i)});
        }
        g.close_for_push_back(false);
        return g;
    }


    auto collapse_array_x(bool transpose = false) {
        std::map<std::set<Index_t>, std::set<Index_t>> equal_class;
        if (transpose) {
            //dict = collapse_identical_elements(row_, col_)
            auto rx = col_.template mutable_unchecked<1>();
            auto ry = row_.template mutable_unchecked<1>();
            size_t n_x = col_.shape(0);
            size_t n_y = row_.shape(0);

        //sort the raw array into a dict, where key is the element of array x
        // and value is the element of array y
        std::map<Index_t, std::set<Index_t>> dict;
        for (size_t i = 0; i < n_x; ++i) {
            auto key = rx(i);
            if (dict.find(key) == dict.end()) {
                //if no such key exists, we create a set and insert as the value
                std::set<Index_t> s;
                s.insert(ry(i));
                dict[key] = s;
            }
            else
                dict[key].insert(ry(i));
        }
        //based on [key, value] pairs of dict, combine its keys if they have the same value.
        //The combined keys will be stored in a new set as the value of equal_class.
        //Whereas the value of dict becomes the key of equal_class.
        
        for (auto&& [k, v] : dict) {
            auto key = v;
            if (equal_class.find(key) == equal_class.end()) {
                std::set<Index_t> s;
                s.insert(k);
                equal_class[key] = s;
            }
            else 
                equal_class[key].insert(k);
        }
        }
        else {
            //dict = collapse_identical_elements(col_, row_);
                    //dict = collapse_identical_elements(row_, col_)
            auto rx = row_.template mutable_unchecked<1>();
            auto ry = col_.template mutable_unchecked<1>();
            size_t n_x = row_.shape(0);
            size_t n_y = col_.shape(0);

        //sort the raw array into a dict, where key is the element of array x
        // and value is the element of array y
        std::map<Index_t, std::set<Index_t>> dict;
        for (size_t i = 0; i < n_x; ++i) {
            auto key = rx(i);
            if (dict.find(key) == dict.end()) {
                //if no such key exists, we create a set and insert as the value
                std::set<Index_t> s;
                s.insert(ry(i));
                dict[key] = s;
            }
            else
                dict[key].insert(ry(i));
        }
        //based on [key, value] pairs of dict, combine its keys if they have the same value.
        //The combined keys will be stored in a new set as the value of equal_class.
        //Whereas the value of dict becomes the key of equal_class.
        
        for (auto&& [k, v] : dict) {
            auto key = v;
            if (equal_class.find(key) == equal_class.end()) {
                std::set<Index_t> s;
                s.insert(k);
                equal_class[key] = s;
            }
            else 
                equal_class[key].insert(k);
        }
        }
        nw::graph::edge_list<nw::graph::directed, Attributes...> g(0);
        g.open_for_push_back();
        for (auto&& [nodes, edges] : equal_class) {
            //because set is stored in sorted order.
            //So the minimum element of the set will reside in the first element.
            //Therefore, this first can be fetched with the help of set.begin() method.
            auto new_edge = *edges.begin();
            for (auto&& neighbor : nodes)
                g.push_back({new_edge, neighbor, 0});
        }
        g.close_for_push_back(false);
        return std::tuple{g, equal_class};
    }
}; //class NWhypergraph

}//namespace hypergraph
}//namespace nw