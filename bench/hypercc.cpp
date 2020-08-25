//
// This file is part of the Graph Standard Library (aka BGL17 aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include "Log.hpp"
#include "algorithms/connected_component.hpp"
#include "common.hpp"
#include "edge_list.hpp"
#include "util/AtomicBitVector.hpp"
#include "util/atomic.hpp"
#include "util/intersection_size.hpp"
#include <docopt.h>

using namespace nw::hypergraph::bench;

static constexpr const char USAGE[] =
    R"(bicc.exe: BGL17 hypergraph connected components benchmark driver.
  Usage:
      bicc.exe (-h | --help)
      bicc.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               input file paths (can have multiples)
      -n NUM                number of trials [default: 1]
      --relabel             relabel the graph or not
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --succession STR      successor/predecessor [default: successor]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";


template<class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph(ExecutionPolicy&& exec, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  life_timer _(__func__);
  //auto n_nbs = adjacency<0>(H);
  //auto e_nbs = adjacency<1>(H);
  edge_list<undirected> two_graph(0);
  two_graph.open_for_push_back();
  auto counter = 0;
  for (size_t i = 0; i < e_nbs.size(); ++i) {
    for (size_t j = i + 1; j < e_nbs.size(); ++j) {
      ++counter;
      size_t count = bgl17::intersection_size(e_nbs[i], e_nbs[j]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
      if (count >= s) {
        two_graph.push_back(i, j);    //         add (i,j) to new edge_list
      }
    }
  }
  std::cout << counter << " intersections performed" << std::endl;
  two_graph.close_for_push_back();
  return two_graph;
}

template<class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv2(ExecutionPolicy&& exec, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  life_timer _(__func__);
  //auto n_nbs = adjacency<0>(H);
  //auto e_nbs = adjacency<1>(H);
  edge_list<undirected> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  auto counter = 0;
  //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
  for (size_t i = 0; i < e_nbs.size(); ++i) { //O(n)
    auto hyperE = i;
    if (hyperedgedegrees[hyperE] < s) continue;
    std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) { //O(average degree of hyperedges)
      auto hyperN = std::get<0>(x);
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& y) { //O(average degree of hypernodes)
        //so we check compid of each hyperedge
        auto anotherhyperE = std::get<0>(y);
        if (hyperE > anotherhyperE) return;
        ++counter;
        //O(average degree of hyperedges)
        size_t count = bgl17::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
        if (count >= s) {
          two_graph.push_back(hyperE, anotherhyperE);    //         add (i,j) to new edge_list
        }
      });
    });
  }
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}

template<class HyperGraph>
auto base_two(HyperGraph& H, size_t s = 1) {
  return base_two(std::execution::seq, H, s);
}

template<class ExecutionPolicy, class HyperNode, class SGraph, class STranspose>
auto base_two(ExecutionPolicy&& exec, HyperNode& hypernodes, SGraph& s_adj, STranspose& s_tran_adj) {
  auto E = ccv1(exec, s_adj);//Afforest(s_adj, s_trans_adj);                 // 3) run whatever on new_adjacency
  auto nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  for (size_t hyperN = 0; hyperN < nhypernodes; ++hyperN) {
    //get the id of the first neighbor (hyperedge) of hypernode hyperN
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  }
  return std::tuple(N, E);
}

template<class ExecutionPolicy, class HyperGraph>
auto to_relabel_graph(ExecutionPolicy&& exec, HyperGraph& aos_a) {

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
  edge_list<undirected> relabel_graph(0);
  relabel_graph.open_for_push_back();
  //we relabel the smaller set (either hypernode or hyperedge)
  if (n_nbs < e_nbs) {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      node = node + e_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }
  else {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      edge = edge + n_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }

  relabel_graph.close_for_push_back();
  return relabel_graph;
}

template<class ExecutionPolicy, class HyperGraph>
auto relabelHyperCC(ExecutionPolicy&& exec, HyperGraph& aos_a) {
  auto relabel_g = to_relabel_graph(exec, aos_a);    // 1) find 2-graph corresponding to s-overlapped hyper edges
  //relabel_g. template lexical_sort_by<0>();
  //relabel_g.uniq();
  //relabel_g.remove_self_loops();

  //relabel_g.stream_edges();

  auto s_adj  = adjacency<0>(relabel_g);        // 2) convert new edge_list to new_adjacency
  //auto s_trans_adj = adjacency<1>(relabel_g);
  //s_adj.stream_indices();

  auto labeling   = //Afforest(s_adj, s_trans_adj); 
  ccv1(exec, s_adj);//

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
  std::vector<vertex_id_t> N, E;
  if (n_nbs < e_nbs) {
    E.assign(labeling.begin(), labeling.begin() + e_nbs);
    N.assign(labeling.begin() + e_nbs, labeling.end());
  }
  else {
    N.assign(labeling.begin(), labeling.begin() + n_nbs);
    E.assign(labeling.begin() + n_nbs, labeling.end());
  }
  return std::tuple(N, E);
}

/*
* baseline has been verified with the Python version
*/
template<typename Execution, typename Graph>
auto baseline(Execution exec, Graph& aos_a) {
  life_timer _(__func__);

  size_t num_hyperedges = aos_a.max()[0] + 1;    // number of hyperedges
  size_t num_hypernodes = aos_a.max()[1] + 1;    // number of hypernodes

  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> E(num_hyperedges, std::numeric_limits<vertex_id_t>::max());

  std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
    auto&& [edge, node] = elt;
    if (E[edge] == std::numeric_limits<vertex_id_t>::max()) E[edge] = edge;
    if (E[edge] == N[node]) return;
    if (N[node] == std::numeric_limits<vertex_id_t>::max()) {
      N[node] = E[edge];
    } else if (N[node] > E[edge]) {
      auto temp = N[node];
      std::replace(exec, N.begin(), N.end(), temp, E[edge]);
      std::replace(exec, E.begin(), E.end(), temp, E[edge]);
    } else if (N[node] < E[edge]) {
      auto temp = E[edge];
      std::replace(exec, N.begin(), N.end(), temp, N[node]);
      std::replace(exec, E.begin(), E.end(), temp, N[node]);
    }
  });
  // for each u in adjacency
  //   if (vertex is unlabeled) label[u] = u
  //   for each neighbor v of u {
  //     if label[u] < label[v] {
  //       label[v] = label[u]
  //     } else if (label[v] < label[u]) {
  //       label[u] = label[v]
  //      }
  //   }
  // }

  // for each hyperedge he in he_adjacency
  //   if (he is unlabeled) label[he] = he
  //   for each hypernode hn in he
  //     for each hyperedge hf in hn
  //       if label[hf] < label[he] update
  //       if label[he] < label[hf] update

  // for every hyper_edge in original list of hyperedges -- or use edge_range(adjacency)
  // get vertex (or edge) number
  // look corresponding component
  // push_back edge to comps[component].push_back(hyper_edge)

  // OR
  // sort aos_a according to component number
  // std::sort( ... , [])
  
  return std::tuple(N, E);
}

template<typename Execution, typename Graph, typename GraphN, typename GraphE>
auto svCC(Execution exec, Graph& aos_a, GraphN& cn, GraphE& ce) {
  life_timer _(__func__);

  size_t num_hyperedges = ce.max() + 1;    // number of hyperedges
  size_t num_hypernodes = cn.max() + 1;    // number of hypernodes
  std::vector<vertex_id_t> N(num_hypernodes), E(num_hyperedges);
  std::iota(E.begin(), E.end(), 0);
  std::iota(N.begin(), N.end(), 0);
  /*
  for (size_t k = 0; k < E.size(); ++k) {
    E[k] = k;
  }
  */

  bool changed = true;
  for (size_t num_iter = 0; num_iter < aos_a.size(); ++num_iter) {
    if (false == changed) break;
    changed = false;

    for (vertex_id_t edge = 0; edge < num_hyperedges; ++edge) {
      for (auto&& elt : ce[edge]) {
        auto&& [node] = elt;

        if (N[node] > E[edge]) {
          N[node] = E[edge];
          while (N[node] != N[N[node]]) {
            N[node] = N[N[node]];
          }
          changed = true;
        } else if (N[node] < E[edge]) {
          E[edge] = N[node];
          while (E[node] != E[E[node]]) {
            E[node] = E[E[node]];
          }
          changed = true;
        }
      }
    }

    for (vertex_id_t node = 0; node < num_hypernodes; ++node) {
      for (auto&& elt : cn[node]) {
        auto&& [edge] = elt;

        if (N[node] > E[edge]) {
          N[node] = E[edge];
          while (N[node] != N[N[node]]) {
            N[node] = N[N[node]];
          }
          changed = true;
        } else if (N[node] < E[edge]) {
          E[edge] = N[node];
          while (E[node] != E[E[node]]) {
            E[node] = E[E[node]];
          }
          changed = true;
        }
      }
    }
    /*
   for (vertex_id_t n = 0; n < num_hyperedges; n++) {
    while (E[n] != E[E[n]]) {
      E[n] = E[E[n]];
    }
  }

  for (vertex_id_t n = 0; n < num_hypernodes; n++) {
    while (N[n] != N[N[n]]) {
      N[n] = N[N[n]];
    }
  }
  */
  }    //for

  return std::tuple{N, E};
}


template<typename GraphN, typename GraphE>
auto lpCC(GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  bgl17::AtomicBitVector   visitedN(num_hypernodes);
  bgl17::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //initial node frontier includes every node
  std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);

  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          bgl17::writeMin(N[hyperN], labelE);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) frontierN.push_back(hyperN);
        }
      });
    });
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          bgl17::writeMin(E[hyperE], labelN);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
      });
    });
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
  }    //while

  return std::tuple{N, E};
}

inline bool updateAtomic(std::vector<vertex_id_t>& dest, std::vector<vertex_id_t>& source, std::vector<vertex_id_t>& prevDest,
                         vertex_id_t d, vertex_id_t s) {    //atomic Update
  auto origID = dest[d];
  return (bgl17::writeMin(dest[d], source[s]) && origID == prevDest[d]);
}
// Verifies CC result by performing a BFS from a vertex in each component
// - Asserts search does not reach a vertex with a different component label
// - If the graph is directed, it performs the search as if it was undirected
// - Asserts every vertex is visited (degree-0 vertex should have own label)
template<typename Execution, typename GraphN, typename GraphE>
auto bfsCC(Execution exec, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> prevE(num_hyperedges);
  std::vector<vertex_id_t> prevN(num_hypernodes);

  bgl17::AtomicBitVector   visitedN(num_hypernodes);
  bgl17::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //frontierE.resize(num_hyperedges);
  //initial node frontier includes every node
  std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);

  /*
  //or use a parallel for loop
  for (vertex_id_t i = 0; i < num_hyperedges; ++i) {
  //std::for_each(exec, std::begin(0), std::end(num_hypernodes), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  }
  //);
  */

  std::vector<size_t> feworkload, fnworkload;

  while (false == (frontierE.empty() && frontierN.empty())) {
    feworkload.push_back(frontierE.size());

    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(hyperedges[hyperE].begin(), hyperedges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          //updateAtomic(N, E, prevN, hyperN, hyperE);
          N[hyperN] = labelE;
          auto old  = N[hyperN];
          while (old == N[hyperN]) {
            auto s = bgl17::cas(N[hyperN], old, labelE);
            if (s) {
              if (visitedN.atomic_get(hyperN) == 0) {
                // visitedN.atomic_set(hyperN);
                frontierN.push_back(hyperN);
              }
              break;
            }
            old = N[hyperN];
            if (N[hyperN] < labelE) {
              break;
            }
          }
        }
        /*
          if (N[hyperN] == E[hyperE]) return;
          else if (E[hyperE] < N[hyperN]) {
            bgl17::writeMin(N[hyperN], E[hyperE]);
            if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) {
              frontierN.push_back(hyperN);
            }
          }
*/
      });
    });
    //reset bitmap for N
    visitedN.clear();
    std::for_each(exec, frontierE.begin(), frontierE.end(), [&](auto& i) {
      //all neighbors of hyperedges are hypernode
      prevE[i] = E[i];
    });
    frontierE.clear();
    fnworkload.push_back(frontierN.size());
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(hypernodes[hyperN].begin(), hypernodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          updateAtomic(E, N, prevE, hyperE, hyperN);
          //while (bgl17::writeMin(E[hyperE], labelN) && prevE[hyperN] == labelE);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
        /*
          if (E[hyperE] == N[hyperN]) return;
          else if (N[hyperN] < E[hyperE]) {
            bgl17::writeMin(E[hyperE], N[hyperN]);
            if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
              frontierE.push_back(hyperE);
            }
          }
          */
      });
    });
    //reset bitmap for E

    std::for_each(exec, frontierN.begin(), frontierN.end(), [&](auto& i) {
      //all neighbors of hyperedges are hypernode
      prevN[i] = N[i];
    });
    visitedE.clear();
    frontierN.clear();
  }    //while

  std::cout << "feworkload:";
  for (auto i : feworkload)
    std::cout << i << " ";
  std::cout << std::endl;
  std::cout << "fnworkload:";
  for (auto i : fnworkload)
    std::cout << i << " ";
  std::cout << std::endl;
  return std::tuple{N, E};
}


template<typename Execution, typename GraphN, typename GraphE>
auto lpCC(Execution exec, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  bgl17::AtomicBitVector   visitedN(num_hypernodes);
  bgl17::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //initial node frontier includes every node
  //std::iota(frontierE.begin(), frontierE.end(), 0);
  //std::iota(E.begin(), E.end(), 0);
  
  //or use a parallel for loop
  std::for_each(exec, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  });
  
  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          bgl17::writeMin(N[hyperN], labelE);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) frontierN.push_back(hyperN);
        }
      });
    });
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          bgl17::writeMin(E[hyperE], labelN);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
      });
    });
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
  }    //while
  return std::tuple{N, E};
}

template<typename Execution, typename GraphN, typename GraphE>
auto lpaNoFrontierCC(Execution exec, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  bgl17::AtomicBitVector visitedN(num_hypernodes);
  bgl17::AtomicBitVector visitedE(num_hyperedges);
  //std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  //frontierN.reserve(num_hypernodes);
  //frontierE.resize(num_hyperedges);
  //initial node frontier includes every node
  //std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);
  std::iota(N.begin(), N.end(), 0);
  /*
  //or use a parallel for loop
  for (vertex_id_t i = 0; i < num_hyperedges; ++i) {
  //std::for_each(exec, std::begin(0), std::end(num_hypernodes), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  }
  //);
  */
  auto nodes  = hypernodes.begin();
  auto edges  = hyperedges.begin();
  auto change = true;
  //while (false == (frontierE.empty() && frontierN.empty())) {
  while (false == change) {
    //std::for_each(exec, frontierE.begin(), frontierE.end(), [&](auto& hyperE) {

    for (vertex_id_t hyperE = 0; hyperE < num_hyperedges; ++hyperE) {
      //find each active hyperedge
      if (visitedE.atomic_get(hyperE)) {
        //push every hypernode's label to all neighbors of hyperedges
        auto labelE = E[hyperE];
        visitedE.atomic_set(hyperE);
        std::for_each(exec, edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
          auto hyperN = std::get<0>(x);
          auto labelN = N[hyperN];
          if (labelE == labelN) {
            //set to unactive
            visitedN.atomic_set(hyperN);
            return;
          } else if (labelE < labelN) {
            //updateAtomic(N, E, prevN, hyperN, hyperE);
            bgl17::writeMin(N[hyperN], labelE);
            //while (bgl17::writeMin(N[hyperN], labelE) && prevN[hyperN] == labelN);
            if (!visitedN.atomic_get(hyperN)) {
              visitedN.atomic_set(hyperN);
              change = true;
              //frontierN.push_back(hyperN);
            }
          } else {
            //once we found all the previous update are all false, we append all the false-modified
            //hypernodes to frontierN
            //safe to write without atomic action
            E[hyperE] = labelN;
            //bgl17::writeMin(E[hyperE], labelN);
            labelE = labelN;
            visitedE.atomic_set(hyperE);
            /*
            std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& y) {
              auto hyperfalseN = std::get<0>(y);
              if (!visitedN.atomic_get(hyperfalseN)){
                visitedN.atomic_set(hyperfalseN);
                frontierN.push_back(hyperfalseN);
              }
            });
            */
          }
        });
      }
    }
    //reset bitmap for N
    //visitedN.clear();
    //frontierE.clear();
    for (vertex_id_t hyperN = 0; hyperN < num_hypernodes; ++hyperN) {
      if (visitedN.atomic_get(hyperN)) {
        //std::for_each(exec, frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
        //all neighbors of hypernodes are hyperedges
        auto labelN = N[hyperN];

        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          auto labelE = E[hyperE];
          if (labelN == E[hyperE]) {
            visitedE.atomic_set(hyperE);
            return;
          } else if (labelN < labelE) {
            bgl17::writeMin(E[hyperE], labelN);
            //updateAtomic(E, N, prevE, hyperE, hyperN);
            //while (bgl17::writeMin(N[hyperN], labelE) && prevN[hyperN] == labelN);
            if (!visitedE.atomic_get(hyperE)) {
              visitedE.atomic_set(hyperE);
              change = true;
              //frontierE.push_back(hyperE);
            }
          } else {
            N[hyperE] = labelE;
            //bgl17::writeMin(N[hyperN], labelE);
            labelN = labelE;
            /*
            std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& y) {
              auto hyperfalseE = std::get<0>(y);
              if (!visitedE.atomic_get(hyperfalseE)){
                visitedE.atomic_set(hyperfalseE);
                frontierE.push_back(hyperfalseE);
              }
            });
            */
          }
        });
      }
    }
    //reset bitmap for E
    //visitedE.clear();
    //frontierN.clear();
  }    //while
  return std::tuple{N, E};
}

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    auto reader = [&](std::string file, bool verbose) {
      auto aos_a   = load_graph<directed>(file);
      auto hyperedgedegrees = aos_a.degrees<0>();

      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        //relabel the column with smaller size
        if (aos_a.max()[0] > aos_a.max()[1]) {
          auto hypernodedegrees = aos_a.degrees<1>();
          std::cout << "relabeling hypernodes..." << std::endl;
          aos_a.relabel_by_degree_bipartite<1>(args["--direction"].asString(), hypernodedegrees);
        }
        else {
          std::cout << "relabeling hyperedges..." << std::endl;
          aos_a.relabel_by_degree_bipartite<0>(args["--direction"].asString(), hyperedgedegrees);
        }
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      adjacency<0> hyperedges(aos_a);
      adjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << aos_a.max()[0] + 1 << " num_hypernodes = " << aos_a.max()[1] + 1 << std::endl;
      return std::tuple(aos_a, hyperedges, hypernodes, hyperedgedegrees);
    };

    auto&& graphs     = reader(file, verbose);
    auto&& aos_a      = std::get<0>(graphs);
    auto&& hyperedges = std::get<1>(graphs);
    auto&& hypernodes = std::get<2>(graphs);
    auto&& hyperedgedegrees = std::get<3>(graphs);

    edge_list<undirected> s_over;
    edge_list<undirected> s_over1 = to_two_graph(std::execution::seq, hyperedges, hypernodes);
    if (std::find(ids.begin(), ids.end(), 6) != ids.end())
    s_over = to_two_graphv2(std::execution::seq, hyperedges, hypernodes, hyperedgedegrees);    // 1) find 2-graph corresponding to s-overlapped hyper edges
    //s_over. template lexical_sort_by<0>();
    //s_over.uniq();
    //s_over.remove_self_loops();

  //s_over.stream_edges();

    auto s_adj  = adjacency<0>(s_over);        // 2) convert new edge_list to new_adjacency
    auto s_trans_adj = adjacency<1>(s_over);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        if (verbose) {
          std::cout << "version " << id << std::endl;
        }

        auto verifier = [&](auto&& result) {
          //TODO
          if (verbose) {
            //print_top_n(graph, comp);
          }
          std::map<vertex_id_t, edge_list<>> comps;
          auto N = std::get<0>(result);
          std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
            auto&& [edge, node] = elt;
            vertex_id_t key     = N[node];
            comps[key].push_back(elt);
          });

          for (auto&& j : comps) {
            auto& [k, v] = j;
            v.close_for_push_back();
          }
          std::cout << comps.size() << " components found" << std::endl;

          if (verify) {
            std::cerr << " v" << id << " failed verification for " << file << " using " << thread << " threads\n";
          }
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };


        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { return baseline(std::execution::seq, aos_a); });
              break;
            case 1:
              record([&] { return baseline(std::execution::par_unseq, aos_a); });
              break;
            case 2:
              record([&] { return lpCC(hypernodes, hyperedges); });
              break;
            case 3:
              record([&] { return svCC(std::execution::seq, aos_a, hypernodes, hyperedges); });
              break;
            case 4:
              record([&] { return bfsCC(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 5:
              record([&] { return lpCC(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 6:
              record([&] { return base_two(std::execution::seq, hypernodes, s_adj, s_trans_adj); });
              break;
            case 7:
              record([&] { return relabelHyperCC(std::execution::seq, aos_a); });
              break;
            case 8:
              record([&] { return relabelHyperCC(std::execution::par_unseq, aos_a); });
              break;
            default:
              std::cout << "Unknown version v" << id << "\n";
          }
        }
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("cc", file, times, header, "Time(s)");
  }

  return 0;
}
