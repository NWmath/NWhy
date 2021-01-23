//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <unordered_map>
#include <fstream>
#include <docopt.h>
#include <edge_list.hpp>
#include <util/parallel_for.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "io/mmio.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hystats.exe: hypergraph stats driver.
  Usage:
      hystats.exe (-h | --help)
      hystats.exe [-f FILE...] [--deg-dist FILE...] [--output FILE...] [-D NUM] [--relabel NUM] [--direction DIR] [-d] [--log FILE] [--log-header]

  Options:
      -h, --help            show this screen
      -f FILE               edge list or matrix market input file paths
      --deg-dist FILE       output degree distribution of edges/nodes to FILE (must have same num of input files)
      --output FILE         output matrix market file (with/without relabeling by degree)
      -D NUM                specify column either [0](edge) or [1](node) [default: 0]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
)";
template<typename T>
auto counting_sort(std::vector<T>& arr) {
  std::unordered_map<T, std::size_t> m;
  for(auto& n : arr) {
    ++m[n];
  }
  return m;
}

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool debug   = args["--debug"].asBool();
  long idx     = args["-D"].asLong();
  long relabelidx     = args["--relabel"].asLong();

  std::vector<std::tuple<std::string, std::string, std::string>> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(std::make_tuple(file, "", ""));
  }
  auto&& dd = args["--deg-dist"].asStringList();
  for (std::size_t i = 0; i < dd.size(); ++i) {
    std::get<1>(files[i]) = dd[i];
  }
  auto&& op = args["--output"].asStringList();
  for (std::size_t i = 0; i < op.size(); ++i) {
    std::get<2>(files[i]) = op[i];
  }
  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& [input, dd_file, output] : files) {
    
    auto reader = [&](std::string file) {
      auto aos_a   = load_graph<directed>(file);
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        // Run relabeling. This operates directly on the incoming edglist.
        if (-1 != relabelidx) {
          nw::hypergraph::relabel_by_degree(hyperedges, hypernodes, idx, args["--direction"].asString());
        }
        auto hyperedge_degrees = hyperedges.degrees();
        auto hypernode_degrees = hyperedges.degrees();
        return std::tuple(hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees);
      }
      else {
        // Run relabeling. This operates directly on the incoming edglist.
        if (-1 != relabelidx) {
          std::cout << "relabeling edge_list by degree..." << std::endl;
          if (1 == relabelidx)
            nw::hypergraph::relabel_by_degree<1>(aos_a, args["--direction"].asString());
          else
            nw::hypergraph::relabel_by_degree<0>(aos_a, args["--direction"].asString());
        }
        adjacency<0> hyperedges(aos_a);
        adjacency<1> hypernodes(aos_a); 
        std::vector<index_t> hyperedge_degrees =  aos_a.degrees<0>();
        std::vector<index_t> hypernode_degrees =  aos_a.degrees<1>();
        return std::tuple(hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees);
      }
    };
    auto&&[ hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees ] = reader(input);

    //hyperedges
    std::size_t M = hyperedges.size();
    std::size_t hyperE_max_degree = 0;
    std::size_t hyperE_total_degree = 0;
    std::for_each(hyperedge_degrees.begin(), hyperedge_degrees.end(), [&](auto deg) {
        hyperE_max_degree = std::max(hyperE_max_degree, deg);
        hyperE_total_degree += deg;
    });
    double hyperE_avg_degree = hyperE_total_degree * 1.0 / M;
    //hypernodes
    std::size_t N = hypernodes.size();
    std::size_t hyperN_max_degree = 0;
    std::size_t hyperN_total_degree = 0;
    std::for_each(hypernode_degrees.begin(), hypernode_degrees.end(), [&](auto deg) {
        hyperN_max_degree = std::max(hyperN_max_degree, deg);
        hyperN_total_degree += deg;
    });
    double hyperN_avg_degree = hyperN_total_degree * 1.0 / N;
    std::cout << "num_hyperedges = " << M << " num_hypernodes = " << N << std::endl;
    std::cout << "d_edge_max= " << hyperE_max_degree << " d_node_max= " << hyperN_max_degree << std::endl;
    std::cout << "d_edge_total= " << hyperE_total_degree << " d_node_total= " << hyperN_total_degree << std::endl;
    std::cout << "d_edge_avg= " << hyperE_avg_degree << " d_node_avg= " << hyperN_avg_degree << std::endl;
    
    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    if ("" != dd_file) {
    std::ofstream deg_dist(dd_file);
    std::cout << "Writing degree distribution file" << std::endl;
    if(deg_dist.is_open()) {
      std::cout << dd_file <<  " " << idx << std::endl;
      if (0 == idx) {
        auto&& deg_counts = counting_sort<std::size_t>(hyperedge_degrees);
        for (auto&& [k, v] : deg_counts) {
          deg_dist << std::to_string(k) << " " << std::to_string(v) << std::endl;
        }
      }
      else if (1 == idx) {
        auto&& deg_counts = counting_sort<std::size_t>(hypernode_degrees);
        for (auto&& [k, v] : deg_counts) {
          deg_dist << std::to_string(k) << " " << std::to_string(v) << std::endl;
        }
      }
      else {
        std::cerr << "unrecognized flag" << std::endl;
      }
      deg_dist.close();
    }
    }
    if ("" != output) {
      std::cout << "Writing mtx file" << std::endl;
      if (0 == idx)
        write_mm_hy(output, hyperedges, M, N);
      else if (1 == idx)
        write_mm_hy(output, hypernodes, N, M);
      else
        std::cerr << "unrecognized flag" << std::endl;
    }
  } //for file

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("cc", file, times, header, "Time(s)");
  }

  return 0;
}
