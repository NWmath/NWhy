/**
 * @file adjoinbc.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#include <execution>
#include <unordered_set>
#include <nwgraph/edge_list.hpp>
#include <nwgraph/adjacency.hpp>
#include <nwgraph/algorithms/betweenness_centrality.hpp>
#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"
#include "io/mmio_hy.hpp"
#include "algorithms/adjoin_x.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace nw::graph;

static constexpr const char USAGE[] =
    R"(adjoinbc.exe: hypergraph betweenness centrality benchmark driver.
  Usage:
      adjoinbc.exe (-h | --help)
      adjoinbc.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

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

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto                     args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;

  std::vector<long> ids     = parse_ids(args["--version"].asStringList());
  std::vector<long> threads = parse_n_threads(args["THREADS"].asStringList());

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
    size_t num_realedges = 0, num_realnodes = 0;
    using vertex_id_t = uint32_t;
    auto&& [g, g_t, iperm]    = graph_reader_adjoin<nw::graph::directedness::undirected, vertex_id_t>(file, num_realedges, num_realnodes);
    std::cout << "num_hyperedges = " << num_realedges << " num_hypernodes = " << num_realnodes << std::endl;
    std::cout << "size of the merged adjacency = " << g.size() << std::endl;

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        auto verifier = [&](auto&& result) {
          auto&& [N, E] = result;
          //TODO use betweenness centrality score of hyperedges/vertices
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        using Graph = nw::graph::adjacency<0>;
        using ExecutionPolicy = decltype(std::execution::par_unseq);
        using score_t = float;
        using accum_t = double; 
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { 
                auto l = nw::graph::exact_brandes_bc<score_t, accum_t, Graph>(g, thread);
                return splitLabeling<ExecutionPolicy, score_t>(std::execution::par_unseq, l, num_realedges, num_realnodes);
              });
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
    log("bc", file, times, header, "Time(s)");
  }

  return 0;
}
