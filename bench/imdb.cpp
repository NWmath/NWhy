//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine, Xu Tony Liu
//

#include <unordered_set>
#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "s_overlap.hpp"
#include "util/slinegraph_helper.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "algorithms/s_connected_components.hpp"
#include <algorithms/betweenness_centrality.hpp>

#include <fstream>
#include <iostream>

#include <deque>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "xtensor/xcsv.hpp"

#include "bfs_edge_range.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(imdb.exe: imdb benchmark driver.
  Usage:
      imdb.exe (-h | --help)
      imdb.exe [--title FILE] [--name FILE] [--principal FILE] [--version ID...] [--loader-version ID] [-B NUM] [-s NUM...] [--relabel NUM] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernal version [default: 4]
      --title FILE          movie title file path
      --name FILE           actor name file path
      --principal FILE      movie to actor file path
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the hypergraph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       hypergraph relabeling direction - ascending/descending [default: ascending]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long num_bins= args["-B"].asLong();
  long loader_version = args["--loader-version"].asLong();
  long idx = args["--relabel"].asLong();
  std::string direction = args["--direction"].asString();

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());
  auto _ = set_n_threads(threads[0]);
  std::vector s_values= parse_ids(args["-s"].asStringList());
  if (s_values.empty()) s_values.push_back(1);

  std::string title_basics_tsv = args["--title"].asString();
  std::string name_basics_tsv = args["--name"].asString();
  std::string title_principals_tsv = args["--principal"].asString();

  nw::graph::edge_list<nw::graph::directedness::directed> edges(0);
  nw::util::timer               t0("load titles");
  std::ifstream                 title_basics_stream(title_basics_tsv);
  auto                          titles     = xt::load_csv<std::string>(title_basics_stream, '\t');
  auto                          titles_shp = titles.shape();
  std::map<std::string, vertex_id_t> titles_map;
  std::map<vertex_id_t, std::string> titles_map_transpose;
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < titles_shp[0]; ++i) {
    if (titles(i, 1) == "movie") {
      titles_map[titles(i, 0)] = i;
      titles_map_transpose[i] = titles(i, 0);
    }
  }
  t0.stop();
  std::cout << t0 << std::endl;

  nw::util::timer               t1("load names");
  std::ifstream                 name_basics_stream(name_basics_tsv);
  auto                          names     = xt::load_csv<std::string>(name_basics_stream, '\t');
  auto                          names_shp = names.shape();
  std::map<std::string, vertex_id_t> names_map;
  // this vector store the PrimaryName of the actors
  std::vector<std::string> names_map_transpose(names_shp[0]);
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < names_shp[0]; ++i) {
    names_map[names(i, 0)] = i;
    names_map_transpose[i] = names(i, 1);
  }
  t1.stop();
  std::cout << t1 << std::endl;

  nw::util::timer t2("load hypergraph");
  std::ifstream title_principals_stream(title_principals_tsv);
  auto          title_principals = xt::load_csv<std::string>(title_principals_stream, '\t');
  auto          shp              = title_principals.shape();

  t2.stop();
  std::cout << t2 << std::endl;

  nw::util::timer t3("build hypergraph");

  
  edges.open_for_push_back();

  size_t title_counter = 0;
  size_t name_counter  = 0;
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < shp[0]; ++i) {
    if (title_principals(i, 3) == "actor" || title_principals(i, 3) == "actress") {

      auto title = title_principals(i, 0);
      auto name  = title_principals(i, 2);
#if 0
      if (name == "nm0837064") {
        auto idx = titles_map[title];
        auto aa  = titles(idx, 2);
        std::cout << title << " " << aa << std::endl;
      }
#endif
      //if the title does not show, then it is not a movie.
      auto it_title = titles_map.find(title);
      if (it_title == titles_map.end()) {
	      //titles_map[title] = title_counter++;
        continue;
      }
      //same if the person is not an actor
      auto it_name = names_map.find(name);
      if (it_name == names_map.end()) {
	      //names_map[name] = name_counter++;
        continue;
      }
      edges.push_back(it_title->second, it_name->second);
    }
  }
  edges.close_for_push_back(false);

  t3.stop();
  std::cout << t3 << std::endl;
  edges.stream_stats();

  nw::util::timer t4("build biadjacencies");

  auto G = nw::graph::adjacency<0>(edges);
  auto H = nw::graph::adjacency<1>(edges);

  t4.stop();
  std::cout << t4 << std::endl;

  nw::util::timer t5("build s_overlap");

  index_t s = *std::min_element(s_values.begin(), s_values.end());
  auto&& degrees = H.degrees();
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  if (1 < s) {
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, H.size()), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          if (degrees[hyperE] < s) continue;
          std::map<size_t, vertex_id_t> K;
          for (auto&& [hyperN] : H[hyperE]) {
            for (auto&& [anotherhyperE] : G[hyperN]) {
              if (degrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (auto&& [anotherhyperE, val] : K) {
            if (val >= s)
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      },
      tbb::auto_partitioner());
  }
  t5.stop();
  std::cout << t5 << std::endl;

  //this dictionary keeps maping from old id (key), to new (value)
  std::unordered_map<vertex_id_t, vertex_id_t> old_to_new;
  nw::util::timer t6("squeeze s_overlap");
  auto&& s_overlap = create_edgelist_with_squeeze<nw::graph::directedness::undirected>(two_graphs, old_to_new);
  t6.stop();
  std::cout << t6 << std::endl;
  std::cout << "linegraph size: " << s_overlap.size() << std::endl;

  //this dictionary keeps maping from new id (key), to old (value)
  std::unordered_map<vertex_id_t, vertex_id_t> new_to_old;
  for (auto& [oldid, newid] : old_to_new) {
    new_to_old[newid] = oldid;
  }

  nw::util::timer t7("build s_overlap adjacency");

  auto L = nw::graph::adjacency<0>(s_overlap);
  t7.stop();
  std::cout << t7 << std::endl;

  // Kevin Bacon is nm0000102
  // David Suchet is nm0837064
  // Kyra Sedgwick is nm0001718
  nw::util::timer t8("compute s-connected components"); 
  auto E = ccv1(L);
  t8.stop();
  std::cout << t8 << std::endl;

  // only verify #cc in the result
  std::unordered_map<vertex_id_t, std::vector<vertex_id_t>> m;
  for (size_t i = 0, e = E.size(); i < e; ++i) {
    auto label = E[i];
    if (m.find(label) == m.end()) {
      std::vector<vertex_id_t> cc;
      cc.push_back(i);
      m[label] = cc;
    }
    else
      m[label].push_back(i);
  }
  size_t numc = 0;
  for (auto& [k, v] : m) {
    if (1 < v.size()) {
      ++numc;
      if (verbose) {
        std::cout << "Here are the " << s << "-connected components:" << std::endl;
        for (auto& i : v) {
          auto oldid = new_to_old[i];
          std::cout << names_map_transpose[oldid] << ", ";
        }
        std::cout << std::endl;
      }
    }
  }
  std::cout << m.size() << " components found" << std::endl;
  std::cout << numc << " non-singleton components found" << std::endl;

  using score_t = float;
  using accum_t = double;
  nw::util::timer t9("compute s-betweenness centrality"); 
  std::vector<score_t> bc = nw::graph::betweenness_brandes<decltype(L), score_t, accum_t>(L, false);
  t9.stop();
  std::cout << t9 << std::endl;
  
  int nbc = bc.size();
  float scale = 1.0;
  if (2 < nbc)
    scale /= ((nbc - 1) * (nbc - 2));
  for (size_t i = 0, e = nbc; i < e; ++i) {
    //normalized
    score_t score = bc[i] * scale;

    if (0.0 != score) {
      auto oldid = new_to_old[i];
      std::cout << names_map_transpose[oldid] << "(" << score << ")" << std::endl;
    }
  }
  return 0;
}
