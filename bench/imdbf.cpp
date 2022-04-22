//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine, Xu Tony Liu
//

#include <docopt.h>
#include <nwgraph/edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "s_overlap.hpp"
#include "util/slinegraph_helper.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "algorithms/s_connected_components.hpp"
#include <nwgraph/algorithms/betweenness_centrality.hpp>

#include <fstream>
#include <iostream>

#include <unordered_map>
#include <string>
#include <vector>

#include "xtensor/xcsv.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(imdb.exe: imdb benchmark driver.
  Usage:
      imdb.exe (-h | --help)
      imdb.exe [--datafolder FILE] [--title FILE] [--name FILE] [--principal FILE] [-B NUM] [--sstart NUM] [--sstop NUM] [--sstep NUM] [--norm] [--styear NUM] [--enyear NUM] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --datafolder FILE     folder path to all datafiles (title.basics.tsv, name.basics.tsv and title.principals.tsv) [default: .]
      --title FILE          movie title file path [default: title.basics.tsv]
      --name FILE           actor name file path [default: name.basics.tsv]
      --principal FILE      movie to actor file path [default: title.principals.tsv]
      --styear NUM          start year filter [default: 1870]
      --enyear NUM          end year filter [default: 2030]
      -B NUM                number of bins [default: 32]
      --sstart NUM          s value of soverlap [default: 130]
      --sstop NUM           s value of soverlap [default: 1]
      --sstep NUM           s value of soverlap [default: 1]
      --norm                normalize centrality scores [default: false]
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

  bool norm_Brande = args["--norm"].asBool();

  std::vector threads = parse_n_threads(args["THREADS"].asStringList());
  auto _ = set_n_threads(threads[0]);
  size_t start_s_value = args["--sstart"].asLong();
  size_t stop_s_value  = args["--sstop"].asLong();
  size_t step_s_value  = args["--sstep"].asLong();

  long start_year= args["--styear"].asLong();
  long end_year  = args["--enyear"].asLong();

  nw::graph::bi_edge_list<nw::graph::directedness::directed> edges(0);
  using vertex_id_t = vertex_id_t<decltype(edges)>;

  std::string datafolder                = args["--datafolder"].asString();
  std::string title_basics_tsv_args     = args["--title"].asString();
  std::string name_basics_tsv_args      = args["--name"].asString();
  std::string title_principals_tsv_args = args["--principal"].asString();  

  std::string title_basics_tsv     = datafolder + "/" + title_basics_tsv_args;
  std::string name_basics_tsv      = datafolder + "/" + name_basics_tsv_args;
  std::string title_principals_tsv = datafolder + "/" + title_principals_tsv_args;

  // std::string title_basics_tsv = args["--title"].asString();
  // std::string name_basics_tsv = args["--name"].asString();
  // std::string title_principals_tsv = args["--principal"].asString();  

  nw::util::timer               t0("load titles");
  std::ifstream                 title_basics_stream(title_basics_tsv);
  auto                          titles     = xt::load_csv<std::string>(title_basics_stream, '\t');
  auto                          titles_shp = titles.shape();
  std::map<std::string, vertex_id_t> titles_map;
  std::map<vertex_id_t, std::string> titles_map_transpose;
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < titles_shp[0]; ++i) {
    if (titles(i, 1) == "movie") {
      if (start_year > 1870 or end_year < 2030) {
        if (titles(i, 5) == "\\N" ){
          continue;
        }
        if (std::stol(titles(i, 5)) < start_year or std::stol(titles(i, 5)) > end_year) {
          continue;
        }
      }
      titles_map[titles(i, 0)] = i;
      titles_map_transpose[i] = titles(i, 0);
    }
  }
  t0.stop();
  std::cout << t0 << std::endl;
  std::cout << titles_map.size() << " titles loaded" << std::endl;

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
  //nw::graph::edge_list<nw::graph::directedness::directed> edges(0);
  edges.open_for_push_back();
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < shp[0]; ++i) {
    if (title_principals(i, 3) == "actor" || title_principals(i, 3) == "actress") {
      auto title = title_principals(i, 0);
      auto name  = title_principals(i, 2);
      //if the title does not show, then it is not a movie.
      auto it_title = titles_map.find(title);
      if (it_title == titles_map.end()) {
        continue;
      }
      //same if the person is not an actor
      auto it_name = names_map.find(name);
      if (it_name == names_map.end()) {
        continue;
      }
      edges.push_back(it_title->second, it_name->second);
    }
  }
  edges.close_for_push_back();
  t3.stop();
  std::cout << t3 << std::endl;
  edges.stream_stats();

  nw::util::timer t4("build biadjacencies");

  auto G = nw::graph::biadjacency<0>(edges);
  auto H = nw::graph::biadjacency<1>(edges);

  t4.stop();
  std::cout << t4 << std::endl;


  
  for (size_t i_s = start_s_value; i_s >= stop_s_value; i_s = i_s - step_s_value) {

    size_t s_value = i_s;

    nw::util::timer t5("build s_overlap");

    auto&& degrees = H.degrees();
    
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    if (1 <= s_value) {
      tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, H.size()), [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();    
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
            if (degrees[hyperE] < s_value) continue;
            std::unordered_map<size_t, size_t> K;
            for (auto&& [hyperN] : H[hyperE]) {
              for (auto&& [anotherhyperE] : G[hyperN]) {
                if (degrees[anotherhyperE] < s_value) continue;
                if (hyperE < anotherhyperE) ++K[anotherhyperE];
              }
            }

            for (auto&& [anotherhyperE, val] : K) {
              if (val >= s_value) {
                two_graphs[worker_index].push_back(
                    std::make_tuple<vertex_id_t, vertex_id_t>(
                        std::forward<vertex_id_t>(hyperE),
                        std::forward<vertex_id_t>(anotherhyperE)));
                // std::cout << names_map_transpose[hyperE] << " and "
                // << names_map_transpose[anotherhyperE] << " has "
                // << val << " collaborations" << std::endl;
              }
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
          std::cout << "Here are the " << s_value << "-connected components:" << std::endl;
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
    std::vector<score_t> bc = nw::graph::brandes_bc<decltype(L), score_t, accum_t>(L, norm_Brande);
    t9.stop();
    std::cout << t9 << std::endl;

    // float nbc = bc.size();
    // float scale = 1.0;
    // if (2 < nbc)
    //   scale /= ((nbc - 1) * (nbc - 2));

    // Get highest three scores
    std::vector<size_t> indices(bc.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::partial_sort(indices.begin(), indices.begin() + 3, indices.end(),
          [&](size_t A, size_t B) {
            return bc[A] > bc[B];
          });

    std::cout << s_value << ":Top 1:" << names_map_transpose[new_to_old[indices[0]]] << ":" << bc[indices[0]] << ":" << s_overlap.size() << ":" << numc << std::endl;
    std::cout << s_value << ":Top 2:" << names_map_transpose[new_to_old[indices[1]]] << ":" << bc[indices[1]] << ":" << s_overlap.size() << ":" << numc << std::endl;
    std::cout << s_value << ":Top 3:" << names_map_transpose[new_to_old[indices[2]]] << ":" << bc[indices[2]] << ":" << s_overlap.size() << ":" << numc << std::endl;    
    // // highest score
    // score_t hscore = 0.0;
    // std::string name;
    // for (size_t i = 0, e = nbc; i < e; ++i) {
    //   //normalized
    //   score_t score = bc[i] * scale;

    //   if (score > hscore) {
    //     hscore = score;
    //     auto oldid = new_to_old[i];
    //     name = names_map_transpose[oldid];
    //     // std::cout << names_map_transpose[oldid] << "(" << score << ")" << std::endl;
    //   }
    // }

    // //print highest score
    // if (hscore > 0) {
    //   std::cout << s_value << ":" << name << ":" << hscore << ":" << s_overlap.size() << ":" << m.size() << ":" << numc << std::endl;
    // }
  }
  return 0;
}
