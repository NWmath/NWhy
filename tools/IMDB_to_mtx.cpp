/**
 * @file IMDB_to_mtx.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#include <map>

#include <docopt.h>
#include "io/loader.hpp"
#include <nwgraph/edge_list.hpp>
#include <nwgraph/io/mmio.hpp>

#include "xtensor/xcsv.hpp"

using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(imdb2mtx.exe: create mtx file from imdb tsv files.
  Usage:
      imdb2mtx.exe (-h | --help)
      imdb2mtx.exe [--datafolder FILE] [--title FILE] [--name FILE] [--principal FILE] [--akas FILE] [--region FILE] [--styear NUM] [--enyear NUM] [-o FILE] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --datafolder FILE     folder path to all datafiles (title.basics.tsv, name.basics.tsv and title.principals.tsv) [default: .]
      --title FILE          movie title file path [default: title.basics.tsv]
      --name FILE           actor name file path [default: name.basics.tsv]
      --principal FILE      movie to actor file path [default: title.principals.tsv]
      --akas FILE           title attributes file [default: title.akas.tsv]
      --region FILE         filter titles by region, cannot be used with mtx [default: ]
      --styear NUM          start year filter, cannot be used with mtx [default: 1870]
      --enyear NUM          end year filter, cannot be used with mtx [default: 2030]
      -o FILE               matrix market output file path
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  std::string output_file = args["-o"].asString();

  long start_year= args["--styear"].asLong();
  long end_year  = args["--enyear"].asLong();
  
  nw::graph::bi_edge_list<nw::graph::directedness::directed> edges(0);
  using vertex_id_t = vertex_id_t<decltype(edges)>;

  std::string datafolder                = args["--datafolder"].asString();
  std::string title_basics_tsv_args     = args["--title"].asString();
  std::string name_basics_tsv_args      = args["--name"].asString();
  std::string title_principals_tsv_args = args["--principal"].asString();  
  std::string title_akas_tsv_args       = args["--akas"].asString();  

  std::string title_basics_tsv     = datafolder + "/" + title_basics_tsv_args;
  std::string name_basics_tsv      = datafolder + "/" + name_basics_tsv_args;
  std::string title_principals_tsv = datafolder + "/" + title_principals_tsv_args;
  std::string title_akas_tsv       = datafolder + "/" + title_akas_tsv_args;

  std::string region = args["--region"].asString();

  // load regions first to filter titles if needed
  std::map<std::string, vertex_id_t> region_map;
  bool check_region = false;

  if (! region.empty()){
    nw::util::timer r1("get regions");
    check_region = true;
    std::ifstream title_akas_stream(title_akas_tsv);
    auto title_akas = xt::load_csv<std::string>(title_akas_stream, '\t');
    auto akas_shp = title_akas.shape();

    for (vertex_id_t ia = 0; ia < akas_shp[0]; ++ia)
    {
      if (title_akas(ia, 1) != "1") {
        continue;
      }
      if (title_akas(ia, 3) != region) {
        continue;
      }
      region_map[title_akas(ia, 0)] = ia;
    }

    r1.stop();
    std::cout << r1 << std::endl;
  } 

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
      if (check_region){
        // if the region does not show, then it's not in that region.
        auto region_title = region_map.find(titles(i, 0));
        if (region_title == region_map.end()){
          continue;
        }
      }
      
      titles_map[titles(i, 0)] = i;
      titles_map_transpose[i] = titles(i, 0);
    }
  }

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

  std::ifstream title_principals_stream(title_principals_tsv);
  auto          title_principals = xt::load_csv<std::string>(title_principals_stream, '\t');
  auto          shp              = title_principals.shape();

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
  edges.stream_stats();


  auto G = nw::graph::biadjacency<0>(edges);
  auto H = nw::graph::biadjacency<1>(edges);

  write_mm<0>(output_file, G, "general");

  return 0;
}