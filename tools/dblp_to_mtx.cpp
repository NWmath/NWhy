/**
 * @file dblp_to_mtx.cpp
 * 
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 * 
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * @author Funso Oje
 * @brief Converts the dblp.xml.gz into mtx format
 * @version 0.1
 * @date 2022-06-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <map>

#include <docopt.h>
#include "io/loader.hpp"
#include <nwgraph/edge_list.hpp>
#include <nwgraph/io/mmio.hpp>

#include <deque>
#include <queue>
#include <string>
#include <vector>

#include "nwgraph/adaptors/bfs_edge_range.hpp"
#include "nwgraph/adjacency.hpp"
#include "nwgraph/containers/compressed.hpp"

#include "util/rapidxml-1.13/rapidxml.hpp"

using namespace nw::hypergraph;
using namespace rapidxml;

static constexpr const char USAGE[] =
        R"(dblp2mtx.exe: create mtx file from dblp.xml.gz file.
  Usage:
      dblp2mtx.exe (-h | --help)
      dblp2mtx.exe [-f FILE] [-o FILE] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      -f FILE               dblp file path [default: ./dblp.xml.gz]
      -o FILE               matrix market output file path [default: ./dblp.mtx]
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

xml_document<> doc;
xml_node<> *root_node = NULL;

int main(int argc, char *argv[])
{
    std::vector <std::string> strings(argv + 1, argv + argc);
    auto args = docopt::docopt(USAGE, strings, true);

    bool verbose = args["--verbose"].asBool();
    bool debug = args["--debug"].asBool();

    std::string input_file = args["-f"].asString();
    std::string output_file = args["-o"].asString();

    nw::graph::bi_edge_list <nw::graph::directedness::directed> edges(0);
    using vertex_id_t = vertex_id_t<decltype(edges)>;

    // read the file
    nw::util::timer t1("reading xml file");
    std::cout << "Reading xml file ...\n";
    std::ifstream dblp_xml_stream(input_file);
    std::vector<char> dblp((std::istreambuf_iterator<char>(dblp_xml_stream)),
                           std::istreambuf_iterator<char>());
    dblp.push_back('\0');
    t1.stop();
    std::cout << t1 << "\n";

    // parse the file
    nw::util::timer t2("parsing xml file");
    doc.parse<0>(&dblp[0]);
    t2.stop();
    std::cout << t2 << "\n";

    nw::util::timer t3("build hypergraph");
    //root node
    root_node = doc.first_node("dblp");

    std::map<std::string, vertex_id_t> names_map;
    std::map<vertex_id_t, std::string> names_map_transpose;

    std::map<std::string, vertex_id_t> titles_map;
    std::map<vertex_id_t, std::string> titles_map_transpose;

//    edges.open_for_push_back();

    int titles_map_size{0};
    int names_map_size{0};
    for (xml_node<> *article = root_node->first_node("article");
         article;
         article = article->next_sibling())
    {
        //find the article title
        xml_node<> *title = article->first_node("title");

        if (!title) //ignore entries with no title
        {
            continue;
        }

        auto article_title = title->value();
        titles_map[article_title] = titles_map_size;
        titles_map_transpose[titles_map_size] = article_title;
        auto it_title = titles_map.find(article_title);

        //find the authors
        for (xml_node<> * author = article->first_node("author");
             author;
             author = author->next_sibling("author"))
        {
            auto author_name = author->value();
            //add name to map if not present
            auto it_name = names_map.find(author_name);
            if (it_name == names_map.end()) {
                // not in map, add it
                names_map[author_name] = names_map_size;
                names_map_transpose[names_map_size] = author_name;
                it_name = names_map.find(author_name);
                ++names_map_size;
            }

//            edges.push_back(it_title->second, it_name->second);

        }
        ++titles_map_size;
    }

//    edges.close_for_push_back();
    t3.stop();
    std::cout << t3 << "\n";
//    edges.stream_stats();

//    auto G = nw::graph::biadjacency<0>(edges);
//    auto H = nw::graph::biadjacency<1>(edges);

//    write_mm<0>(output_file, G, "general");

    return 0;

}