//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Funso Oje
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

#include "util/rapidxml-1.13/rapidxml.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace rapidxml;


static constexpr const char USAGE[] =
        R"(dblp.exe: dblp benchmark driver.
  Usage:
      dblp.exe (-h | --help)
      dblp.exe [-f FILE] [-o FILE] [--mtx FILE] [--sstart NUM] [--sstop NUM] [--sstep NUM] [-B NUM] [-t NUM] [--topn NUM] [--norm] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      -f FILE               dblp file path [default: ./dblp.xml.gz]
      -o FILE               matrix market output file path [default: ./dblp.mtx]
      --mtx FILE            mtx file for imdb [default: dblp.mtx]
      --sstart NUM          start s value of soverlap [default: 490]
      --sstop NUM           stop s value of soverlap [default: 1]
      --sstep NUM           step s value of soverlap [default: 1]
      -t NUM                number of threads [default: 1]
      -B NUM                number of bins [default: 32]
      --topn NUM            print top n centrality scores for each s [default: 3]
      --norm                normalize centrality scores [default: false]
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
    long num_bins = args["-B"].asLong();
    bool norm_Brande = args["--norm"].asBool();
    std::string mtx_file = args["--mtx"].asString();

    size_t threads = args["-t"].asLong();
    auto _ = set_n_threads(threads);
    size_t start_s_value = args["--sstart"].asLong();
    size_t stop_s_value  = args["--sstop"].asLong();
    size_t step_s_value  = args["--sstep"].asLong();

    size_t topn = args["--topn"].asLong();

    std::string input_file = args["-f"].asString();
    std::string output_file = args["-o"].asString();

    nw::graph::bi_edge_list<nw::graph::directedness::directed> edges(0);
    using vertex_id_t = vertex_id_t<decltype(edges)>;

    nw::util::timer m1("load mtx");

    edges = load_graph(mtx_file);

    m1.stop();
    std::cout << m1 << std::endl;

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

    std::map <std::string, vertex_id_t> names_map;
    std::map <vertex_id_t, std::string> names_map_transpose;

    std::map <std::string, vertex_id_t> titles_map;
    std::map <vertex_id_t, std::string> titles_map_transpose;

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

        //find the authors
        for (xml_node<> *author = article->first_node("author");
             author;
             author = author->next_sibling("author"))
        {
            auto author_name = author->value();
            //add name to map if not present
            auto it_name = names_map.find(author_name);
            if (it_name == names_map.end())
            {
                // not in map, add it
                names_map[author_name] = names_map_size;
                names_map_transpose[names_map_size] = author_name;
                it_name = names_map.find(author_name);
                ++names_map_size;
            }

        }
    }

    t3.stop();
    std::cout << t3 << "\n";

    edges.stream_stats();

    nw::util::timer t4("build biadjacencies");

    auto G = nw::graph::biadjacency<0>(edges);
    auto H = nw::graph::biadjacency<1>(edges);

    t4.stop();
    std::cout << t4 << "\n";


    for (size_t i_s = start_s_value; i_s >= stop_s_value; i_s = i_s - step_s_value)
    {

        size_t s_value = i_s;

        nw::util::timer t5("build s_overlap");

        auto &&degrees = H.degrees();

        // std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
        size_t M = edges.size();
        using linegraph_t = std::vector <std::vector<std::tuple < vertex_id_t, vertex_id_t>>>;
        linegraph_t two_graphs(num_bins);
        using container_t = std::unordered_map<size_t, size_t>;
        if (1 < s_value)
        {
            {
                map::to_two_graph_map_blocked<container_t>(
                        std::forward<linegraph_t>(two_graphs), H, G, degrees,
                        s_value, threads);
            }
        } else
        {
            {
                efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs),
                                                H, G, M / threads, 0, M);
            }
        }

        t5.stop();
        // std::cout << t5 << std::endl;

        //this dictionary keeps maping from old id (key), to new (value)
        std::unordered_map <vertex_id_t, vertex_id_t> old_to_new;
        nw::util::timer t6("squeeze s_overlap");
        auto &&s_overlap = create_edgelist_with_squeeze<nw::graph::directedness::undirected>(two_graphs, old_to_new);
        t6.stop();
        std::cout << t6 << std::endl;
        // std::cout << "linegraph size: " << s_overlap.size() << std::endl;

        //this dictionary keeps maping from new id (key), to old (value)
        std::unordered_map <vertex_id_t, vertex_id_t> new_to_old;
        for (auto&[oldid, newid]: old_to_new)
        {
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

        // std::cout << L.size() << std::endl;

        // only verify #cc in the result
        std::unordered_map <vertex_id_t, std::vector<vertex_id_t>> m;
        for (size_t i = 0, e = E.size(); i < e; ++i)
        {
            auto label = E[i];
            if (m.find(label) == m.end())
            {
                std::vector <vertex_id_t> cc;
                cc.push_back(i);
                m[label] = cc;
            } else
                m[label].push_back(i);
        }

        size_t numc = m.size();
        if (verbose)
        {
            for (auto&[k, v]: m)
            {
                std::cout << "Here are the " << s_value << "-connected components:" << std::endl;
                for (auto &i: v)
                {
                    auto oldid = new_to_old[i];
                    std::cout << names_map_transpose[oldid] << ", ";
                }
                std::cout << std::endl;
            }
        }


        using score_t = float;
        using accum_t = double;
        nw::util::timer t9("compute s-betweenness centrality");


        std::vector <score_t> bc = nw::graph::exact_brandes_bc<score_t, accum_t>(L, threads, std::execution::par_unseq,
                                                                                 std::execution::par_unseq,
                                                                                 norm_Brande);

        t9.stop();
        std::cout << s_value << ":" << t9 << std::endl;

        // Get highest three scores
        std::vector <size_t> indices(bc.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::partial_sort(indices.begin(), indices.begin() + topn, indices.end(),
                          [&](size_t A, size_t B) {
                              return bc[A] > bc[B];
                          });

        for (size_t tt = 0; tt < topn; ++tt)
        {
            std::cout << s_value << ":Top " << tt + 1 << ":" << names_map_transpose[new_to_old[indices[tt]]]
                      << ":" << bc[indices[tt]] << ":" << indices[tt] << ":" << E[indices[tt]] << std::endl;
        }

        std::cout << s_value << ":Stats:" << s_overlap.size() << ":" << numc << std::endl;


    }
    return 0;
}
