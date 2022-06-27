/**
 * @file hyper_betweenness_centrality.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#pragma once
#include <stack>
#include <queue>
#include <nwgraph/util/timer.hpp>
#include <nwgraph/util/AtomicBitVector.hpp>
#include <nwgraph/util/atomic.hpp>
#include <nwgraph/adjacency.hpp>

namespace nw {
namespace hypergraph {

template <class score_t = float, class accum_t = size_t, adjacency_list_graph GraphN, adjacency_list_graph GraphE>
std::tuple<std::vector<score_t>, std::vector<score_t>> betweenness_centrality(const GraphN& vertices, const GraphE& edges, bool normalize = true) {
    using vertex_id_t = GraphN::vertex_id_type;
    size_t nvertices = num_vertices(vertices, 0);
    size_t nedges = num_vertices(edges, 0); //or num_vertices(vertices, 1);
    std::vector<score_t> scoresN(nvertices, 0.0), scoresE(nedges, 0.0);
    std::stack<vertex_id_t> SE, SN;
    std::queue<vertex_id_t> QN, QE;
    std::vector<accum_t> npathsN(nvertices, 0.0), npathsE(nedges, 0.0);
    std::vector<score_t> depthsE(nedges, -1), depthsN(nvertices, -1);

    auto forward_propagation = [](
        const auto& g,
        std::stack<vertex_id_t>& S,
        std::queue<vertex_id_t>& Q,
        std::queue<vertex_id_t>& next,
        std::vector<std::list<size_t>>& paths,
        const std::vector<score_t>& depths_to_read,
        std::vector<score_t>& depths_to_write,
        const std::vector<accum_t>& npaths_to_read,
        std::vector<accum_t>& npaths_to_write
    ) -> void {
        //forward propagation for every element (node or hyperedge) in Q
        while (!Q.empty()) {
            vertex_id_t element = Q.front();
            Q.pop();
            S.push(element);
            for (auto&& e : g[element]){
                vertex_id_t neighbor = target(g, e);
                if (-1 == depths_to_write[neighbor]) {
                    //if neighbor has not been touched
                    depths_to_write[neighbor] = depths_to_read[element] + 1;
                    next.push(neighbor);
                }
                if (depths_to_write[neighbor] == depths_to_read[element] + 1)
                {
                    npaths_to_write[neighbor] += npaths_to_read[element];
                    paths[neighbor].push_back(element);
                }
            }
        }
    };
    for(vertex_id_t sourceE = 0; sourceE < nedges; ++sourceE) {
        npathsE.assign(nedges, 0);
        npathsN.assign(nvertices, 0);
        
        depthsE.assign(nedges, -1);
        depthsN.assign(nvertices, -1);
        std::vector<std::list<size_t>> pathsE(nedges), pathsN(nvertices);

        npathsE[sourceE] = 1;     
        depthsE[sourceE] = 0;
        QE.push(sourceE);

        while (QE.empty() && QN.empty()) {
            //forward propagation for every hyperedge in QE
            forward_propagation(edges, SE, QE, QN, pathsE, depthsE, depthsN, npathsE, npathsN);
            forward_propagation(vertices, SN, QN, QE, pathsN, depthsN, depthsE, npathsN, npathsE);
            /*
            while (!QE.empty()) {
                vertex_id_t hyperE = QE.front();
                QE.pop();
                SE.push(hyperE);
                for (auto&& e : edges[hyperE]){
                    vertex_id_t hyperN = target(edges, e);
                    if (-1 == depthsN[hyperN]) {
                        //if hyperN has not been touched
                        depthsN[hyperN] = depthsE[hyperE] + 1;
                        QN.push(hyperN);
                    }
                    if (depthsN[hyperN] == depthsE[hyperE] + 1)
                    {
                        npathsN[hyperN] += npathsE[hyperE];
                        pathsN[hyperN].push_back(hyperE);
                    }
                }
            }
            //forward propagation for every node in frontierN
            while(!QN.empty()) {
                vertex_id_t hyperN = QN.front();
                QN.pop();
                SN.push(hyperN);
                for (auto&& e : vertices[hyperN]){
                    vertex_id_t hyperE = target(vertices, e);
                    if (-1 == depthsE[hyperE]) {
                        //if hyperE has not been touched
                        depthsE[hyperE] = depthsN[hyperN] + 1;
                        QE.push(hyperE);
                    }
                    if (depthsE[hyperE] == depthsN[hyperN] + 1)
                    {
                        npathsE[hyperE] += npathsN[hyperN];
                        pathsE[hyperE].push_back(hyperN);
                    }
                }
            }
            */
        }//while

        std::vector<score_t> deltaE(nedges, 0.0);
        while(!SE.empty()) {
            vertex_id_t hyperE = SE.top();
            SE.pop();
            for (auto it = pathsE[hyperE].begin(); it != pathsE[hyperE].end(); ++it) {
                deltaE[*it] += static_cast<score_t>(npathsE[*it]) / static_cast<score_t>(npathsE[hyperE]) * (1 + deltaE[hyperE]);
            }
            if (hyperE != sourceE) {
                scoresE[hyperE] += deltaE[hyperE];
            }
        }
        std::vector<score_t> deltaN(nvertices, 0.0);
        while(!SN.empty()) {
            vertex_id_t hyperN = SN.top();
            SN.pop();
            for (auto it = pathsN[hyperN].begin(); it != pathsN[hyperN].end(); ++it) {
                deltaN[*it] += static_cast<score_t>(npathsN[*it]) / static_cast<score_t>(npathsN[hyperN]) * (1 + deltaN[hyperN]);
            }
            if (hyperN != sourceE) {
                scoresN[hyperN] += deltaN[hyperN];
            }
        }
    }
    return {scoresN, scoresE};
}

}//namespace hypergraph
}//namespace nw