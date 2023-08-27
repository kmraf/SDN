#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 
#include <set>
#include <algorithm>
#include <cctype>
#include <limits.h>

#include "parse_and_wright.hpp"
#include "other_algs.hpp"
#include "dist_func.hpp"

int main(int argc, char* argv[])
{
    if (argc != 3) return 1;

    int K = std::stoi(argv[2]);
    if ((K != 1) && (K != 2)) return 2;
    std::string if_name = argv[1];

    //read graph from file
    std::vector<std::vector<double>> weight_matr;
    const std::string of_name("Claranet_topo.csv");

    if (!fill_matrix_write_csv(if_name, of_name, weight_matr)) return 3;

    // k1
    //
    //
    if (K == 1) {

        //build sptrees
        std::vector<std::vector<double>> buf;
        buf.resize(weight_matr.size());
        for (int i = 0; i < buf.size(); i++) {
            buf[i].resize(buf.size());
        }

        std::vector<std::vector<double>> buf2;
        buf2.resize(weight_matr.size());
        for (int i = 0; i < buf2.size(); i++) {
            buf2[i].resize(buf2.size());
        }

        std::vector<std::vector<int>> edges;
        std::vector<int> cur_edge;
        for (int i = 0; i < weight_matr.size(); i++) {
            for (int j = i; j < weight_matr.size(); j++) {
                if (static_cast<int>(weight_matr[i][j])) {
                    cur_edge.push_back(i);
                    cur_edge.push_back(j);
                    edges.push_back(cur_edge);
                    cur_edge.clear();
                }
            }
        }

        int best_num = 2;
        //read current sptree in cur_buf
        std::vector<std::vector<double>> *cur_buf = &buf;

        double all_min_dist = INT_MAX;
        double cur_min_dist;
        int root;
        int cur_root;
        bool end;


        while (true) {

            if (best_num == 2) {
                cur_buf = &buf;
            }
            else {
                cur_buf = &buf2;
            }
            end = !gen_next_stree(*cur_buf, edges);
            if (end) break;

            //assign weights for stree
            for (int i = 0; i < (*cur_buf).size(); i++) {
                for (int j = 0; j < (*cur_buf).size(); j++) {
                    if (static_cast<int>((*cur_buf)[i][j]))
                        (*cur_buf)[i][j] = weight_matr[i][j];
                }
            }

            cur_min_dist = INT_MAX;
            double dij_res;
            for (int i = 0; i < (*cur_buf).size(); i++) {
                dij_res = dijkstra_alg((*cur_buf), i, nullptr);
                if (dij_res < cur_min_dist) {
                    cur_min_dist = dij_res;
                    cur_root = i;
                }
            }

            if (cur_min_dist < all_min_dist) {
                if (best_num == 2) {
                    best_num = 1;
                }
                else {
                    best_num = 2;
                }
                all_min_dist = cur_min_dist;
                root = cur_root;
            }

            //form_routes_file(rts_k1_fname, *cur_buf, root);
        }
        cur_buf == nullptr;

        std::vector<std::vector<double>> &best_stree = (best_num == 1) ? buf : buf2;
        const std::string rts_k1_fname = "Claranet_routes.csv";
        if (!form_routes_file(rts_k1_fname, best_stree, root)) return 4;

    }
    //mst and K2
    //
    //
    else {
        std::vector<std::vector<double>> mst;
        prim_mst(weight_matr, mst);

        int mst_min_dist = INT_MAX;
        int mst_root;
        double dij_res;
        for (int i = 0; i < mst.size(); i++) {
            dij_res = dijkstra_alg(mst, i, nullptr);
            if (dij_res < mst_min_dist) {
                mst_min_dist = dij_res;
                mst_root = i;
            }
        }

        const std::string rts_k2_fname = "Claranet_routes_K2_K1.csv";
        if (!form_routes_file(rts_k2_fname, mst, mst_root)) return 5;
    }
    
    return 0;
}