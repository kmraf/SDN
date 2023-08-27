#pragma once

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cctype>
#include <limits.h>

void print_matr(const std::vector < std::vector<double>> &matr);

bool next_combination(std::vector<int> &a, int n);

bool dfs(const std::vector<std::vector<double>> &matr, int vertex, std::set<int>&visited, int parent);
bool has_cycle(const std::vector<std::vector<double>> &matr);

bool gen_next_stree(std::vector<std::vector<double>> &buf, const std::vector<std::vector<int>> &edges);

int min_distance(int total_v, std::vector<double> &dist, std::vector<bool> &spt_set);
double dijkstra_alg(const std::vector<std::vector<double>> &matr, int src, std::vector<double> *dists);
void prim_mst(std::vector<std::vector<double>> matr, std::vector<std::vector<double>> &buf);