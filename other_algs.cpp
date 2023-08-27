#include <string>
#include <vector>
#include<set>  
#include <algorithm>
#include <cctype>
#include <limits.h>
#include <iostream>

//for debugging, bgg
void print_matr(const std::vector < std::vector<double>> &matr) {

    int node_cnt = matr.size();
    for (int i = 0; i < node_cnt; i++) {
        for (int j = 0; j < node_cnt; j++) {
            std::cout << matr[i][j] << "   ";
        }
        std::cout << std::endl;
    }
}

//gerenate combinations for spanning trees generation
bool next_combination(std::vector<int> &a, int n) {
    int k = (int)a.size();
    for (int i = k - 1; i >= 0; --i)
        if (a[i] < n - k + i + 1) {
            ++a[i];
            for (int j = i + 1; j < k; ++j)
                a[j] = a[j - 1] + 1;
            return true;
        }
    return false;
}

//dfs and has_cycle for spanning trees generation
bool dfs(const std::vector<std::vector<double>> &matr, int vertex, std::set<int>&visited, int parent) {
    visited.insert(vertex);
    for (int v = 0; v < matr.size(); v++) {
        if ((int)matr[vertex][v]) {
            if (v == parent)
                continue;
            if (visited.find(v) != visited.end())
                return true;
            if (dfs(matr, v, visited, vertex))
                return true;
        }
    }
    return false;
}

bool has_cycle(const std::vector<std::vector<double>> &matr) {
    std::set<int> visited;
    for (int v = 0; v < matr.size(); v++) {
        if (visited.find(v) != visited.end())
            continue;
        if (dfs(matr, v, visited, -1)) {
            return true;
        }
    }
    return false;
}

//the previous combination of edges is stored inside 
bool gen_next_stree(std::vector<std::vector<double>> &buf, const std::vector<std::vector<int>> &edges) {

    //sptree contains |V| - 1 edges
    static std::vector<int> prev_comb;

    static bool end = false;
    if (end) return false;

    int N = edges.size();
    int K = buf.size() - 1;
    if (prev_comb.size() == 0) {
        for (int i = 1; i <= K; i++) {
            prev_comb.push_back(i);
        }
    }

    while (true) {
        //prepare buf from prev call
        for (int i = 0; i < buf.size(); i++) {
            for (int j = 0; j < buf.size(); j++) {
                buf[i][j] = 0;
            }
        }

        int row;
        int col;

        for (int i = 0; i < K; i++) {
            row = edges[prev_comb[i] - 1][0];
            col = edges[prev_comb[i] - 1][1];
            buf[row][col] = 1;
            buf[col][row] = 1;
        }

        if (!has_cycle(buf)) {
            if (!next_combination(prev_comb, N)) {
                end = true;
            }
            return true;
        }
        else {
            if (!next_combination(prev_comb, N)) return false;
        }
    }
}

//function for dijkstra alg
int min_distance(int total_v, std::vector<double> &dist, std::vector<bool> &spt_set) {
    double min = INT_MAX;
    int min_index;
    for (int v = 0; v < total_v; v++)
        if (spt_set[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
    return min_index;
}

double dijkstra_alg(const std::vector<std::vector<double>> &matr, int src, std::vector<double> *dists) {
    std::vector<double> dist;
    dist.resize(matr.size());
    std::vector<bool> spt_set;
    spt_set.resize(matr.size());
    for (int i = 0; i < matr.size(); i++)
        dist[i] = INT_MAX, spt_set[i] = false;
    dist[src] = 0;
    for (int count = 0; count < matr.size() - 1; count++) {
        int u = min_distance(matr.size(), dist, spt_set);
        spt_set[u] = true;
        for (int v = 0; v < matr.size(); v++)
            if (!spt_set[v] && (matr[u][v] > 0) && dist[u] != INT_MAX && dist[u] + matr[u][v] < dist[v]) dist[v] = dist[u] + matr[u][v];
    }

    if (dists)  *dists = dist;
    return *std::max_element(dist.begin(), dist.end());
}

void prim_mst(std::vector<std::vector<double>> matr, std::vector<std::vector<double>> &buf) {

    int n = matr.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!static_cast<int>(matr[i][j]))
                matr[i][j] = INT_MAX;
        }
    }

    buf.resize(n);
    for (int i = 0; i < n; i++) {
        buf[i].resize(n);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            buf[i][j] = 0;
        }
    }

    std::vector<bool> used(n);
    std::vector<double> min_e(n, INT_MAX);
    std::vector<int> sel_e(n, -1);
    min_e[0] = 0;
    for (int i = 0; i < n; ++i) {
        int v = -1;
        for (int j = 0; j < n; ++j) {
            if (!used[j] && (v == -1 || min_e[j] < min_e[v]))
                v = j;
        }

        used[v] = true;
        if (sel_e[v] != -1) {
            int a = v;
            int b = sel_e[v];
            double c = matr[v][sel_e[v]];
            buf[v][sel_e[v]] = matr[v][sel_e[v]];
        }

        for (int to = 0; to < n; ++to) {
            if (matr[v][to] < min_e[to]) {
                min_e[to] = matr[v][to];
                sel_e[to] = v;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (static_cast<int>(buf[i][j]))
                buf[j][i] = buf[i][j];
        }
    }
}