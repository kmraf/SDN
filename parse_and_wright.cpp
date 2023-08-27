#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include "other_algs.hpp"
#include "dist_func.hpp"

//topology file parsing (gml) and write in csv file
bool fill_matrix_write_csv(const std::string &if_name, const std::string &of_name, std::vector<std::vector<double>> &weight_matr) {

    //read file in string
    std::ifstream fin(if_name);
    if (!fin) return false;

    std::string file_data;
    char cur;
    while (fin.get(cur)) {
        file_data += cur;
    }
    fin.close();

    int node_cnt = 0;
    std::vector<std::vector<double>> coords;
    double longitude;
    double latitude;
    std::string s_longitude;
    std::string s_latitude;
    std::string s_label;
    size_t prev_pos = 0;

    //for labels, longs and lats
    std::vector<std::vector<std::string>> lab_lon_lat;

    while ((prev_pos = file_data.find("node [", prev_pos)) != std::string::npos) {

        s_label.clear();
        prev_pos = file_data.find("label ", prev_pos);
        prev_pos += std::string("label ").length();
        while (file_data[prev_pos] != '\n') {
            s_label += file_data[prev_pos];
            prev_pos++;
        }

        s_longitude.clear();
        prev_pos = file_data.find("Longitude ", prev_pos);
        prev_pos += std::string("Longitude ").length();
        while ((file_data[prev_pos] >= '0') && (file_data[prev_pos] <= '9') || (file_data[prev_pos] == '-') || (file_data[prev_pos] == '.')) {
            s_longitude += file_data[prev_pos];
            prev_pos++;
        }
        longitude = std::stod(s_longitude);

        s_latitude.clear();
        prev_pos = file_data.find("Latitude ", prev_pos);
        prev_pos += std::string("Latitude ").length();
        while ((file_data[prev_pos] >= '0') && (file_data[prev_pos] <= '9') || (file_data[prev_pos] == '-') || (file_data[prev_pos] == '.')) {
            s_latitude += file_data[prev_pos];
            prev_pos++;
        }
        latitude = std::stod(s_latitude);

        coords.push_back(std::vector<double>());
        coords[node_cnt].push_back(longitude);
        coords[node_cnt].push_back(latitude);

        lab_lon_lat.push_back(std::vector<std::string>());
        lab_lon_lat[node_cnt].push_back(s_label);
        lab_lon_lat[node_cnt].push_back(s_longitude);
        lab_lon_lat[node_cnt].push_back(s_latitude);

        node_cnt++;
    }
    s_latitude.clear();
    s_longitude.clear();

    //form matrix;
    weight_matr.resize(node_cnt);
    for (int i = 0; i < node_cnt; i++) {
        weight_matr[i].resize(node_cnt);
        for (int j = 0; j < node_cnt; j++) {
            weight_matr[i][j] = 0;
        }
    }

    //find edges in data
    prev_pos = 0;
    int src;
    int trg;
    std::string s_src;
    std::string s_trg;
    while ((prev_pos = file_data.find("edge [", prev_pos)) != std::string::npos) {

        s_src.clear();
        prev_pos = file_data.find("source ", prev_pos);
        prev_pos += std::string("source ").length();
        while ((file_data[prev_pos] >= '0') && (file_data[prev_pos] <= '9')) {
            s_src += file_data[prev_pos];
            prev_pos++;
        }
        src = std::stoi(s_src);

        s_trg.clear();
        prev_pos = file_data.find("target ", prev_pos);
        prev_pos += std::string("target ").length();
        while ((file_data[prev_pos] >= '0') && (file_data[prev_pos] <= '9')) {
            s_trg += file_data[prev_pos];
            prev_pos++;
        }
        trg = std::stoi(s_trg);

        weight_matr[src][trg] = 1;
        weight_matr[trg][src] = 1;
    }
    s_src.clear();
    s_trg.clear();


    //add weigths in mattrix
    double lat1;
    double lon1;
    double lat2;
    double lon2;

    for (int i = 0; i < node_cnt; i++) {
        for (int j = 0; j < node_cnt; j++) {
            if ((static_cast<int>(weight_matr[i][j]) == 1) && (j > i)) {
                lon1 = coords[i][0];
                lat1 = coords[i][1];
                lon2 = coords[j][0];
                lat2 = coords[j][1];

                weight_matr[i][j] = dist_from_coord(lat1, lon1, lat2, lon2);
            }
            if ((static_cast<int>(weight_matr[i][j]) == 1) && (j < i)) {
                weight_matr[i][j] = weight_matr[j][i];
            }
        }
    }
    coords.clear();

    //write in csv file
    std::ofstream fout(of_name);
    if (!fout) return false;

    fout << "Node1 (id)" << "," << "Node1 (label)" << "," << "Node1 (longitute)" << "," << "Node1 (latitude)" << ",";
    fout << "Node2 (id)" << "," << "Node2 (label)" << "," << "Node2 (longitute)" << "," << "Node2 (latitude)" << ",";
    fout << "Distance (km)" << "," << "Delay (mks)" << "\n";

    for (int i = 0; i < weight_matr.size(); i++) {
        for (int j = 0; j < weight_matr.size(); j++) {
            if (static_cast<int>(weight_matr[i][j])) {

                fout << std::to_string(i + 1) << ",";
                fout << lab_lon_lat[i][0] << ",";
                fout << lab_lon_lat[i][1] << ",";
                fout << lab_lon_lat[i][2] << ",";

                fout << std::to_string(j + 1) << ",";
                fout << lab_lon_lat[j][0] << ",";
                fout << lab_lon_lat[j][1] << ",";
                fout << lab_lon_lat[j][2] << ",";

                fout << std::to_string(weight_matr[i][j]) << ",";
                fout << std::to_string(4.8 * weight_matr[i][j]) << "\n";
            }
        }
    }
    fout.close();

    return true;
}

//function to write sptree
bool form_routes_file(const std::string &of_name, std::vector<std::vector<double>> &matr, int root) {

    const int n = matr.size();

    std::ofstream fout(of_name);
    if (!fout) return false;

    std::vector<double> dists;
    dists.resize(n);
    dijkstra_alg(matr, root, &dists);

    std::vector<std::vector<int>> routes;
    routes.resize(n);
    std::vector<double> costs(n, 0);

    //for double comp
    double eps = 0.000001;
    int cur_vert;
    int prev_vert;
    for (int i = 0; i < n; i++) {
        if (i == root) continue;

        cur_vert = i;
        while (cur_vert != root) {
            for (int j = 0; j < n; j++) {

                if (static_cast<int>(matr[cur_vert][j]) && fabs(dists[cur_vert] - (dists[j] + matr[cur_vert][j])) < eps) {
                    costs[i] += matr[cur_vert][j];
                    routes[i].insert(routes[i].begin(), cur_vert);
                    cur_vert = j;
                    break;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) {
        routes[i].insert(routes[i].begin(), root);
    }

    fout << "Node1 (id)" << "," << "Node2 (id)" << "," << "Path type" << "," << "Path" << ",";
    fout << "Delay (mks)" << "\n";
    for (int i = 0; i < n; i++) {
        if (i == root) continue;

        fout << std::to_string(root) << "," << std::to_string(i) << "," << "main" << ",";
        fout << '"' << "[";
        for (int j = 0; j < routes[i].size(); j++) {
            fout << std::to_string(routes[i][j]);

            if (j == routes[i].size() - 1) fout << " ";
            else fout << ", ";
        }
        fout << "]" << '"' << ",";
        fout << std::to_string(costs[i]) << "\n";
    }

    fout.close();
    return true;
}