#pragma once
#include <string>

bool fill_matrix_write_csv(const std::string &if_name, const std::string &of_name, std::vector<std::vector<double>> &weight_matr);
bool form_routes_file(const std::string &of_name, std::vector<std::vector<double>> &matr, int root);