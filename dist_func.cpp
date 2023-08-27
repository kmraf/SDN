#include <cmath> 

//use the haversine formula to find the distance
double dist_from_coord(double lat1, double lon1, double lat2, double lon2) {

    const double Pi = 3.14159265;
    const int R = 6371;

    double delta_lon = (lon2 - lon1)* Pi / 180;
    double delta_lat = (lat2 - lat1) * Pi / 180;
    double c_lat1 = cos(lat1 * Pi / 180);
    double c_lat2 = cos(lat2 * Pi / 180);
    double asin_arg_pow2 = sin(delta_lat / 2) * sin(delta_lat / 2) + c_lat2 * c_lat1 * sin(delta_lon / 2) * sin(delta_lon / 2);
    return R * 2 * atan2(sqrt(asin_arg_pow2), sqrt(1 - asin_arg_pow2));
}