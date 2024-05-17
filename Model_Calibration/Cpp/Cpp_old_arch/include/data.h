#ifndef _DATA_H
#define _DATA_H
#include "main.h"



class DataNoRegressor{
public:
    double x_max, y_max;
	int n_x, n_y;
	int nb_weeks, nb_years;
    std::vector<double> durations;
    ulong T;
	ulong R;
	ulong C;
    xt::xarray<int> nb_observations;
    xt::xarray<int> nb_arrivals;

    

};


class DataRegressor{
public:

};

#endif