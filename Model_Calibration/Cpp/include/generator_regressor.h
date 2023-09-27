#ifndef _GENERATOR_REGRESSOR_H
#define _GENERATOR_REGRESSOR_H

#include "main.h"




class GeneratorRegressor{
public:
	GRBEnv& env;
	double x_max, y_max;
	int n_x, n_y;
	int nb_weeks, nb_years;
	double slot_duration;
	std::vector<double> durations;
	ulong T;
	ulong R;
	ulong C;
	ulong D;
	std::vector<std::vector<std::pair<int, int>>> groups;
	std::vector<std::vector<int>> which_group;
	std::vector<int> daily_obs;
	int nb_obs;
	ulong nb_holidays_years;
	ulong nb_land_types;
	ulong nb_regressors;
	std::vector<std::pair<bool, int>> is_holidays;
	std::vector<Location> regions;

	std::vector<std::vector<double>> u;

	double alpha, sigma, beta, beta_bar, weight;
	int max_iter;
	std::vector<double> weights;
	xt::xarray<int> nb_observations;
	xt::xarray<int> nb_arrivals;

	xt::xarray<double> regressors;

	std::vector<std::vector<int>> neighbors;
	std::vector<int> type;
	xt::xarray<double> distance;
	xt::xarray<std::vector<int>> sample;
	xt::xarray<double> theoretical_beta;
	xt::xarray<double> l_bounds;


    std::vector<double> obs_days;
    std::vector<double> obs_before;

	GeneratorRegressor(GRBEnv& env);
	GeneratorRegressor(GRBEnv& env, std::string calls_path, std::string neighbors_path, std::string info_path);
	GeneratorRegressor(GRBEnv& env, xt::xarray<int>& N, xt::xarray<int>& M, xt::xarray<double>& reg);
	~GeneratorRegressor() = default;

	xt::xarray<double> oracle_gradient_model2(xt::xarray<double>& x);
	double oracle_objective_model2(xt::xarray<double>& x);
	xt::xarray<double> projection_regressors(xt::xarray<double>& x);
	// std::vector<double> projected_gradient_armijo_boundary(xt::xarray<double>& x,
	// 	xt::xarray<double>& x, double alpha, double sigma, double beta_bar);
	std::vector<double> projected_gradient_armijo_feasible(xt::xarray<double>& 
		x_beta);


	CrossValidationResult cross_validation(double proportion, std::vector<double>& group_weights);
	void test();
	void calibrate();
	bool is_neighbor(int r, int s);

	double average_difference(xt::xarray<double>& x);

	double mat_prod(xt::xarray<double>& a, xt::xarray<double>& b);

	void comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps);	

	void write_params(xt::xarray<double>& x);

	double get_population(int r, std::vector<bool>& is_blue);

	double inner_integral(int i, int k);

	void print_vars(xt::xarray<double>& x, std::string prefix);
};

xt::xarray<double> laspated_reg(xt::xarray<int>& N, xt::xarray<int>& M, xt::xarray<double>& reg, xt::xarray<double>& x);

#endif