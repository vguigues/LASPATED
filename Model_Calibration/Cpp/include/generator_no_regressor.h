#ifndef _GENERATOR_NO_REGRESSOR_H
#define _GENERATOR_NO_REGRESSOR_H

#include "main.h"


class GeneratorNoRegressor{
public:
	double x_max, y_max;
	int n_x, n_y;
	int nb_weeks, nb_years;
	double slot_duration;
	std::vector<double> durations;
	ulong T;
	ulong R;
	ulong C;
	ulong D;
	std::vector<std::vector<int>> groups;
	std::vector<int> which_group;
	std::vector<int> daily_obs;
	ulong nb_holidays_years;
	ulong nb_land_types;
	ulong nb_regressors;
	std::vector<std::pair<bool, int>> is_holidays;
	std::vector<Location> regions;
	bool full_neighbors;
	int neighbor_factor;
	bool read_covariates;
	bool use_simulation;

	double sigma, beta, beta_bar, weight;
	std::vector<double> weights;
	int max_iter;

	xt::xarray<int> nb_observations;
	xt::xarray<int> nb_arrivals;

	xt::xarray<double> regressors;
	xt::xarray<double> alpha;

	std::vector<std::vector<int>> neighbors;
	std::vector<int> type_region;
	xt::xarray<double> distance;
	xt::xarray<int> sample;
	xt::xarray<double> theoretical_lambda;
	xt::xarray<double> l_bounds;

	GeneratorNoRegressor();
	GeneratorNoRegressor(std::string calls_path, 
		std::string neighbors_path, std::string info_path);
	GeneratorNoRegressor(xt::xarray<int>& N, xt::xarray<int>& M, std::vector<double>& a_durations, 
		std::vector<std::vector<int>>& a_groups, std::vector<double>& a_weights, 
		xt::xarray<double>& alphas, xt::xarray<double>& a_distance, std::vector<int>& a_type, 
		std::vector<std::vector<int>>& a_neighbors);
	~GeneratorNoRegressor() = default;

	xt::xarray<double> oracle_gradient_model(xt::xarray<double>& x);
	double oracle_objective_model(xt::xarray<double>& x);
	// xt::xarray<double> projection_regressors(xt::xarray<double>& x_beta);
	// std::vector<double> projected_gradient_armijo_boundary(xt::xarray<double>& x_beta,
	// 	xt::xarray<double>& x_delta, double alpha, double sigma, double beta_bar);
	std::vector<double> projected_gradient_armijo_feasible(xt::xarray<double>& 
		x);


	CrossValidationResult cross_validation(double proportion, std::vector<double>& alphas, 
		std::vector<double>& group_weights);
	void test();
	void calibrate();
	bool is_neighbor(int r, int s);
	void write_cv_results(CrossValidationResult& cv_result);

	double average_difference(xt::xarray<double>& x);

	double mat_prod(xt::xarray<double>& a, xt::xarray<double>& b);

	void comp_wise_max(xt::xarray<double>& z, double eps);	

	void write_params(xt::xarray<double>& x);

	void print_var(xt::xarray<double>& x, std::string prefix);
};

xt::xarray<double> laspated_no_reg(xt::xarray<int>& N, xt::xarray<int>& M, std::vector<double>& a_durations, 
		std::vector<std::vector<int>>& a_groups, std::vector<double>& a_weights, 
		xt::xarray<double>& alphas, xt::xarray<double>& a_distance, std::vector<int>& a_type, 
		std::vector<std::vector<int>>& a_neighbors, xt::xarray<double>& x);


#endif