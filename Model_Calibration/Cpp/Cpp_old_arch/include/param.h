#ifndef _PARAM_H
#define _PARAM_H

#include "boost/program_options.hpp"
#include <fstream>
#include <iostream>

class Param{
public:
	Param();
	Param(const Param& p);
	explicit Param(boost::program_options::variables_map vm);

	std::string generator_folder, method, model;
	bool debug;
	int max_iter;
	double sigma, beta_bar, cv_proportion;
	double EPS = 0.001;
	std::string weights_file;
	std::vector<double> weights_list;
	int type_proj_gradient;


	Param& operator=(const Param& p);
	friend std::ostream& operator<<(std::ostream& out, const Param &s);
	~Param();
};

#endif