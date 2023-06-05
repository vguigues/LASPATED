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

	std::string generator_folder, model;
	bool debug;
	int max_iter;
	double sigma, beta_bar;
	double EPS = 0.001;


	Param& operator=(const Param& p);
	friend std::ostream& operator<<(std::ostream& out, const Param &s);
	~Param();
};

#endif