#include "../include/param.h"

Param::Param(){}

Param::Param(const Param &p): 
	generator_folder(p.generator_folder),
	model(p.model),
	debug(p.debug),
	max_iter(p.max_iter),
	EPS(p.EPS),
	sigma(p.sigma),
	beta_bar(p.beta_bar){}

Param::Param(boost::program_options::variables_map vm): 
	generator_folder{vm["generator_folder"].as<std::string>()},
	model{vm["model"].as<std::string>()},
	debug{vm["debug"].as<bool>()},
	max_iter{vm["max_iter"].as<int>()},
	EPS{vm["EPS"].as<double>()},
	sigma{vm["sigma"].as<double>()},
	beta_bar{vm["beta_bar"].as<double>()}{
}

Param& Param::operator=(const Param& p){
	if(this == &p)
		return *this;

	generator_folder = p.generator_folder;
	model = p.model;
	debug = p.debug;
	max_iter = p.max_iter;
	EPS = p.EPS;
	sigma = p.sigma;
	beta_bar = p.beta_bar;
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Param &p){
	out << "Parameters:\n";
	if(p.debug){
		out << "DEBUG MODE\n";
	}
	out << "Generator Folder: " << p.generator_folder << "\n";
	out << "Max Iter = " << p.max_iter << "\n";
	out << "EPS = " << p.EPS << "\n";
	out << "sigma = " << p.sigma << "\n";
	out << "beta_bar = " << p.beta_bar << "\n";
	return out;
}

Param::~Param(){
}