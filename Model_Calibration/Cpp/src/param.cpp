#include "../include/param.h"

Param::Param(){}

Param::Param(const Param &p): 
	generator_folder(p.generator_folder),
	model(p.model),
	debug(p.debug),
	max_iter(p.max_iter),
	EPS(p.EPS),
	sigma(p.sigma),
	beta_bar(p.beta_bar),
	weights_file(p.weights_file),
	weights_list(p.weights_list){}

Param::Param(boost::program_options::variables_map vm): 
	generator_folder{vm["generator_folder"].as<std::string>()},
	model{vm["model"].as<std::string>()},
	debug{vm["debug"].as<bool>()},
	max_iter{vm["max_iter"].as<int>()},
	EPS{vm["EPS"].as<double>()},
	sigma{vm["sigma"].as<double>()},
	beta_bar{vm["beta_bar"].as<double>()}
	weights_file{vm["weights_file"].as<std::string>()},
	weights_list{vm["weights_list"].as<std::vector<double>>()}{
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
	weights_file = p.weights_file;
	weights_list = p.weights_list;
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