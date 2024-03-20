#include "../include/param.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
Param::Param(){}

Param::Param(const Param &p): 
	generator_folder(p.generator_folder),
	method(p.method),
	model(p.model),
	debug(p.debug),
	max_iter(p.max_iter),
	EPS(p.EPS),
	sigma(p.sigma),
	beta_bar(p.beta_bar),
	cv_proportion(p.cv_proportion),
	weights_file(p.weights_file),
	weights_list(p.weights_list),
	type_proj_gradient(p.type_proj_gradient){}

Param::Param(boost::program_options::variables_map vm): 
	generator_folder{vm["generator_folder"].as<std::string>()},
	method{vm["method"].as<std::string>()},
	model{vm["model"].as<std::string>()},
	debug{vm["debug"].as<bool>()},
	max_iter{vm["max_iter"].as<int>()},
	EPS{vm["EPS"].as<double>()},
	sigma{vm["sigma"].as<double>()},
	beta_bar{vm["beta_bar"].as<double>()},
	cv_proportion{vm["cv_proportion"].as<double>()},
	weights_file{vm["weights_file"].as<std::string>()},
	type_proj_gradient{vm["type_proj_gradient"].as<int>()}{
	
	// fmt::print("Tokens: {}\n", vm["weights_list"].as< std::vector<std::string> >());
	try{
		weights_list = vm["weights_list"].as< std::vector<double> >();
	}catch(boost::wrapexcept<boost::bad_any_cast>& e){}

	if(weights_list.size() == 0  && weights_file == ""){
		if(model == "no_reg"){
			weights_list = {0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
			110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
			270,280,290,300,310,320,330,340,350,360,370,380,390,400};
		}else{
			weights_list = {0};
		}
	}else if(weights_file != "" && weights_list.size() == 0){
		std::ifstream weights_arq(weights_file, std::ios::in);
		while(!weights_arq.eof()){
			double w;
			weights_arq >> w;
			weights_list.push_back(w);
		}
		weights_arq.close();
	}
}

Param& Param::operator=(const Param& p){
	if(this == &p)
		return *this;

	generator_folder = p.generator_folder;
	model = p.model;
	method = p.method;
	debug = p.debug;
	max_iter = p.max_iter;
	EPS = p.EPS;
	sigma = p.sigma;
	beta_bar = p.beta_bar;
	cv_proportion = p.cv_proportion;
	weights_file = p.weights_file;
	weights_list = p.weights_list;
	type_proj_gradient = p.type_proj_gradient;
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