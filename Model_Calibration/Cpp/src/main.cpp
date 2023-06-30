#include "../include/generator_no_regressor.h"
#include "../include/generator_regressor.h"
#include "../include/main.h"

//ILOSTLBEGIN
namespace po = boost::program_options;
std::string config_file;
Param g_params;

int main(int argc, char* argv[]){
	srand(time(NULL));
	
	g_params = load_params(argc,argv);

	std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);


	if(g_params.model == "no_reg"){
		// GeneratorNoRegressor gen;
		GeneratorNoRegressor gen(calls_path, neighbors_path, info_path);
		auto t0 = std::chrono::high_resolution_clock::now();
		gen.test();
		auto dt = std::chrono::high_resolution_clock::now();
		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
			/ pow(10,9);
		fmt::print("Run time final = {}\n", run_time_old);
	}else if(g_params.model == "reg"){
		GRBEnv env;
		GeneratorRegressor gen(env);
		// GeneratorRegressor gen(env, calls_path, neighbors_path, info_path);
		auto t0 = std::chrono::high_resolution_clock::now();
		gen.test();
		auto dt = std::chrono::high_resolution_clock::now();
		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
			/ pow(10,9);
		fmt::print("Run time final = {}\n", run_time_old);
	}

}


void version(){
	std::cout << "Spatio Temporal Calibration of medical emergency events.\n" <<
		"Version number 0.1 - Developed by: Victor Hugo (UFRJ-PESC-LABOTIM).\n";
}

Param load_params(int argc, char** argv){
	// Declare a group of options that will be 
	// allowed only on command line
	po::options_description generic("Generic Options");
	generic.add_options()
    	("version,v", "Prints version message.")
    	("help,h", "Prints help message.")
    	("file,f", po::value<std::string>(&config_file),
                  "configuration file path.")
    ;
    // Declare a group of options that will be 
	// allowed both on command line and in
	// config file
	po::options_description config("Configuration");
	config.add_options()
	    ("generator_folder,G", po::value<std::string>()->default_value("calibration/Rect10x10"),
	    	"scenario generator file path.")
		("model", po::value<std::string>()->default_value("no_reg"),
			"Type of data. reg|no_reg are data with and without regressors.")
	    ("debug,d", po::value<bool>()->default_value(false),
	    	"Runs in debug mode. default = false")
		("max_iter", po::value<int>()->default_value(100),
			"Max number of iterations in armijo procedures. default = 100")
		("EPS", po::value<double>()->default_value(0.001),
			"Tolerance")
		("sigma", po::value<double>()->default_value(0.5),
			"Sigma parameter for armijo procedures. default = 0.5")
		("beta_bar", po::value<double>()->default_value(1),
			"Beta_bar parameter for armijo procedures. default = 1.0")
		("weights_file", po::value<std::string>()->default_value(""),
			"path for file containing weights to be used. default = \"\"")
		("weights_list", po::value<std::vector<double>>()->multitoken()->default_value({1}),
			"list of weights to be used. default = {1}")
	;

	po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
	config_file_options.add(config);

	po::options_description visible("Allowed Options");
	visible.add(generic).add(config);

	po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).run(), vm);
    notify(vm);

	std::ifstream ifs(config_file.c_str());
    if (vm.count("file") && !ifs){
        std::cout << "Unable to open file: " << config_file << "\n";
        exit(1);
    }
    else if(vm.count("file") && ifs){
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }

    if (vm.count("help")) {
        std::cout << visible << "\n";
    }
    if (vm.count("version")) {
    	//TODO: ler versão do SVN?
        version();
    }

    return Param(vm);
}


// double travel_time = osrm_helper.get_duration(std::make_pair(-22.900434, -43.277936),
// 	std::make_pair(-22.859857, -43.307751));
// std::cout << "Travel Time OSRM: " << travel_time << "\n";


	// std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);


	// GRBEnv env;
	// GeneratorRegressor gen(env, calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_old = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T regressor old(s): {}\n", run_time_quick);

	// QuickGeneratorRegressor qr_gen(env,calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// qr_gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_quick = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T regressor quick(s): {}\n", run_time_quick);

	// Generator gen(calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_old = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T No regressor (s): {}\n", run_time_old);


	// QuickGenerator q_gen(calls_path,neighbors_path, info_path);
	// auto t0 = std::chrono::high_resolution_clock::now();
	// q_gen.test();
	// auto dt = std::chrono::high_resolution_clock::now();
	// double run_time_quick = std::chrono::duration_cast<chrono::nanoseconds>(dt - t0).count() 
	// 	/ pow(10,9);

	// fmt::print("T no regressor quick(s): {}\n", run_time_quick);
	// fmt::print("T old(s): {}, T new(s): {}\n", run_time_old, run_time_quick);




	
	// GRBEnv env;
	// GeneratorRegressor gen(env, calls_path,neighbors_path, info_path);
	// gen.test();

	// fmt::print("Finished instance {}\n",g_params.generator_folder);