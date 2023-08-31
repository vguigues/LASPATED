#include "../include/generator_no_regressor.h"
#include "../include/generator_regressor.h"
#include "../include/main.h"

//ILOSTLBEGIN
namespace po = boost::program_options;
std::string config_file;
Param g_params;


void test_laspated_reg(){
	using namespace std;

	int x_max = 10; int y_max = 10;
	int n_x = 10; int n_y = 10;
	ulong R = n_x*n_y;
	ulong C = 1;
	ulong T = 4;
	ulong D = 15;
	ulong nb_regressors = 3;

	int nb_weeks = 20;
	int nb_years = floor(nb_weeks/52);
	int nb_obs = nb_weeks*7;

	xt::xarray<double> theoretical_beta = xt::zeros<double>({C,D,T,nb_regressors});

	for(int d = 0; d < 7; ++d){
		theoretical_beta(0,d,1,0) = 0.05;
		theoretical_beta(0,d,3,0) = 0.05;

		theoretical_beta(0,d,0,1) = 6;
		theoretical_beta(0,d,1,1) = 18;
		theoretical_beta(0,d,2,1) = 6;
		theoretical_beta(0,d,3,1) = 18;
		
		theoretical_beta(0,d,0,2) = 3;
		theoretical_beta(0,d,1,2) = 6;
		theoretical_beta(0,d,2,2) = 3;
		theoretical_beta(0,d,3,2) = 6;
	}

	for(int d = 7; d < 15; ++d){
		theoretical_beta(0,d,1,0) = 0.1;
		theoretical_beta(0,d,3,0) = 0.1;
		
		theoretical_beta(0,d,0,1) = 12;
		theoretical_beta(0,d,1,1) = 36;
		theoretical_beta(0,d,2,1) = 12;
		theoretical_beta(0,d,3,1) = 36;
		
		theoretical_beta(0,d,0,2) = 6;
		theoretical_beta(0,d,1,2) = 12;
		theoretical_beta(0,d,2,2) = 6;
		theoretical_beta(0,d,3,2) = 12;
	}

	xt::xarray<double> regressors = xt::zeros<double>({nb_regressors, R});
	std::default_random_engine gen;
	std::uniform_real_distribution<double> rnd(0,1);
    int count = 0;
    for(int i = 0; i < n_y / 2; ++i){
        for(int j = 0; j < n_x/2; ++j){
            regressors(0, count) = 50 + 50*rnd(gen);
			regressors(1, count) = 0.5;
			regressors(2, count) = 0.25;
            ++count;
        } 

        for(int j = 0; j < n_x/2; ++j){
            regressors(0, count) = 50*rnd(gen);
			regressors(1, count) = 0.25;
			regressors(2, count) = 0.5;
            ++count;
        }
    }

    for(int i = n_y/2; i < n_y; ++i){
        for(int j = 0; j < n_x/2; ++j){
            regressors(0, count) = 50*rnd(gen);
			regressors(1, count) = 0.25;
			regressors(2, count) = 0.5;
            ++count;
        } 

        for(int j = 0; j < n_x/2; ++j){
            regressors(0, count) = 50 + 50*rnd(gen);
			regressors(1, count) = 0.5;
			regressors(2, count) = 0.25;
            ++count;
        }
    }

	xt::xarray<vector<int>> sample = xt::zeros<vector<int>>({C,D,T,R});
	xt::xarray<int> nb_observations = xt::zeros<int>({C,D,T,R});
	xt::xarray<int> nb_arrivals = xt::zeros<int>({C,D,T,R});

	int max_obs = 0;
	for(int index = 0; index < nb_weeks*7; ++index){
		int day = index % 7;
		for(int c = 0; c < C; ++c){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rate = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += theoretical_beta(c,day,t,j)*regressors(j,r);
					}
					poisson_distribution<int> pd(rate);
					int this_nb_arrival = pd(gen);
					sample(c,day,t,r).push_back(this_nb_arrival);
                    ++nb_observations(c,day,t,r);
					nb_arrivals(c,day,t,r) += this_nb_arrival;
				}
			}
		}
	}

	xt::xarray<double> lambda0 = xt::ones<double>({C,D,T,nb_regressors});
	xt::xarray<double> lambda = laspated_reg(nb_observations, nb_arrivals, regressors, lambda0);
	ofstream arq("x_reg.txt", std::ios::out);
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rate = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += lambda(c,d,t,j)*regressors(j,r);
					}
					arq << c << " " << d << " " << t << " " << r << " " << rate << "\n";
				}
			}
		}
	}
	arq.close();
	fmt::print("Lambdas saved at x_reg.txt\n");
}

void test_laspated_no_reg(){
	using namespace std;

	int x_max = 10;
	int y_max = 10;
	int n_x = 10;
	int n_y = 10;
	int nb_groups = 4;
	ulong R = n_x*n_y;
	ulong C = 1;
	ulong T = 4*7;


	vector<vector<int>> groups = vector<vector<int>>(nb_groups,vector<int>());
	for(int i = 0; i < nb_groups; ++i){
		groups[i].push_back(i);
	}

	for(int i = 0; i < nb_groups; ++i){
		for(int j = 0; j < 6; ++j){
			groups[i].push_back(groups[i][j]+4);
		}
	}

	xt::xarray<double> theoretical_lambda = xt::zeros<double>({T,R,C});
	vector<int> type_region = vector<int>(R,-1);
	int count = 0;
	for(int i = 0; i < n_y/2; ++i){
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 0;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.1;
			}
			++count;
		}

		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 1;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.5;
			}
			++count;
		}
	}

	for(int i = n_y/2; i < n_y; ++i){
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 1;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.5;
			}
			++count;
		}
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 0;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.1;
			}
			++count;
		}
	}

	int nb_weeks = 10;
	int nb_years = 1;
	int nb_observations_total = nb_weeks*nb_years;
	vector<double> durations = vector<double>(T,6);
	std::random_device rd;
    std::mt19937 gen(rd());


	int max_obs = nb_observations_total;
	xt::xarray<int> sample = xt::zeros<int>({T,R,C,static_cast<ulong>(max_obs)});
	xt::xarray<int> nb_observations = xt::zeros<int>({C,R,T});
	xt::xarray<int> nb_arrivals = xt::zeros<int>({C,R,T});
	for(int t = 0; t < T; ++t){
		for(int r = 0; r < R; ++r){
			for(int c = 0; c < C; ++c){
				nb_observations(c,r,t) = nb_observations_total;
				double sum = 0;
				poisson_distribution<int> pd(theoretical_lambda(t,r,c)*durations[t]);
				for(int k = 0; k < nb_observations_total; ++k){
					sample(t,r,c,k) = pd(gen);
					sum += sample(t,r,c,k);
					nb_arrivals(c,r,t) += sample(t,r,c,k);
				}
			}
		}
	}

	vector<vector<int>> neighbors = vector<vector<int>>(R, vector<int>());
	xt::xarray<double> distance = GRB_INFINITY*xt::ones<double>({R, R});
	for(int r = 0; r < R; ++r){
		int Xi = (r+1) % n_x;
		int Yi = -1;
		if (Xi==0){
			Yi= (r+1)/n_x;
		}else{
			Yi=1+(r+1-Xi)/n_x;
		}
		if (Xi==0){
			Xi=n_x;
		}

		if (Yi-1>=1){
			neighbors[r].push_back((Yi-2)*n_x+Xi - 1);     
		}

		if (Yi+1<=n_y){
			neighbors[r].push_back(Yi*n_x+Xi - 1);     
		}

		if (Xi-1>=1){
			neighbors[r].push_back((Yi-1)*n_x+Xi-1 - 1);     
		}

		if (Xi+1<=n_x){
			neighbors[r].push_back((Yi-1)*n_x+Xi+1 - 1);
		}

		if ((Yi+1<=n_y)&&(Xi-1>=1)){
			neighbors[r].push_back(Yi*n_x+Xi-1 - 1);   
		}

		if ((Yi-1>=1)&&(Xi-1)>=1){
			neighbors[r].push_back((Yi-2)*n_x+Xi-1 - 1);    
		}

		if ((Yi+1<=n_y)&&(Xi+1)<=n_x){
			neighbors[r].push_back(Yi*n_x+Xi+1 - 1);     
		}

		if ((Yi-1>=1)&&(Xi+1)<=n_x){
			neighbors[r].push_back((Yi-2)*n_x+Xi+1 - 1);
		}

	}

	vector<double> absc(n_x, 0);
	for(int i = 0; i < n_x; ++i){
		absc[i] = x_max/(2*n_x) + (x_max / n_x)*i;
	}
	vector<double> ord(n_x, 0);
	for(int j = 0; j < n_y; ++j){
		ord[j] = y_max / (2*n_y) + (y_max / n_y)*j;
	}

	for(int r = 0; r < R; ++r){
		for(int s: neighbors[r]){
			distance(r,s) = 1;
		}
	}

	xt::xarray<double> alphas = xt::ones<double>({R,R});
	vector<double> weights = vector<double>(nb_groups, 1);

	xt::xarray<double> lambda0 = xt::ones<double>({C,R,T});
	xt::xarray<double> lambda = laspated_no_reg(nb_observations, nb_arrivals, durations, groups,
		weights, alphas, distance, type_region, neighbors, lambda0);

	ofstream arq("x_no_reg.txt", std::ios::out);
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				arq << c << " " << r << " " << t << " " << lambda(c,r,t) << "\n";
			}
		}
	}
	arq.close();
	fmt::print("Lambdas saved at x_no_reg.txt\n");
}


int main(int argc, char* argv[]){
	srand(time(NULL));
	using namespace std;
	g_params = load_params(argc,argv);

	test_laspated_reg();
	// test_laspated_no_reg();
	// std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);

	// bool is_custom_data = g_params.generator_folder != "";
	// if(g_params.method == ""){
	// 	if(g_params.model == "no_reg"){
	// 		if(is_custom_data){
	// 			std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// 			std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// 			std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);
	// 			GeneratorNoRegressor gen(calls_path, neighbors_path, info_path);
	// 			auto t0 = std::chrono::high_resolution_clock::now();
	// 			gen.calibrate();
	// 			auto dt = std::chrono::high_resolution_clock::now();
	// 			double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 				/ pow(10,9);
	// 			fmt::print("Run time = {}\n", run_time_old);
	// 		}else{
	// 			GeneratorNoRegressor gen;
	// 			auto t0 = std::chrono::high_resolution_clock::now();
	// 			gen.test();
	// 			auto dt = std::chrono::high_resolution_clock::now();
	// 			double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 				/ pow(10,9);
	// 			fmt::print("Run time = {}\n", run_time_old);
	// 		}
	// 	}else if(g_params.model == "reg"){
	// 		GRBEnv env;
	// 		if(is_custom_data){
	// 			std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// 			std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// 			std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);
	// 			GeneratorRegressor gen(env, calls_path, neighbors_path, info_path);
	// 			auto t0 = std::chrono::high_resolution_clock::now();
	// 			gen.calibrate();
	// 			auto dt = std::chrono::high_resolution_clock::now();
	// 			double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 				/ pow(10,9);
	// 			fmt::print("Run time = {}\n", run_time_old);
	// 		}else{
	// 			GeneratorRegressor gen(env);
	// 			auto t0 = std::chrono::high_resolution_clock::now();
	// 			gen.calibrate();
	// 			auto dt = std::chrono::high_resolution_clock::now();
	// 			double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 				/ pow(10,9);
	// 			fmt::print("Run time = {}\n", run_time_old);
	// 		}
	// 	}
	// }

	// if(g_params.model == "no_reg" && g_params.method == "calibration"){
	// 	if(is_custom_data){
	// 		std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// 		std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// 		std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);
	// 		GeneratorNoRegressor gen(calls_path, neighbors_path, info_path);
	// 		auto t0 = std::chrono::high_resolution_clock::now();
	// 		gen.calibrate();
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}else{
	// 		GeneratorNoRegressor gen;
	// 		auto t0 = std::chrono::high_resolution_clock::now();
	// 		gen.calibrate();
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}
	// }else if(g_params.model == "no_reg" && g_params.method == "cross_validation"){
	// 	if(is_custom_data){
	// 		std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// 		std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// 		std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);
	// 		GeneratorNoRegressor gen(calls_path, neighbors_path, info_path);
	// 		auto t0 = std::chrono::high_resolution_clock::now();

	// 		double proportion = g_params.cv_proportion;
	// 		auto weights = g_params.weights_list;
	// 		auto alpha = weights;
	// 		auto result = gen.cross_validation(proportion, alpha, weights);
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		gen.write_cv_results(result);
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}else{
	// 		GeneratorNoRegressor gen;
	// 		double proportion = g_params.cv_proportion;
	// 		vector<double> test_weights = g_params.weights_list;
	// 		vector<double> alphas = test_weights;
	// 		auto t0 = std::chrono::high_resolution_clock::now();
	// 		auto result = gen.cross_validation(proportion, alphas, test_weights);
	// 		gen.write_cv_results(result);
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}
	// }else if(g_params.model == "reg"){
	// 	GRBEnv env;
	// 	if(is_custom_data){
	// 		std::string calls_path = fmt::format("{}/calls.dat", g_params.generator_folder);
	// 		std::string neighbors_path = fmt::format("{}/neighbors.dat", g_params.generator_folder);
	// 		std::string info_path = fmt::format("{}/info.dat", g_params.generator_folder);
	// 		GeneratorRegressor gen(env, calls_path, neighbors_path, info_path);
	// 		auto t0 = std::chrono::high_resolution_clock::now();
	// 		gen.calibrate();
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}else{
	// 		GeneratorRegressor gen(env);
	// 		auto t0 = std::chrono::high_resolution_clock::now();
	// 		gen.calibrate();
	// 		auto dt = std::chrono::high_resolution_clock::now();
	// 		double run_time_old = std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() 
	// 			/ pow(10,9);
	// 		fmt::print("Run time = {}\n", run_time_old);
	// 	}
	// }

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
	    ("generator_folder,G", po::value<std::string>()->default_value(""),
	    	"Directory with observations data. Must have a calls.dat, info.dat and neighbors.dat files. default=\"\" (use of simulated data)")
		("method", po::value<std::string>()->default_value("calibration"),
			"Type of function. calibration|cross_validation")
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
		("cv_proportion", po::value<double>()->default_value(0.2),
			"Proportion of data allocated to training set. Real value between 0 and 1. default = 0.2")
		("weights_file", po::value<std::string>()->default_value(""),
			"path for file containing weights to be used. default = \"\"")
		("weights_list", po::value<std::vector<double>>()->multitoken(),
			"list of weights to be used.")
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

	try{
		auto params = Param(vm);
	}catch(boost::wrapexcept<boost::bad_any_cast>& e){
		std::cerr << "Erro: " << e.what() << "\n";
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