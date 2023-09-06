#include "../include/generator_no_regressor.h"

using namespace std;

//Hardcoded test
GeneratorNoRegressor::GeneratorNoRegressor(){
	x_max = y_max = 10;
	n_x = n_y = 10;
	R = n_x*n_y;
	C = 1;
	int nb_groups = 2;
	T = 4*7;

	which_group = vector<int>(T, -1);
	groups = vector<vector<int>>(nb_groups,vector<int>());
	for(int t = 0; t < T; ++t){
		groups[t % nb_groups].push_back(t);
		which_group[t] = t % nb_groups;
	}

	// for(int i = 0; i < nb_groups; ++i){
	// 	fmt::print("Groups {}: ", i+1);
	// 	for(auto j: groups[i]){
	// 		fmt::print("{} ", j+1);
	// 	}
	// 	fmt::print("\n");
	// }
	// cin.get();

	type_region = vector<int>(R,-1);
	int count = 0;
	for(int i = 0; i < n_y/2; ++i){
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 0;
			++count;
		}

		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 1;
			++count;
		}
	}

	for(int i = n_y/2; i < n_y; ++i){
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 1;
			++count;
		}
		for(int j = 0; j < n_x/2; ++j){
			type_region[count] = 0;
			++count;
		}
	}



	nb_weeks = 1000;
	nb_years = 1;
	int nb_observations_total = nb_weeks*nb_years;
	int max_obs = nb_observations_total;
	durations = vector<double>(T,1);
	std::default_random_engine gen(600);
	// std::random_device rd;
    // std::mt19937 gen(rd());


	sample = xt::zeros<int>({T,R,C,static_cast<ulong>(max_obs)});
	nb_observations = nb_observations_total*xt::ones<int>({C,R,T});
	nb_arrivals = xt::zeros<int>({C,R,T});
	uniform_real_distribution<double> rand_coord(0,1);
	for(int t = 0; t < T; ++t){
		if(t % 2  != 0){
			double lambda_bar_b = 75;
			double lambda_b = 625;
			double lambda_bar_r = 20;
			double lambda_r = 250;

			for(int n = 0; n < nb_observations_total; ++n){
				poisson_distribution<int> pd_Nb(lambda_b);
				int N = pd_Nb(gen);
				for(int k = 0; k < N; ++k){
					bool contde = true;
					while(contde){
						bool contd = true;
						int x = 0,y = 0;
						while(contd){
							x = 9*rand_coord(gen);
							y = 9*rand_coord(gen);

							if((x < 5 && y >= 5 && y < 10) || (x >= 5 && x < 10 && y < 5)){
								contd = false;
							}
						}
						double us = lambda_bar_b*rand_coord(gen);
						if(us <= 5*(x+y)){
							contde = 0;
							int ind_x = floor(x);
							int ind_y = floor(y);
							int zone_nb = ind_y*10+ind_x;
							sample(t,zone_nb,0,n) += 1;
						}
					}
				}
				poisson_distribution<int> pd_Nr(lambda_r);
				N = pd_Nr(gen);
				for(int k = 0; k < N; ++k){
					bool contde = true;
					while(contde){
						bool contd = true;
						int x = 0,y = 0;
						while(contd){
							x = 9*rand_coord(gen);
							y = 9*rand_coord(gen);

							if((x < 5 && y < 5) || (x >= 5 && x < 10 && y >= 5 && y < 10)){
								contd = false;
							}
						}
						double us = lambda_bar_r*rand_coord(gen);
						if(us <= x+y){
							contde = 0;
							int ind_x = floor(x);
							int ind_y = floor(y);
							int zone_nr = ind_y*10+ind_x;
							sample(t,zone_nr,0,n) += 1;
						}
					}
				}
			}
		}else{
			double lambda_bar_b = 15;
			double lambda_b = 125;
			double lambda_bar_r = 100;
			double lambda_r = 1250;

			for(int n = 0; n < nb_observations_total; ++n){
				poisson_distribution<int> pd_Nb(lambda_b);
				int N = pd_Nb(gen);
				for(int k = 0; k < N; ++k){
					bool contde = true;
					while(contde){
						bool contd = true;
						int x = 0,y = 0;
						while(contd){
							x = 9*rand_coord(gen);
							y = 9*rand_coord(gen);

							if((x < 5 && y >= 5 && y < 10) || (x >= 5 && x < 10 && y < 5)){
								contd = false;
							}
						}
						double us = lambda_bar_b*rand_coord(gen);
						if(us <= x+y){
							contde = 0;
							int ind_x = floor(x);
							int ind_y = floor(y);
							int zone_nb = ind_y*10+ind_x;
							sample(t,zone_nb,0,n) += 1;
						}
					}
				}
				poisson_distribution<int> pd_Nr(lambda_r);
				N = pd_Nr(gen);
				for(int k = 0; k < N; ++k){
					bool contde = true;
					while(contde){
						bool contd = true;
						int x = 0,y = 0;
						while(contd){
							x = 9*rand_coord(gen);
							y = 9*rand_coord(gen);

							if((x < 5 && y < 5) || (x >= 5 && x < 10 && y >= 5 && y < 10)){
								contd = false;
							}
						}
						double us = lambda_bar_b*rand_coord(gen);
						if(us <= 5*(x+y)){
							contde = 0;
							int ind_x = floor(x);
							int ind_y = floor(y);
							int zone_nb = ind_y*10+ind_x;
							sample(t,zone_nb,0,n) += 1;
						}
					}
				}
			}
		}
	}


	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				for(int n = 0; n < nb_observations_total; ++n){
					nb_arrivals(c,r,t) += sample(t,r,c,n);
				}
			}
		}
	}


	vector<bool> is_blue(100,false);
	int u = 0;
	for(int k = 0; k < 5; ++k){
		u += 5;
		for(int i = 0; i < 5; ++i){
			is_blue[u+i] = true;
		}
	}

	for(int k = 0; k < 5; ++k){
		for(int i = 0; i < 5; ++i){
			is_blue[u+i] = true;
		}
		u += 10;
	}

	theoretical_lambda = xt::zeros<double>({T,R,C});
	xt::xarray<double> integrated_lambda = xt::zeros<double>({T,R});

	for(int r = 1; r <= R; ++r){
		int i = r % 10;
		int j = floor((r-1) / 10) + 1;

		if(i == 0){
			i = 10;
		}

		for(int t = 0; t < T; ++t){
			if(t % 2 != 0){
				if(is_blue[r-1]){
					theoretical_lambda(t,r-1, 0) = 5*(i+j-1);
					integrated_lambda(t,r-1) = 2.5*(pow(i,2)-pow(i-1,2)) + 2.5*(pow(j,2) - pow(j-1,2));
				}else{
					theoretical_lambda(t,r-1, 0) = i+j-1;
					integrated_lambda(t,r-1) = 0.5*(pow(i,2)-pow(i-1,2)) + 0.5*(pow(j,2) - pow(j-1,2));
				}
			}else{
				if(is_blue[r-1]){
					theoretical_lambda(t,r-1, 0) = i+j-1;
					integrated_lambda(t,r-1) = 0.5*(pow(i,2)-pow(i-1,2)) + 0.5*(pow(j,2) - pow(j-1,2));
				}else{
					theoretical_lambda(t,r-1, 0) = 5*(i+j-1);
					integrated_lambda(t,r-1) = 2.5*(pow(i,2)-pow(i-1,2)) + 2.5*(pow(j,2) - pow(j-1,2));
				}
			}
		}
	}

	xt::xarray<double> estimated = xt::zeros<double>({C,R,T});
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				estimated(c,r,t) = nb_arrivals(c,r,t) / (nb_observations(c,r,t)*durations[t]);
			}
		}
	}


	xt::xarray<double> est = xt::zeros<double>({2,4});
	for(int t = 0; t < 4; ++t){
		int nb_type0 = 0;
		int nb_type1 = 0;
		for(int r = 0; r < R; ++r){
			if(type_region[r] == 0){
				est(0,t) += nb_arrivals(0,r,t);
				nb_type0 += nb_observations_total;
			}else{
				est(1,t) += nb_arrivals(0,r,t);
				nb_type1 += nb_observations_total;
			}
		}
		// fmt::print("t = {}, est1 = {}, est2 = {}\n",t+1, est(0,t), est(1,t));
		est(0,t) =  est(0,t) / (nb_type0 * durations[t]);
		est(1,t) = 	est(1,t) / (nb_type1 * durations[t]);
		// fmt::print("t = {}, est1 = {}, est2 = {}\n",t+1, est(0,t), est(1,t));
	}



	neighbors = vector<vector<int>>(R, vector<int>());
	distance = GRB_INFINITY*xt::ones<double>({R, R});
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
		int x_r = r % n_x;
		int y_r = r / n_y;

		for(int s: neighbors[r]){
			int x_s = s % n_x;
			int y_s = s / n_y;

			// distance(r,s) = sqrt(pow(x_r - x_s, 2) + pow(y_r - y_s, 2));
			// if(x_r == x_s || y_r == y_s){
			// 	distance(r,s) = 1;
			// }else{
			// 	distance(r,s) = pow(2,0.5);
			// }
			// distance(r,s) = sqrt(pow(absc[x_r] - absc[x_s], 2) + pow(ord[y_r]-ord[y_s], 2));
			distance(r,s) = 1;
		}
		// fmt::print("Neighbors {} (type_region = {}): ", r+1, type_region[r]+1);
		// for(auto s: neighbors[r]){
		// 	fmt::print("{} ", s+1);
		// }
		// fmt::print("\n");
	}
	// cin.get();


	l_bounds = xt::zeros<double>({C,R,T});
	g_params.EPS = pow(10,-3);
	alpha = 1;
	weight = 1;
	weights = vector<double>(nb_groups, 1);
	sigma = 0.5;
	beta_bar = 1;
	max_iter = 30;
	std::cout << "Initialized No Regressor\n";
}



GeneratorNoRegressor::GeneratorNoRegressor(std::string calls_path, 
		std::string neighbors_path, std::string info_path){
	auto info_arq = ifstream(info_path, ios::in);
	info_arq >> T >> D >> R >> C >> nb_regressors >> nb_holidays_years;
	slot_duration = 24 / T;
	fmt::print("info: {} {} {} {} {} {}\n", T,D,R,C,nb_regressors, nb_holidays_years);
	daily_obs = std::vector<int>(D, 0);
	for(int d = 0; d < D; ++d){
		info_arq >> daily_obs[d];
	}
	fmt::print("daily obs: {}\n", daily_obs);
	info_arq.close();
	nb_land_types = nb_regressors - 2;
	unsigned long max_obs = *max_element(daily_obs.begin(), daily_obs.end());
	unsigned long min_obs = *min_element(daily_obs.begin(), daily_obs.end());

	xt::xarray<int> nb_observations_file = xt::zeros<int>({C,D,T,R});
	xt::xarray<int> sample_file = xt::zeros<int>({C,D,T,R, max_obs});
	xt::xarray<int> nb_arrivals_file = xt::zeros<int>({C,D,T,R});

	auto calls_arq = ifstream(calls_path, ios::in);
	std::string aux_str;
	// is_holidays = std::vector<std::pair<bool, int>>(max_obs, make_pair(false,-1));
	do{
		std::getline(calls_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		int t, d, r, c, j, h, val;
		ss >> t >> d >> r >> c >> j >> val >> h;
		// fmt::print("calls {} {} {} {} {} {} {}\n",t,d,r,c,j,val,h);
		sample_file(c,d,t,r,j) = val;
		nb_arrivals_file(c,d,t,r) += val;
		nb_observations_file(c,d,t,r) += 1;
	}while(true);
	calls_arq.close();

	xt::xarray<double> estimated = xt::zeros<double>({C,D,T,R});
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					nb_observations(c,d,t,r) = daily_obs[d];
					estimated(c,d,t,r) = nb_arrivals(c,d,t,r) / nb_observations(c,d,t,r);
				}
			}
		}
	}


	type_region = std::vector<int>(R,-1);
	xt::xarray<double> regressors = xt::zeros<double>({nb_regressors, R});
	neighbors = std::vector<vector<int>>(R, std::vector<int>());
	distance = xt::zeros<double>({R, R});
	regions = std::vector<Location>(R, null_location);
	auto neighbors_arq = ifstream(neighbors_path, ios::in);
	double pop1,pop2;
	while(true){
		int ind, terrain_type, s; 
		double lat, longi, dist;
		std::getline(neighbors_arq, aux_str);
		if(aux_str == "END"){
			break;
		}
		std::istringstream ss(aux_str);
		ss >> ind >> lat >> longi >> terrain_type;
		type_region[ind] = terrain_type;
		regions[ind] = make_pair(lat, longi);
		for(int j = 0; j < nb_land_types; ++j){
			ss >> regressors(j,ind);
		}
		ss >> pop1 >> pop2;
		regressors(nb_regressors - 2, ind) = pop1;
		regressors(nb_regressors - 1, ind) = pop1 + pop2;
		while(ss >> s >> dist){
			distance(ind,s) = dist;
			neighbors[ind].push_back(s);
		}
	}
	neighbors_arq.close();

	sample = xt::zeros<int>({7*T,R,C,min_obs});
	nb_observations = xt::zeros<int>({C,R,7*T});
	nb_arrivals = xt::zeros<int>({C,R,7*T});
	durations = vector<double>(7*T,0.5);

	for(int r = 0; r < R; ++r){
		for(int c = 0; c < C; ++c){
			for(int d = 0; d < D; ++d){
				for(int t = 0; t < T; ++t){
					int index = d*T + t;
					nb_observations(c,r,index) = nb_observations_file(c,d,t,r);
					nb_arrivals(c,r,index) = nb_arrivals_file(c,d,t,r);
					for(int j = 0; j < min_obs; ++j){
						sample(index,r,c,j) = sample_file(c,d,t,r,j);
					}
				}
			}
		}
	}

	T = 7*T;
	which_group = vector<int>(T, 0);
	groups = vector<vector<int>>(T,vector<int>());
	for(int t = 0; t < T; ++t){
		groups[t].push_back(t);
		which_group[t] = t;
		// fmt::print("which_group {} = {}\n", t+1, which_group[t]+1);
	}
	// cin.get();
	sigma = 0.5;
	beta_bar = 1;
	max_iter = 30;
	weight = 0.1;
	alpha = 0.1*xt::ones<double>({R,R});
	g_params.EPS = 0.001;
	durations = vector<double>(T,0.5);
	std::cout << "Initialized No Regressor Real Data\n";
}



GeneratorNoRegressor::GeneratorNoRegressor(xt::xarray<int>& N, xt::xarray<int>& M, std::vector<double>& a_durations, 
		std::vector<std::vector<int>>& a_groups, std::vector<double>& a_weights, 
		xt::xarray<double>& alphas, xt::xarray<double>& a_distance, std::vector<int>& a_type, 
		std::vector<std::vector<int>>& a_neighbors){
	if(N.dimension() !=  3){ //N should be C,R,T
		fmt::print("Error: N has {} dimensions but should be 3. Problem was not set.\n",N.dimension());
		exit(1);
	}
	if(M.dimension() != 3){ //M also should be C,R,T
		fmt::print("Error: M has {} dimensions but should be 3. Problem was not set.\n",M.dimension());
		exit(1);
	}

	C = N.shape(0);
	R = N.shape(1);
	T = N.shape(2);

	nb_observations = N;
	nb_arrivals = M;
	durations  = a_durations;
	alpha = alphas;
	weights = a_weights;
	distance = a_distance;
	neighbors = a_neighbors;
	type_region = a_type;
	groups = a_groups;
	int nb_groups = groups.size();

	which_group = vector<int>(T, 0);
	for(int g = 0; g < nb_groups; ++g){
		for(auto i: groups[g]){
			which_group[i] = g;
		}
	}

	g_params.EPS = pow(10,-3);
	sigma = 0.5;
	beta_bar = 1;
	max_iter = 30;
}


xt::xarray<double> laspated_no_reg(xt::xarray<int>& N, xt::xarray<int>& M, std::vector<double>& a_durations, 
		std::vector<std::vector<int>>& a_groups, std::vector<double>& a_weights, 
		xt::xarray<double>& alphas, xt::xarray<double>& a_distance, std::vector<int>& a_type, 
		std::vector<std::vector<int>>& a_neighbors, xt::xarray<double>& x){
	
	GeneratorNoRegressor gen(N,M,a_durations, a_groups, a_weights, alphas,  a_distance, a_type, a_neighbors);

	if(x.dimension() != 3){
		fmt::print("Error: x has {} dimensions but must be 3.\n", x.dimension());
		exit(1);
	}

	if(x.shape(0) != gen.C || x.shape(1) != gen.R || x.shape(2) != gen.T){
		std::string shape_x = fmt::format("({},{},{})", x.shape(0), x.shape(1), x.shape(2));
		fmt::print("Error: x has shape {}, but expected is ({},{},{}).\n", shape_x,gen.C, gen.R, gen.T);
		exit(1);
	}

	auto lambda = x;
	auto f_val = gen.projected_gradient_armijo_feasible(lambda);
	return lambda;

}



void GeneratorNoRegressor::test(){
	double epsilon = g_params.EPS;
	// vector<double> test_weights = {0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
	// 	110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
	// 	270,280,290,300,310,320,330,340,350,360,370,380,390,400};
	vector<double> test_weights = g_params.weights_list;
	vector<double> test_alphas = test_weights;
	xt::xarray<double> x = xt::ones<double>({C,R,T});
	ofstream err_arq("err_no_reg.txt", std::ios::out);
	ofstream week_arq(fmt::format("week_no_reg.txt"), std::ios::out);
	double min_err = 10e100;
	double min_w = -1;
	for(int i = 0; i < test_weights.size(); ++i){
		x = epsilon*xt::ones<double>({C,R,T});
		alpha = test_alphas[i]*xt::ones<double>({R,R});
		// alpha = 0;
		weights = vector<double>(groups.size(), test_weights[i]);
		// weights = vector<double>(groups.size(), 0);
		auto f_val = projected_gradient_armijo_feasible(x);

		if(g_params.generator_folder == ""){
			double err = average_difference(x);
			fmt::print("weight {} err {}\n", test_weights[i], err);
			if(err < min_err){
				min_err = err;
				min_w = i;
			}
			err_arq << fmt::format("{:.7f}",err) << "\n";
		}
	}
	if(g_params.generator_folder == ""){
		fmt::print("Wrote weight results at err_no_reg.txt\n");
	}

	test_alphas = test_weights;
	auto result = cross_validation(0.2, test_alphas, test_weights);
	write_cv_results(result);
	fmt::print("Cross validation time = {}\n",result.cpu_time);
	fmt::print("Cross validation weight = {}\n",result.weight);
	double best_w = result.weight;
	x = result.lambda;
	auto f_val = projected_gradient_armijo_feasible(x);
	// fmt::print("best_weight err = {}\n", average_difference(x));
	// err_arq << fmt::format("{:.3f}\n{:.7f}\n", best_w, average_difference(x));
	// err_arq.close();
	ofstream x_arq(fmt::format("x_no_reg_d{}.txt", durations[0]), std::ios::out);
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				x_arq << fmt::format("{} {} {} = {}\n",c,r,t, result.lambda(c,r,t));
			}
		}
	}
	x_arq.close();
	fmt::print("Wrote cross_validation intensities at {}\n", fmt::format("x_no_reg_d{}.txt", durations[0]));
}

void GeneratorNoRegressor::calibrate(){
	double epsilon = g_params.EPS;
	vector<double> test_weights = g_params.weights_list;
	vector<double> alphas = test_weights;
	xt::xarray<double> x = xt::ones<double>({C,R,T});
	fmt::print("Running projected gradient for weights {}\n", test_weights);
	ofstream err_arq(fmt::format("err_no_diag_g{}_obs{}_alpha.txt", groups.size(), nb_weeks), std::ios::out);
	ofstream week_arq(fmt::format("week_no_reg.txt"), std::ios::out);
	double min_err = 10e100;
	int min_w = -1;
	for(int i = 0; i < test_weights.size(); ++i){
		x = epsilon*xt::ones<double>({C,R,T});
		alpha = alphas[i]*xt::ones<double>({R,R});
		// alpha = xt::zeros<double>({R,R});
		// alpha = 0;
		weights = vector<double>(groups.size(), test_weights[i]);
		// weights = vector<double>(groups.size(), 0);
		auto f_val = projected_gradient_armijo_feasible(x);
		// fmt::print("alpha = {}, weight = {}, diff = {}\n",alpha, weights[0], average_difference(x));
		
		if(g_params.generator_folder == ""){
			double err = average_difference(x);
			if(err < min_err){
				min_err = err;
				min_w = i;
			}
			err_arq << fmt::format("{:.7f}",err) << "\n";
			fmt::print("w({}) = {}, err = {}\n", i, test_weights[i], err);
		}
	}
	alpha = alphas[min_w];
	if(g_params.generator_folder == ""){
		fmt::print("Best weight = {} {}\n", alpha, min_w);
		fmt::print("Wrote weight results at err_no_reg.txt\n");
	}
	
	weights = vector<double>(groups.size(), test_weights[min_w]);
	x = epsilon*xt::ones<double>({C,R,T});
	auto f_val = projected_gradient_armijo_feasible(x);

	ofstream x_arq(fmt::format("x_no_reg.txt"), std::ios::out);
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				x_arq << fmt::format("{} {} {} = {}\n",c,r,t, x(c,r,t));
			}
		}
	}
	x_arq.close();
	fmt::print("Wrote intensities at x_no_reg.txt\n");
}


void GeneratorNoRegressor::write_cv_results(CrossValidationResult& cv_result){
	ofstream arq("cv_x_no_reg.txt", std::ios::out);
	arq << "Best weight = "<< cv_result.weight << ", cpu time (s) = " << cv_result.cpu_time << "\n";
	xt::xarray<double> x = g_params.EPS*xt::ones<double>(cv_result.lambda.shape());
	auto f_val = projected_gradient_armijo_feasible(x);
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				arq << c << " " << r << " " << t << " " << x(c,r,t) << "\n";
			}
		}
	}
	arq.close();
	fmt::print("Wrote intensities at cv_x_no_reg.txt\n");
	ofstream weekly_total("weekly_total_no_reg.txt", std::ios::out);
	for(int t = 0; t < T; ++t){
		double sum = 0;
		for(int c = 0; c < C; ++c){
			for(int r = 0; r < R; ++r){
				sum += x(c,r,t);
			}
		}
		weekly_total << sum << "\n";
	}
	weekly_total.close();

	ofstream weekly0("weekly_0_no_reg.txt", std::ios::out);
	for(int t = 0; t < T; ++t){
		double sum = 0;
		for(int r = 0; r < R; ++r){
			sum += x(0,r,t);
		}
		weekly0 << sum << "\n";
	}
	weekly0.close();

	ofstream weekly1("weekly_1_no_reg.txt", std::ios::out);
	for(int t = 0; t < T; ++t){
		double sum = 0;
		for(int r = 0; r < R; ++r){
			sum += x(1,r,t);
		}
		weekly1 << sum << "\n";
	}
	weekly1.close();

	ofstream weekly2("weekly_2_no_reg.txt", std::ios::out);
	for(int t = 0; t < T; ++t){
		double sum = 0;
		for(int r = 0; r < R; ++r){
			sum += x(1,r,t);
		}
		weekly2 << sum << "\n";
	}
	weekly2.close();

}

std::vector<double> GeneratorNoRegressor::projected_gradient_armijo_feasible(
	xt::xarray<double>& x){
	double eps = 0.001;
	xt::xarray<double> z = xt::zeros<double>(x.shape());
	// fmt::print("alpha = {}\n", alpha);
	int k = 0;
	std::vector<double> f_val;
	double b_param = 2;
	double beta_k = b_param;
	f_val.reserve(max_iter);
	int j = 0;
	
	while(k < max_iter){
		double fold = oracle_objective_model(x);
		xt::xarray<double> gradient = oracle_gradient_model(x);
		xt::xarray<double> z = x - beta_k*gradient;
		comp_wise_max(z,eps);

		bool stop = false;
		j = 0;
		xt::xarray<double> diff_aux = x-z;
		double rhs = mat_prod(gradient, diff_aux);
		// fmt::print("fold = {:.6f} rhs = {:.6f}\n", fold, rhs);
		double f = GRB_INFINITY;
		xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
		while(!stop){
			z_aux = x + (1/pow(2,j))*(z - x);
			f = oracle_objective_model(z_aux);
			// fmt::print("\tj = {} f = {:.6f} test = {:.6f}\n", j, f, fold -(sigma/pow(2,j))*rhs);
			if(f <= fold - (sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}
		f_val.push_back(f);
		x = z_aux;
		// print_var(x,fmt::format("final x iter {}", k));
		// fmt::print("k = {}\n", k+1);
		beta_k = b_param / pow(2,j);
		++k;
		// cin.get();
	}
	return f_val;
}

double GeneratorNoRegressor::mat_prod(xt::xarray<double>& a, xt::xarray<double>& b){
	double sum = 0;
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				sum += a(c,r,t)*b(c,r,t);
			}
		}
	}
	return sum;
}



xt::xarray<double> GeneratorNoRegressor::oracle_gradient_model(xt::xarray<double>& x){
	
	xt::xarray<double> gradient = xt::zeros<double>(x.shape());

	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				double current_lambda = x(c,r,t);
				double grad_component = nb_observations(c,r,t)*durations[t] - 
					(nb_arrivals(c,r,t) / current_lambda);
				// double prev_comp = grad_component;
				// double sum_neighbors = 0;
				for(int s: neighbors[r]){
					if(type_region[r] == type_region[s])
					// if(true)
					{
						grad_component += 2*alpha(r,s)*(x(c,r,t) - x(c,s,t)) / (distance(r,s));
					}
				}
				gradient(c,r,t) = grad_component;
			}
		}
	}
	
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				double sum_groups = 0;
				auto& group = groups[which_group[t]];
				for(int j = 0; j < group.size(); ++j){
					int tp = group[j];
					if(tp != t){
						gradient(c,r,t) += 2*weights[which_group[t]]*(x(c,r,t)-x(c,r,tp));
					}
				}
			}
		}
	}
	
	return gradient;
}


double GeneratorNoRegressor::oracle_objective_model(xt::xarray<double>& x){
	double f = 0;
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				double current_lambda = x(c,r,t);
				f += nb_observations(c,r,t)*current_lambda*durations[t] - 
					nb_arrivals(c,r,t)*log(current_lambda*durations[t]);
				for(int s: neighbors[r]){
					if(type_region[r] == type_region[s])
					// if(true)
					{
						f += (0.5*alpha(r,s))*pow(x(c,r,t)- x(c,s,t), 2) / distance(r,s);
					}
					// f += (0.5*alpha(r,s))*pow(x(c,r,t)- x(c,s,t), 2) / distance(r,s);
				}
			}
		}
	}
	//sum_{c \in C}sum{r \in R}sum_{t \in T}\sum_{t' \neq t \in which_group[t]}
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int grindex = 0; grindex < groups.size(); ++grindex){
				auto& group = groups[grindex];
				for(int j = 0; j < group.size(); ++j){
					int t = group[j];
					for(int itp = 0; itp < group.size(); ++itp){
						int tp = group[itp];
						if(tp != t){
							f += (0.5*weights[grindex])*(pow(x(c,r,t) - 
								x(c,r,tp), 2));
						}
					}
				}
			}
		}
	}

	return f;
}


double GeneratorNoRegressor::average_difference(xt::xarray<double>& x){
	vector<double> difference_l2;
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				difference_l2.push_back(abs(theoretical_lambda(t,r,c)-x(c,r,t)) / theoretical_lambda(t,r,c));
			}
		}
	}
	// fmt::print("Error: {}\n", difference_l2);
	return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) / difference_l2.size();
}



void GeneratorNoRegressor::comp_wise_max(xt::xarray<double>& z, double eps){
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				z(c,r,t) = max(z(c,r,t), eps);
			}
		}
	}
}


bool GeneratorNoRegressor::is_neighbor(int r, int s){
	return r != s;
}

void GeneratorNoRegressor::print_var(xt::xarray<double>& x,
	std::string prefix){
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				if(abs(x(c,r,t)) > g_params.EPS){
					fmt::print("nonzero: c = {}, r = {}, t = {} (type_region = {}), val  = {}, theoretical = {}\n",c+1,r+1,t+1,
						type_region[r]+1, x(c,r,t), theoretical_lambda(t,r,c));
				}
			}
		}
	}
	// fmt::print("{} c = 0, r = 84, t = 27, val = {}\n", prefix, x(0,84,27));
	// fmt::print("{} c = 0, r = 85, t = 27, val = {}\n", prefix, x(0,85,27));
	cin.get();
}


CrossValidationResult GeneratorNoRegressor::cross_validation(double proportion, vector<double>& alphas, 
		vector<double>& group_weights){
	auto t0 = std::chrono::high_resolution_clock::now();
	int nb_observations_total = sample.shape(3);
	int nb_groups = groups.size();
	double min_loss = GRB_INFINITY;
	double cpu_time = 0;
	int nb_in_block = floor(nb_observations_total * proportion);
	xt::xarray<int> initial_nb_obs = nb_observations;
	xt::xarray<int> initial_nb_arrivals = nb_arrivals;
	fmt::print("cross validation weights.size = {}\n", group_weights.size());
	double best_alpha = GRB_INFINITY;
	double best_weight = GRB_INFINITY;
	fmt::print("Running cross validation with proportion = {} and weights = {}\n", proportion, group_weights);
	for(int index_alpha = 0; index_alpha < alphas.size(); ++index_alpha){
		double likelihood = 0;
		alpha = alphas[index_alpha]*xt::ones<double>({R,R});
		weights = vector<double>(groups.size(), group_weights[index_alpha]);
		fmt::print("Testing weight = {}\n", group_weights[index_alpha]);
		for(int index_cross = 0; index_cross < floor(1/proportion); ++index_cross){
			xt::xarray<int> nb_observations_current = xt::zeros<int>({C,R,T});
			xt::xarray<int> nb_calls_current = xt::zeros<int>({C,R,T});
			for(int index = index_cross*nb_in_block; index < (index_cross+1)*nb_in_block; ++index){
				for(int c = 0; c < C; ++c){
					for(int r = 0; r < R; ++r){
						for(int t = 0; t < T; ++t){
							++nb_observations_current(c,r,t);
							nb_calls_current(c,r,t) += sample(t,r,c,index); 
						}
					}
				}
			}
			xt::xarray<double> x = g_params.EPS*xt::ones<double>({C,R,T});
			nb_observations = nb_observations_current;
			nb_arrivals = nb_calls_current;
			auto f_val = projected_gradient_armijo_feasible(x);
			xt::xarray<int> nb_calls_remaining = xt::zeros<int>({C,R,T});
			for(int index = 0; index < index_cross*nb_in_block; ++index){
				for(int c = 0; c < C; ++c){
					for(int r = 0; r < R; ++r){
						for(int t = 0; t < T; ++t){
							nb_calls_remaining(c,r,t) += sample(t,r,c,index);
						}
					}
				}					
			}
			for(int index = (index_cross+1)*nb_in_block; index < nb_observations_total; ++index){
				for(int c = 0; c < C; ++c){
					for(int r = 0; r < R; ++r){
						for(int t = 0; t < T; ++t){
							nb_calls_remaining(c,r,t) += sample(t,r,c,index);
						}
					}
				}						
			}
			double f = 0;
			for(int c = 0; c < C; ++c){
				for(int r = 0; r < R; ++r){
					for(int t = 0; t < T; ++t){
						double current_lambda = x(c,r,t);
						f += (nb_observations_total-nb_in_block)*current_lambda*durations[t] - 
							nb_calls_remaining(c,r,t)*log(current_lambda);
					}
				}
			}
			// fmt::print("f test set = {}\n",f);
			likelihood += f;
		}
		likelihood = likelihood / floor(1/proportion);
		// fmt::print("Likelihood current_alpha {} = {}\n", alphas[index_alpha], likelihood);
		if(likelihood < min_loss){
			min_loss = likelihood;
			best_alpha = alphas[index_alpha];
			best_weight = group_weights[index_alpha];
		}
	}
	alpha = best_alpha;
	weights = vector<double>(groups.size(), best_weight);
	xt::xarray<int> nb_observations_current = xt::zeros<int>({C,R,T});
	xt::xarray<int> nb_calls_current = xt::zeros<int>({C,R,T});
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				for(int index = 0; index < nb_observations_total; ++index){
					nb_calls_current(c,r,t) += sample(t,r,c,index);
					++nb_observations_current(c,r,t);
				}
			}
		}
	}

	nb_observations =  nb_observations_current;
	nb_arrivals = nb_calls_current;
	xt::xarray<double> x = g_params.EPS*xt::ones<double>({C,R,T});
	auto f_val = projected_gradient_armijo_feasible(x);
	auto dt = std::chrono::high_resolution_clock::now();
	cpu_time = std::chrono::duration_cast<std::chrono::seconds>(dt - t0).count();
	nb_observations = initial_nb_obs;
	nb_arrivals = initial_nb_arrivals;
	fmt::print("best_weight = {}\n", best_weight);
	return {cpu_time, best_weight, x};
}



