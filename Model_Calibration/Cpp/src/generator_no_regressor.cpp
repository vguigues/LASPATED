#include "../include/generator_no_regressor.h"

using namespace std;

//Hardcoded test
GeneratorNoRegressor::GeneratorNoRegressor(){
	x_max = y_max = 10;
	n_x = n_y = 10;
	R = n_x*n_y;
	C = 1;
	int nb_groups = 4;
	T = 4*7;

	which_group = vector<int>(T, 0);
	for(int t = 0; t < T; ++t){
		which_group[t] = t % nb_groups;
		// fmt::print("which_group {} = {}\n", t+1, which_group[t]+1);
	}
	// cin.get();
	groups = vector<vector<int>>(nb_groups,vector<int>());
	for(int i = 0; i < nb_groups; ++i){
		groups[i].push_back(i);
		// fmt::print("groups {} = ", i+1);
		// for(auto j: groups[i]){
		// 	fmt::print("{} ",j+1);
		// }
		// fmt::print("\n");
	}


	for(int i = 0; i < nb_groups; ++i){
		for(int j = 0; j < 6; ++j){
			groups[i].push_back(groups[i][j]+4);
		}
		// fmt::print("groups {} = ", i+1);
		// for(auto j: groups[i]){
		// 	fmt::print("{} ",j+1);
		// }
		// fmt::print("\n");
	}
	// cin.get();

	// for(int i = 0; i < nb_groups; ++i){
	// 	fmt::print("Groups {}: ", i+1);
	// 	for(auto j: groups[i]){
	// 		fmt::print("{} ", j+1);
	// 	}
	// 	fmt::print("\n");
	// }
	// cin.get();

	theoretical_lambda = xt::zeros<double>({T,R,C});
	type_region = vector<int>(R,-1);
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

	// for(int r = 0; r < R; ++r){
	// 	fmt::print("Region {}, type_region = {}\n", r+1, type_region[r]+1);
	// }
	// cin.get();

	nb_weeks = 10;
	nb_years = 1;
	int nb_observations_total = nb_weeks*nb_years;
	int max_obs = nb_observations_total;
	sample = xt::zeros<int>({T,R,C,static_cast<ulong>(max_obs)});
	durations = vector<double>(T,6);
	// std::default_random_engine gen(600);
	std::random_device rd;
    std::mt19937 gen(rd());


	// static_cast<ulong>(7*nb_weeks)
	nb_observations = xt::zeros<int>({C,R,T});
	nb_arrivals = xt::zeros<int>({C,R,T});
	// ifstream sample_arq("sample6.txt",std::ios::in);
	// int st, sr, sp, sk, sval;
	// while(!sample_arq.eof()){
	// 	sample_arq >> st >> sr >> sp >> sk >> sval;
	// 	sample(st-1,sr-1,sp-1,sk-1) = sval;
	// 	nb_arrivals(sp-1,sr-1,st-1) += sample(st-1,sr-1,sp-1,sk-1);
	// 	++nb_observations(sp-1,sr-1,st-1);
	// }
	// sample_arq.close();
	xt::random::seed(600);
	for(int t = 0; t < T; ++t){
		for(int r = 0; r < R; ++r){
			for(int c = 0; c < C; ++c){
				nb_observations(c,r,t) = nb_observations_total;
				double sum = 0;
				poisson_distribution<int> pd(theoretical_lambda(t,r,c)*durations[t]);
				for(int k = 0; k < nb_observations_total; ++k){
					sample(t,r,c,k) = pd(gen);
					// sample(t,r,c,k) = floor(theoretical_lambda(t,r,c)*durations[t]) + 1;
					sum += sample(t,r,c,k);
					nb_arrivals(c,r,t) += sample(t,r,c,k);
				}
				// fmt::print("mean = {}, theoretical = {}\n", static_cast<double>(sum) / nb_observations_total,
				// 	theoretical_lambda(t,r,c)*durations[t]);
				// cin.get();
			}
		}
	}

	// for(int c = 0; c < C; ++c){
	// 	for(int t = 0; t < T; ++t){
	// 		for(int r = 0; r < R; ++r){
	// 			fmt::print("c = {}, r = {} (type_region {}), t = {}, lambda = {}, obs = {}, arr = {}\n",
	// 				c+1,r+1, type_region[r]+1,t+1, theoretical_lambda(t,r,c), nb_observations(c,r,t), nb_arrivals(c,r,t));
	// 		}
	// 	}
	// }
	// cin.get();

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
		est(0,t) =  est(0,t) / nb_type0;
		est(1,t) = 	est(1,t) / nb_type1;
		// fmt::print("t = {}, est1 = {}, est2 = {}\n",t+1, est(0,t), est(1,t));
	}
	// cin.get();

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

		// fmt::print("r = {}, Xi = {}, Yi = {}\n",r,Xi, Yi);

		// int Xi = 1 + (r % n_x);
		// int Yi = r / n_y;

		// vector<pair> candidates{{Xi,Yi-1}, {Xi,Yi+1},{Xi-1,Yi},
		// 	{Xi+1,Yi}, {Xi-1,Yi+1}, {Xi-1,Yi-1}, {Xi+1,Yi+1}, {Xi+1,Yi-1}};

		// for(auto& candidate: candidates){
		// 	int x_i = candidate.first;
		// 	int y_i = candidate.second;
		// 	if(x_i >= 0 && x_i < n_x && y_i >= 0 && y_i < n_y){
		// 		neighbors[r].push_back(y_i*n_x + x_i);
		// 		distance[r][s] = 1;
		// 	}

		// }


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
		vector<int> neigh;
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
			neigh.push_back(s);
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
	fmt::print("info: {} {} {} {} {} {}\n", T,D,R,C,nb_land_types, nb_holidays_years);
	daily_obs = std::vector<int>(D, 0);
	for(int d = 0; d < D; ++d){
		info_arq >> daily_obs[d];
	}
	fmt::print("daily obs: {}\n", daily_obs);
	info_arq.close();
	nb_land_types = nb_regressors - 2;
	nb_regressors = 1 + nb_land_types;
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


	nb_regressors = nb_land_types;
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
		regressors(nb_regressors-1, ind) = pop1 + pop2;
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

	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int d = 0; d < D; ++d){
				for(int t = 0; t < T; ++t){
					int index = d*t + t;
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
	alpha = 0.1;
	durations = vector<double>(T,0.5);
	std::cout << "Initialized No Regressor Real Data\n";
}


void GeneratorNoRegressor::test(){
	double epsilon = 0.001;
	// x = theoretical_lambda;
	// vector<double> test_weights = {0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
	// 	110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
	// 	270,280,290,300,310,320,330,340,350,360,370,380,390,400};
	vector<double> test_weights = {0.1};
	vector<double> alphas = test_weights;
	xt::xarray<double> x = xt::ones<double>({C,R,T});
	// ofstream err_arq(fmt::format("err_no_reg_d{}.txt", durations[0]), std::ios::out);
	ofstream week_arq(fmt::format("week_no_reg_d{}.txt", durations[0]), std::ios::out);
	for(int i = 0; i < test_weights.size(); ++i){
		x = epsilon*xt::ones<double>({C,R,T});
		alpha = alphas[i];
		// alpha = 0;
		weights = vector<double>(groups.size(), test_weights[i]);
		// weights = vector<double>(groups.size(), 0);
		auto f_val = projected_gradient_armijo_feasible(x);
		// fmt::print("alpha = {}, weight = {}, diff = {}\n",alpha, weights[0], average_difference(x));
		vector<double> week_lambda(T, 0);
		for(int c = 0; c < C; ++c){
			for(int r = 0; r < R; ++r){
				for(int t = 0; t < T; ++t){
					week_lambda[t] += x(c,r,t);
				}
			}
		}
		for(auto lam: week_lambda){
			week_arq << fmt::format("{:.7f}\n",lam);
		}
		// err_arq << fmt::format("{:.7f}",average_difference(x)) << "\n";f
	}
	fmt::print("End test_weights\n");
	test_weights = {0,1,10,100,1000,10000};
	alphas = test_weights;
	auto result = cross_validation(0.2, alphas, test_weights);
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
				x_arq << fmt::format("{} {} {} = {}\n",c+1,r+1,t+1, result.lambda(c,r,t));
			}
		}
	}
	x_arq.close();
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
					if(type_region[r] == type_region[s]){
						grad_component += 2*alpha*(x(c,r,t) - x(c,s,t)) / (distance(r,s));
						// sum_neighbors  +=  2*alpha*(x(c,r,t) - x(c,s,t)) / (distance(r,s));
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
					if(type_region[r] == type_region[s]){
						f += (0.5*alpha)*pow(x(c,r,t)- x(c,s,t), 2)/distance(r,s);
					}
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

	for(int index_alpha = 0; index_alpha < alphas.size(); ++index_alpha){
		double likelihood = 0;
		alpha = alphas[index_alpha];
		weights = vector<double>(groups.size(), group_weights[index_alpha]);
		fmt::print("Testing weight/alpha = {} / {}\n", group_weights[index_alpha], alpha);
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
	return {cpu_time, best_weight, x};
}



