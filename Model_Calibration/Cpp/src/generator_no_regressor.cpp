#include "../include/generator_no_regressor.h"

using namespace std;

//Hardcoded test
GeneratorNoRegressor::GeneratorNoRegressor(){
	x_max = y_max = 10;
	n_x = n_y = 10;
	R = n_x*n_y;
	C = 1;
	int nb_groups = 4;
	T = nb_groups*7;

	which_group = vector<int>(T, 0);
	for(int t = 0; t < T; ++t){
		which_group[t] = t % nb_groups;
	}

	groups = vector<vector<int>>(nb_groups,vector<int>());
	for(int i = 0; i < nb_groups; ++i){
		groups[i].push_back(i);
	}


	for(int i = 0; i < nb_groups; ++i){
		for(int j = 0; j < 6; ++j){
			groups[i].push_back(groups[i][j]+4);
		}
		// fmt::print("group {}: {}\n",i,groups[i]);
	}
	// cin.get();

	theoretical_lambda = xt::zeros<double>({T,R,C});
	type = vector<int>(R,0);
	int count = 0;
	for(int i = 0; i < n_y/2; ++i){
		for(int j = 0; j < n_x/2; ++j){
			type[count] = 0;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.1;
			}
			++count;
		}

		for(int j = 0; j < n_x/2; ++j){
			type[count] = 1;
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
			type[count] = 1;
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
			type[count] = 0;
			for(int k = 0; k < groups[0].size(); ++k){
				theoretical_lambda(groups[0][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[1].size(); ++k){
				theoretical_lambda(groups[1][k],count,0) = 0.1;
			}
			for(int k = 0; k < groups[2].size(); ++k){
				theoretical_lambda(groups[2][k],count,0) = 0.5;
			}
			for(int k = 0; k < groups[3].size(); ++k){
				theoretical_lambda(groups[3][k],count,0) = 0.1;
			}
			++count;
		}
	}


	nb_weeks = 100000;
	nb_years = 1;
	int nb_observations_total = nb_weeks*nb_years;
	int max_obs = nb_observations_total;
	int nb_vars = 1;
	sample = xt::zeros<int>({T,R,C,static_cast<ulong>(max_obs)});
	durations = vector<double>(T,1);
	std::default_random_engine gen(600);


	// static_cast<ulong>(7*nb_weeks)
	nb_observations = xt::zeros<int>({C,R,T});
	nb_arrivals = xt::zeros<int>({C,R,T});
	// xt::random::seed(600);
	for(int t = 0; t < T; ++t){
		for(int r = 0; r < R; ++r){
			for(int c = 0; c < C; ++c){
				nb_observations(c,r,t) = nb_observations_total;
				double sum = 0;
				for(int k = 0; k < nb_observations_total; ++k){
					poisson_distribution<int> pd(theoretical_lambda(t,r,c)*durations[t]);
					sample(t,r,c,k) = pd(gen);
					// sample(t,r,c,k) = floor(theoretical_lambda(t,r,c)*durations[t]) + 1;
					sum += sample(t,r,c,k);
					// if(theoretical_lambda(t,r,c)*durations[t] != 3){
					// 	fmt::print("Sample {} {} {} {} rate = {}, sample = {}\n",t,r,c,k,theoretical_lambda(t,r,c)*durations[t],
					// 		sample(t,r,c,k));
					// 	cin.get();
					// }
					nb_arrivals(c,r,t) += sample(t,r,c,k);
				}
				// fmt::print("mean = {}\n", static_cast<double>(sum) / nb_observations_total);
				// cin.get();
				++nb_vars;
			}
		}
	}

	// for(int c = 0; c < C; ++c){
	// 	for(int t = 0; t < T; ++t){
	// 		for(int r = 0; r < R; ++r){
	// 			fmt::print("c = {}, r = {} (type {}), t = {}, lambda = {}, obs = {}, arr = {}\n",
	// 				c,r, type[r],t, theoretical_lambda(t,r,c), nb_observations(c,r,t), nb_arrivals(c,r,t));
	// 		}
	// 	}
	// }
	// cin.get();

	xt::xarray<double> est = xt::zeros<double>({2,4});
	for(int t = 0; t < groups.size(); ++t){
		int nb_type0 = 0;
		int nb_type1 = 0;
		for(int r = 0; r < R; ++r){
			if(type[r] == 0){
				est(0,t) += nb_arrivals(0,r,t);
				nb_type0 += nb_observations_total;
			}else{
				est(1,t) += nb_arrivals(0,r,t);
				nb_type1 += nb_observations_total;
			}
		}
		fmt::print("t = {}, est0 = {}, est1 = {}\n",t, est(0,t), est(1,t));
		est(0,t) =  est(0,t) / nb_type0;
		est(1,t) = 	est(1,t) / nb_type1;
		fmt::print("t = {}, est0 = {}, est1 = {}\n",t, est(0,t), est(1,t));
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

		fmt::print("r = {}, Xi = {}, Yi = {}\n",r,Xi, Yi);

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
			distance(r,s) = sqrt(pow(absc[x_r] - absc[x_s], 2) + pow(ord[y_r]-ord[y_s], 2));
			neigh.push_back(s);
		}
		fmt::print("Neighbors {}: {}\n", r, neigh);
	}
	// cin.get();


	l_bounds = xt::zeros<double>(theoretical_lambda.shape());
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
	// auto info_arq = ifstream(info_path, ios::in);
	// info_arq >> T >> G >> R >> P >> nb_land_types >> nb_holidays_years; 
	// slot_duration = 24 / T;
	// daily_obs = std::vector<int>(G, 0);
	// for(int g = 0; g < G; ++g){
	// 	info_arq >> daily_obs[g];
	// }
	// info_arq.close();

	// int max_obs = *max_element(daily_obs.begin(), daily_obs.end());

	// xt::xarray<double>::shape_type shape_d = {T,G,R,P};
	// xt::xarray<int>::shape_type shape_i = {T,G,R,P};
	// nb_observations = xt::zeros<int>(shape_i);

	// shape_i = {T,G,R,P,static_cast<ulong>(max_obs)};
	// sample = xt::zeros<int>(shape_i);

	// shape_i = {T, nb_holidays_years, R, P};
	// nb_observations_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,nb_holidays_years,G,R,P};
	// nb_arrivals_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,G,R,P};
	// nb_arrivals_no_holidays = xt::zeros<int>(shape_i);

	// auto calls_arq = ifstream(calls_path, ios::in);
	// std::string aux_str;
	// is_holidays = std::vector<std::pair<bool, int>>(max_obs, make_pair(false,-1));
	// do{
	// 	std::getline(calls_arq, aux_str);
	// 	if(aux_str == "END"){
	// 		break;
	// 	}
	// 	std::istringstream ss(aux_str);
	// 	int t, g, r, p, j, h, val;
	// 	ss >> t >> g >> r >> p >> j >> val >> h;
	// 	sample(t,g,r,p,j) = val;
	// 	if(h == -1){
	// 		nb_arrivals_no_holidays(t,g,r,p) += sample(t,g,r,p,j);
	// 	}else{
	// 		is_holidays[j] = make_pair(true, h);
	// 		nb_observations_holidays(t,h,r,p) += 1;
	// 		nb_arrivals_holidays(t,h,g,r,p) += sample(t,g,r,p,j);
	// 	}
	// 	nb_observations(t,g,r,p) += 1;
	// }while(true);

	// calls_arq.close();


	// neighbors = std::vector<vector<int>>(R, std::vector<int>());
	// xt::xarray<double>::shape_type dist_shape = {R, R};
	// regions = std::vector<Location>(R, null_location);
	// distance = xt::xarray<double>(dist_shape);
	// type = std::vector<int>(R,-1);
	// xt::xarray<double>::shape_type reg_shape = {R, nb_land_types};
	// regressors = xt::zeros<double>(reg_shape);
	// auto neighbors_arq = ifstream(neighbors_path, ios::in);
	// while(true){
	// 	int ind, terrain_type, s; 
	// 	double lat, longi, dist;
	// 	std::getline(neighbors_arq, aux_str);
	// 	if(aux_str == "END"){
	// 		break;
	// 	}
	// 	std::istringstream ss(aux_str);
	// 	ss >> ind >> lat >> longi >> terrain_type;
	// 	type[ind] = terrain_type;
	// 	regions[ind] = make_pair(lat, longi);
	// 	for(int j = 0; j < nb_land_types; ++j){
	// 		ss >> regressors(ind,j);
	// 	}
	// 	while(ss >> s >> dist){
	// 		distance(ind,s) = dist;
	// 		neighbors[ind].push_back(s);
	// 	}
	// }
	// neighbors_arq.close();
	std::cout << "Initialized\n";
}


void GeneratorNoRegressor::test(){
	double epsilon = g_params.EPS;
	// x = theoretical_lambda;
	vector<double> test_weights = {0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
		110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
		270,280,290,300,310,320,330,340,350,360,370,380,390,400};
	fmt::print("test_weights.size() = {}\n", test_weights.size());
	// vector<double> test_weights = {0.05};
	vector<double> alphas = test_weights;
	xt::xarray<double> x = epsilon*xt::ones<double>(theoretical_lambda.shape());
	ofstream err_arq(fmt::format("err_no_regressors.txt"), std::ios::out);
	for(int i = 0; i < test_weights.size(); ++i){
		x = epsilon*xt::ones<double>(theoretical_lambda.shape());
		alpha = alphas[i];
		weights = vector<double>(groups.size(), test_weights[i]);
		auto f_val = projected_gradient_armijo_feasible(x);
		fmt::print("alpha = {}, weight = {}, diff = {}\n",alpha, weights[0], average_difference(x));
		err_arq << average_difference(x) << "\n";
		// cin.get();
	}
	err_arq.close();
	fmt::print("End test_weights\n");
	cin.get();
	// ofstream x_arq(fmt::format("x_a{}_w{}.txt", alpha, weights[0]), std::ios::out);
	// vector<double> difference_l2;
	// for(int c = 0; c < C; ++c){
	// 	for(int r = 0; r < R; ++r){
	// 		for(int t = 0; t < T; ++t){
	// 			difference_l2.push_back(abs(theoretical_lambda(c,r,t)-x(c,r,t)) / theoretical_lambda(c,r,t));
	// 			err_arq << difference_l2.back() << "\n";
	// 			x_arq << c << " " << r << " " << t << " " << x(c,r,t) << "\n";
	// 		}
	// 	}
	// }
	// fmt::print("Average difference = {}\n", accumulate(difference_l2.begin(), difference_l2.end(), 0.0) / 
	// 	difference_l2.size());
	// err_arq.close();
	// x_arq.close();
	// cin.get();

	auto result = cross_validation(0.2, alphas, test_weights);
	fmt::print("Cross validation time = {}\n",result.cpu_time);
	fmt::print("Cross validation weight = {}\n",result.weight);
	x = result.lambda;
	double best_w = result.weight;
	auto f_val = projected_gradient_armijo_feasible(x);
	fmt::print("Err best weight = {}\n", average_difference(x));
	ofstream x_arq(fmt::format("x_noreg_w{}.txt",best_w), std::ios::out);
	vector<double> difference_l2;
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				difference_l2.push_back(abs(theoretical_lambda(c,r,t)-x(c,r,t)) / theoretical_lambda(c,r,t));
				err_arq << difference_l2.back() << "\n";
				x_arq << c << " " << r << " " << t << " " << x(c,r,t) << "\n";
			}
		}
	}
	x_arq.close();

	// double best_w = -1;
	// double best_a = -1;
	// double min_diff = GRB_INFINITY;
	// for(auto w: test_weights){
	// 	weights = vector<double>(groups.size(),w);
	// 	for(auto a: alphas){
	// 		alpha = a;
	// 		xt::xarray<double> x = 1*epsilon*xt::ones<double>(theoretical_lambda.shape());
	// 		auto f_val2 = projected_gradient_armijo_feasible(x);
	// 		double diff = average_difference(x);
	// 		if(diff < min_diff){
	// 			min_diff = diff;
	// 			best_w = w;
	// 			best_a = a;
	// 		}
	// 		fmt::print("weight = {}, alpha = {}, Diff = {}\n", w, alpha, min_diff);
	// 		x = 1*epsilon*xt::ones<double>(theoretical_lambda.shape());
	// 		alpha = best_a;
	// 		weights = vector<double>(groups.size(),best_w);
	// 		auto f_val = projected_gradient_armijo_feasible(x);
	// 		ofstream err_arq(fmt::format("err_a{}_w{}.txt", best_a,best_w), std::ios::out);
	// 		ofstream x_arq(fmt::format("x_a{}_w{}.txt", best_a,best_w), std::ios::out);
	// 		vector<double> difference_l2;
	// 		for(int c = 0; c < C; ++c){
	// 			for(int r = 0; r < R; ++r){
	// 				for(int t = 0; t < T; ++t){
	// 					difference_l2.push_back(abs(theoretical_lambda(c,r,t)-x(c,r,t)) / theoretical_lambda(c,r,t));
	// 					err_arq << difference_l2.back() << "\n";
	// 					x_arq << c << " " << r << " " << t << " " << x(c,r,t) << "\n";
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	//iter 100: 5.548103476539734
	// fmt::print("Min diff = {}, best w = {}, best alpha = {}\n", min_diff, best_w, best_a);
	// write_params(x);
}

std::vector<double> GeneratorNoRegressor::projected_gradient_armijo_feasible(
	xt::xarray<double>& x){
	using xt::linalg::dot;
	double eps = g_params.EPS;
	xt::xarray<double> z = xt::zeros<double>(x.shape());

	int k = 0;
	std::vector<double> f_val;
	double b_param = 2;
	double beta_k = b_param;
	f_val.reserve(max_iter);
	int j = 0;
	while(k < max_iter){
		double fold = oracle_objective_model(x);
		xt::xarray<double> gradient = oracle_gradient_model(x);
		xt::xarray<double> z = x-beta_k*gradient;
		comp_wise_max(z,eps);

		bool stop = false;
		j = 0;
		xt::xarray<double> diff_aux = x-z;
		double rhs = mat_prod(gradient, diff_aux);
		// fmt::print("fold = {} rhs = {}\n", fold, rhs);
		double f = GRB_INFINITY;
		xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
		while(!stop){
			z_aux = x + (1/pow(2,j))*(z - x);
			f = oracle_objective_model(z_aux);
			// fmt::print("\t f = {} test {}\n", f, fold -(sigma/pow(2,j))*rhs);
			if(f <= fold-(sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}
		f_val.push_back(f);
		x = z_aux;
		// fmt::print("k = {}, f = {}, j = {}, diff = {}\n", k, f, j, average_difference(x));
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
				double prev_comp = grad_component;
				double sum_neighbors = 0;
				for(int s: neighbors[r]){
					if(type[r] == type[s]){
						grad_component += 2*alpha*(x(c,r,t) - x(c,s,t)) / (distance(r,s));
						sum_neighbors  +=  2*alpha*(x(c,r,t) - x(c,s,t)) / (distance(r,s));
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
				f += nb_observations(c,r,t)*current_lambda*durations[t] - nb_arrivals(c,r,t)*log(current_lambda*durations[t]);
				for(int s: neighbors[r]){
					if(type[r] == type[s]){
						f += (0.5*alpha)*pow(x(c,r,t)- x(c,s,t), 2)/distance(r,s);
					}
				}
			}
		}
	}
	//for each c, each r, and each pair of elements inside each group
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int grindex = 0; grindex < groups.size(); ++grindex){
				auto& group = groups[grindex];
				for(int j = 0; j < group.size(); ++j){
					int t = group[j];
					for(int itp = 0; itp < group.size(); ++itp){
						int tp = group[itp];
						if(tp != t){
							f += (0.5*weights[grindex])*(pow(x(c,r,t)- x(c,r,tp), 2));
						}
					}
				}
			}
		}
	}

	return f;


	// double f = 0;
	// for(int c = 0; c < C; ++c){
	// 	for(int d = 0; d < D; ++d){
	// 		for(int t = 0; t < T; ++t){
	// 			for(int r = 0; r < R; ++r){
	// 				double rates = 0;
	// 				for(int j = 0; j < nb_regressors; ++j){
	// 					rates += x(c,d,t,j)*regressors(j,r);
	// 				}
	// 				if(rates < 0.0000001){
	// 					rates = 0.0000001;
	// 				}
	// 				f += nb_observations(c,d,t,r)*rates - 
	// 					nb_arrivals(c,d,t,r)*log(rates);
	// 			}
	// 		}
	// 	}
	// }


	// for(int g = 0; g < groups.size(); ++g){
	// 	auto& group = groups[g];
	// 	for(int c = 0; c < C; ++c){
	// 		for(auto& e1: group){
	// 			int d1 = e1.first;
	// 			int t1 = e1.second;
	// 			for(auto& e2: group){
	// 				if(e1 != e2){
	// 					int d2 = e2.first;
	// 					int t2 = e2.second;
	// 					for(int j = 0; j < nb_regressors; ++j){
	// 						f += (weights[g]/2) * pow(
	// 							(x(c,d1,t1,j)/durations[t1]) - 
	// 							(x(c,d2,t2,j)/durations[t2]),
	// 							2);
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	return f;
}


double GeneratorNoRegressor::average_difference(xt::xarray<double>& x){
	vector<double> difference_l2;
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				difference_l2.push_back(abs(theoretical_lambda(c,r,t)-x(c,r,t)) / theoretical_lambda(c,r,t));
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

void GeneratorNoRegressor::print_var(xt::xarray<double>& x){
	for(int c = 0; c < C; ++c){
		for(int r = 0; r < R; ++r){
			for(int t = 0; t < T; ++t){
				if(abs(x(c,r,t)) > g_params.EPS){
					fmt::print("nonzero: c = {}, r = {}, t = {}, val  = {}\n",c,r,t,
						x(c,r,t));
				}
			}
		}
	}
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
							nb_calls_current(c,r,t) += sample(c,r,t,index); 
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
							nb_calls_remaining += sample(c,r,t,index);
						}
					}
				}					
			}
			for(int index = (index_cross+1)*nb_in_block; index < nb_observations_total; ++index){
				for(int c = 0; c < C; ++c){
					for(int r = 0; r < R; ++r){
						for(int t = 0; t < T; ++t){
							nb_calls_remaining += sample(c,r,t,index);
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
			likelihood += f;
		}
		likelihood = likelihood / floor(1/proportion);
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
					nb_calls_current(c,r,t) += sample(c,r,t);
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


void GeneratorNoRegressor::write_params(xt::xarray<double>& x){
	// ofstream out_file(fmt::format("{}/xRegressorT{}G{}I{}P{}K{}J{}_alpha{}.txt",
	// 	g_params.generator_folder,T,G,R,P, nb_holidays_years, nb_land_types,alpha), ios::out);

	// for(int t = 0; t < T; ++t){
	// 	for(int k = 0; k < nb_holidays_years; ++k){
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				out_file << t << " " << k << " " << r << " " << p;
	// 				out_file << " " << x_delta(t,k,r,p) << "\n";
	// 			}
	// 		}
	// 	}
	// }

	// for(int t = 0; t < T; ++t){
	// 	for(int g = 0; g < G; ++g){
	// 		for(int p = 0; p < P; ++p){
	// 			for(int j = 0; j < nb_land_types; ++j){
	// 				out_file << t << " " << g << " " << p << " " << j;
	// 				out_file << " " << x(t,g,p,j) << "\n";
	// 			}
	// 		}
	// 	}
	// }

	// out_file << "END";
	// out_file.close();

	// ofstream plot(fmt::format("{}/params_plot_reg_{}.txt", g_params.generator_folder, alpha), 
	// 	ios::out);
	// for(int k = 0; k < nb_holidays_years; ++k){
	// 	plot << k << " ";
	// 	for(int t = 0; t < T; ++t){
	// 		double sum = 0;
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				sum += x_delta(t,k,r,p);
	// 			}
	// 		}
	// 		plot << sum << " ";
	// 	}
	// 	plot << "\n";
	// }
	// for(int g = 0; g < G; ++g){
	// 	plot  << g << " ";
	// 	for(int t = 0; t < T; ++t){
	// 		double sum = 0;
	// 		for(int p = 0; p < P; ++p){
	// 			for(int j = 0; j < nb_land_types; ++j){
	// 				sum += x(t,g,p,j);
	// 			}
	// 		}
	// 		plot << sum << " ";
	// 	}
	// 	plot << "\n";
	// }
	// plot.close();
}


	// for(int j = 0; j < nb_weeks*nb_years*G; ++j){
	// 	is_holidays.push_back(make_pair(false,0));
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+1) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[((year)*nb_weeks)*G + day] = make_pair(true,0);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+4) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(11 + (year)*nb_weeks)*G + day] = make_pair(true,1);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(23 + (year)*nb_weeks)*G + day] = make_pair(true,2);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(30 + (year)*nb_weeks)*G + day] = make_pair(true,3);
	// }


	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+7) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(40 + (year)*nb_weeks)*G + day] = make_pair(true,4);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+3) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(50 + (year)*nb_weeks)*G + day] = make_pair(true,5);
	// }

	// double GeneratorNoRegressor::cross_validation(xt::xarray<double>& x, 
// 	xt::xarray<double>& x_delta, double sigma, double beta_bar, 
// 		double proportion){
// 	double best_alpha = -1;
// 	// std::vector<double> alphas{0.001,0.01,0.05,0.1,0.5,1,2,5,10,50,100,1000};
// 	// // std::vector<double> alphas{0.001,0.5,1};
// 	// double max_likelihood = GRB_INFINITY;
// 	// int nb_obs = static_cast<int>(sample.shape(4));
// 	// int nb_in_block = nb_obs*proportion;
// 	// xt::xarray<double>::shape_type shape = {T,G,R,P};

// 	// double beta_tilde = 1;
// 	// double beta_hat = 2;

// 	// for(double alpha: alphas){
// 	// 	double likelihood = 0;
// 	// 	for(int ind_cross = 0; ind_cross < floor(1/proportion); ++ind_cross){
// 	// 		nb_observations = nb_in_block*xt::ones<int>(shape);
// 	// 		nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
// 	// 		nb_arrivals_holidays = xt::zeros<int>(nb_arrivals_holidays.shape());
// 	// 		nb_arrivals_no_holidays = xt::zeros<int>(nb_arrivals_no_holidays.shape());

// 	// 		for(int ind = ind_cross*nb_in_block; ind < (ind_cross+1)*nb_in_block; ++ind){
// 	// 			for(int t = 0; t < T; ++t){
// 	// 				for(int g = 0; g < G; ++g){
// 	// 					for(int r = 0; r < R; ++r){
// 	// 						for(int p = 0; p < P; ++p){
// 	// 							if(is_holidays[ind].first){
// 	// 								int k = is_holidays[ind].second;
// 	// 								nb_observations_holidays(t,k,r,p) += 1;
// 	// 								nb_arrivals_holidays(t,k,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}else{
// 	// 								nb_arrivals_no_holidays(t,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}
// 	// 						}
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}

// 	// 		auto f_val = projected_gradient_armijo_feasible(x, x_delta, alpha, 
// 	// 			sigma, beta_tilde, beta_hat);
// 	// 		xt::xarray<int> nb_observations_remaining = (nb_obs-nb_in_block)*xt::ones<int>(shape);
// 	// 		xt::xarray<int> nb_observations_holidays_remaining = xt::zeros<int>(
// 	// 			nb_observations_holidays.shape());
// 	// 		xt::xarray<int> nb_arrivals_holidays_remaining = xt::zeros<int>(nb_arrivals_holidays.shape());
// 	// 		xt::xarray<int> nb_arrivals_no_holidays_remaining = xt::zeros<int>(
// 	// 			nb_arrivals_no_holidays.shape());

// 	// 		for(int ind = 0; ind < ind_cross*nb_in_block; ++ind){
// 	// 			for(int t = 0; t < T; ++t){
// 	// 				for(int g = 0; g < G; ++g){
// 	// 					for(int r = 0; r < R; ++r){
// 	// 						for(int p = 0; p < P; ++p){
// 	// 							if(is_holidays[ind].first){
// 	// 								int k = is_holidays[ind].second;
// 	// 								nb_observations_holidays_remaining(t,k,r,p) += 1;
// 	// 								nb_arrivals_holidays_remaining(t,k,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}else{
// 	// 								nb_arrivals_no_holidays_remaining(t,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}
// 	// 						}
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}

// 	// 		for(int ind = (ind_cross+1)*nb_in_block; ind < nb_obs; ++ind){
// 	// 			for(int t = 0; t < T; ++t){
// 	// 				for(int g = 0; g < G; ++g){
// 	// 					for(int r = 0; r < R; ++r){
// 	// 						for(int p = 0; p < P; ++p){
// 	// 							if(is_holidays[ind].first){
// 	// 								int k = is_holidays[ind].second;
// 	// 								nb_observations_holidays_remaining(t,k,r,p) += 1;
// 	// 								nb_arrivals_holidays_remaining(t,k,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}else{
// 	// 								nb_arrivals_no_holidays_remaining(t,g,r,p) += sample(t,g,r,p,ind);
// 	// 							}
// 	// 						}
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}

// 	// 		double f = 0;
// 	// 		xt::xarray<double> rates = xt::zeros<double>(shape);
// 	// 		for(int t = 0; t < T; ++t){
// 	// 			for(int g = 0; g < G; ++g){
// 	// 				for(int r = 0; r < R; ++r){
// 	// 					for(int p = 0; p < P; ++p){
// 	// 						for(int j = 0; j < nb_land_types; ++j){
// 	// 							rates(t,g,r,p) += x(t,g,p,j)*regressors(r,j);
// 	// 						}
// 	// 						f += nb_observations_remaining(t,g,r,p)*rates(t,g,r,p) - 
// 	// 							nb_arrivals_no_holidays_remaining(t,g,r,p)*log(rates(t,g,r,p));
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 		for(int t = 0; t < T; ++t){
// 	// 			for(int k = 0; k < nb_holidays_years; ++k){
// 	// 				for(int r = 0; r < R; ++r){
// 	// 					for(int p = 0; p < P; ++p){
// 	// 						f += nb_observations_holidays_remaining(t,k,r,p)*x_delta(t,k,r,p);
// 	// 						for(int g = 0; g < G; ++g){
// 	// 							f -= nb_arrivals_holidays_remaining(t,k,g,r,p)*log(rates(t,g,r,p)*x_delta(t,k,r,p));
// 	// 						}
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}

// 	// 		likelihood += f;
// 	// 	}
// 	// 	fmt::print("Alpha = {}, likelihood = {}\n", alpha, likelihood);
// 	// 	if(likelihood < max_likelihood){
// 	// 		max_likelihood = likelihood;
// 	// 		best_alpha = alpha;
// 	// 	}
// 	// }
	
// 	// nb_observations = nb_obs*xt::ones<int>(shape);
// 	// nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
// 	// nb_arrivals_holidays = xt::zeros<int>(nb_arrivals_holidays.shape());
// 	// nb_arrivals_no_holidays = xt::zeros<int>(nb_arrivals_no_holidays.shape());
// 	// for(int ind = 0; ind < nb_obs; ++ind){
// 	// 	for(int t = 0; t < T; ++t){
// 	// 		for(int g = 0; g < G; ++g){
// 	// 			for(int r = 0; r < R; ++r){
// 	// 				for(int p = 0; p < P; ++p){
// 	// 					if(is_holidays[ind].first){
// 	// 						int k = is_holidays[ind].second;
// 	// 						nb_observations_holidays(t,k,r,p) += 1;
// 	// 						nb_arrivals_holidays(t,k,g,r,p) += sample(t,g,r,p,ind);
// 	// 					}else{
// 	// 						nb_arrivals_no_holidays(t,g,r,p) += sample(t,g,r,p,ind);
// 	// 					}
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 	}
// 	// }
// 	// auto f_val_best = projected_gradient_armijo_feasible(x, x_delta, best_alpha, 
// 	// 	sigma, beta_tilde, beta_hat);

// 	return best_alpha;
// }


// std::vector<double> GeneratorNoRegressor::projected_gradient_armijo_boundary(
// 	xt::xarray<double>& x, xt::xarray<double>& x_delta, double alpha, 
// 	double sigma, double beta_bar){

// 	std::vector<double> f_val;
// 	// xt::xarray<double> z_beta = xt::zeros<double>(x.shape());
// 	// xt::xarray<double> z_delta = xt::zeros<double>(x_delta.shape());
// 	// int k = 0;
// 	// double eps = 0.1;
// 	// f_val.reserve(max_iter);
// 	// while(k < max_iter){
// 	// 	double fold = oracle_objective_model2(x, x_delta, alpha);
// 	// 	xt::xarray<double> gradient_beta, gradient_delta;
// 	// 	tie(gradient_beta, gradient_delta) = oracle_gradient_model2(x, x_delta, 
// 	// 		alpha);
// 	// 	bool stop = false;
// 	// 	int j = 0;
// 	// 	double f = 0;
// 	// 	do{
// 	// 		xt::xarray<double> x_aux = x - (beta_bar/pow(2,j))*gradient_beta;
// 	// 		xt::xarray<double> x_delta_aux = x_delta - (beta_bar/pow(2,j))*gradient_delta;
// 	// 		xt::xarray<double> z_beta, z_delta;
// 	// 		tie(z_beta, z_delta) = projection_regressors(x_aux, x_delta_aux, eps);
// 	// 		f = oracle_objective_model2(z_beta, z_delta, alpha);
// 	// 		double rhs = fold - sigma*(xt::sum(gradient_beta*(x-z_beta))() + 
// 	// 			xt::sum(gradient_delta*(x_delta-z_delta))());
// 	// 		fmt::print("k = {}, j = {}, f = {}, rhs = {}\n", k, j, f, rhs);
// 	// 		if(f <= rhs){
// 	// 			stop = true;
// 	// 		}else{
// 	// 			++j;
// 	// 		}
// 	// 	}while(!stop);
// 	// 	f_val.push_back(f);
// 	// 	x = z_beta;
// 	// 	x_delta = z_delta;
// 	// 	++k;
// 	// }

// 	return f_val;
// }


// xt::xarray<double> GeneratorNoRegressor::projection_regressors(xt::xarray<double>& x){

// 	xt::xarray<GRBVar> y_beta(x.shape());

// 	GRBModel model(env);
// 	stringstream name;

// 	for(int c = 0; c < C; ++c){
// 		for(int d = 0; d < D; ++d){
// 			for(int t = 0; t < T; ++t){
// 				for(int j = 0; j < nb_regressors; ++j){
// 					name << "yb_" << c << "_" << d << "_" << t << "_" << j;
// 					double ub = pow(10,6);
// 					y_beta(c,d,t,j) = model.addVar(l_bounds(c,d,t,j), ub,
// 						0,GRB_CONTINUOUS, name.str());
// 					name.str("");
// 				}
// 			}
// 		}
// 	}

// 	GRBQuadExpr obj = 0;
// 	for(int c = 0; c < C; ++c){
// 		for(int d = 0; d < D; ++d){
// 			for(int t = 0; t < T; ++t){
// 				for(int j = 0; j < nb_regressors; ++j){
// 					obj += 0.5*y_beta(c,d,t,j)*y_beta(c,d,t,j) -
// 						x(c,d,t,j)*y_beta(c,d,t,j);
// 				}
// 			}
// 		}
// 	}
// 	try{
// 		model.setObjective(obj, GRB_MINIMIZE);
// 	}catch(GRBException& ex){
// 		cout << ex.getMessage() << "\n";
// 	}

// 	for(int c = 0; c < C; ++c){
// 		for(int d = 0; d < D; ++d){
// 			for(int t = 0; t < T; ++t){
// 				for(int r = 0; r < R; ++r){
// 					GRBLinExpr con1 = 0;
// 					for(int j = 0; j < nb_regressors; ++j){
// 						con1 += y_beta(c,d,t,j)*regressors(j,r);
// 					}
// 					name << "con1_" << c << "_" << d << "_" << t << "_" << r;
// 					model.addConstr(con1, GRB_GREATER_EQUAL, g_params.EPS, name.str());
// 					name.str("");
// 					con1 = 0;
// 				}
// 			}
// 		}
// 	}
// 	model.update();
// 	// model.write("test.lp");
// 	model.set(GRB_IntParam_OutputFlag,0);
// 	model.set(GRB_IntParam_NumericFocus, 3);
// 	model.set(GRB_IntParam_DualReductions, 0);

// 	model.optimize();

// 	auto status = model.get(GRB_IntAttr_Status);

// 	xt::xarray<double> beta_val(y_beta.shape());

// 	for(int c = 0; c < C; ++c){
// 		for(int d = 0; d < D; ++d){
// 			for(int t = 0; t < T; ++t){
// 				for(int j = 0; j < nb_regressors; ++j){
// 					beta_val(c,d,t,j) = y_beta(c,d,t,j).get(GRB_DoubleAttr_X);
// 				}
// 			}
// 		}
// 	}
// 	return beta_val;

// }