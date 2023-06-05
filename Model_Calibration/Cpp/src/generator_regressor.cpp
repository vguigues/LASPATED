#include "../include/generator_regressor.h"

using namespace std;

//Hardcoded test
GeneratorRegressor::GeneratorRegressor(GRBEnv& env): env(env){
	x_max = y_max = 10;
	n_x = n_y = 10;
	R = n_x*n_y;
	C = 1;
	T = 4;
	D = 15;

	which_group = vector<vector<int>>(D, 
		vector<int>(T, 0));

	for(int d = 0; d < 7; ++d){
		which_group[d][0] = 0;
		which_group[d][1] = 1;
		which_group[d][2] = 0;
		which_group[d][3] = 1;
	}

	for(int d = 7; d < 15; ++d){
		which_group[d][0] = 2;
		which_group[d][1] = 3;
		which_group[d][2] = 2;
		which_group[d][3] = 3;
	}

	vector<pair<int,int>> aux_group;
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,0));
	}
	groups.push_back(aux_group); //groups[0]
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,1));
	}
	groups.push_back(aux_group); //groups[1]
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,2));
	}
	groups[0].insert(groups[0].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,3));
	}
	groups[1].insert(groups[1].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();
	for(int i = 7; i < 15; ++i){
		aux_group.push_back(make_pair(i,0));
	}
	groups.push_back(aux_group); //groups[2]
	aux_group.clear();
	for(int i = 7; i < 15; ++i){
		aux_group.push_back(make_pair(i,1));
	}
	groups.push_back(aux_group); //groups[3]
	aux_group.clear();
	for(int i = 7; i < 15; ++i){
		aux_group.push_back(make_pair(i,2));
	}
	groups[2].insert(groups[2].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();
	for(int i = 7; i < 15; ++i){
		aux_group.push_back(make_pair(i,3));
	}
	groups[3].insert(groups[3].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();

	fmt::print("Groups\n");
	for(int i = 0; i < groups.size(); ++i){
		fmt::print("{}: {}\n", i, groups[i]);
	}
	fmt::print("Which Group:\n");
	for(int d = 0; d < D; ++d){
		fmt::print("{}\n", which_group[d]);
	}
	
	nb_weeks = 20;
	nb_years = floor(nb_weeks/52);
	nb_obs = nb_weeks*7;
	durations = vector<double>(T,6);
	nb_holidays_years = 8;
	is_holidays = vector<pair<bool, int>>(nb_weeks*7, make_pair(false,-1));
	vector<int> days_h;
	for(int i = 0; i < 7; ++i){
		days_h.push_back(i);
	}
	vector<int> index_h;
	for(int i = 7; i < 15; ++i){
		index_h.push_back(i);
	}

	for(int year = 0; year < nb_years; ++year){
		for(int k = 0; k < days_h.size(); ++k){
			is_holidays[year*52*7+days_h[k]] = make_pair(true,k);
		}
	}

	for(int k = 0; k < days_h.size(); ++k){
		if(nb_years*52*7 + days_h[k] <= nb_weeks*7){
			is_holidays[nb_years*52*7+days_h[k]] = make_pair(true,k);
		}
	}

	nb_land_types = 2;
	nb_regressors = 1 + nb_land_types;

	theoretical_beta = xt::zeros<double>({C,D,T,nb_regressors});
	regressors = xt::zeros<double>({nb_regressors, R});

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

	std::default_random_engine gen(600);
	std::uniform_real_distribution<double> rnd(0,1);
    int count = 0;
    type = vector<int>(R, -1);
    for(int i = 0; i < n_y / 2; ++i){
        for(int j = 0; j < n_x/2; ++j){
            type[count] = 0;
            regressors(0, count) = 50 + 50; //50 + 50*rnd(gen)
			regressors(1, count) = 0.5;
			regressors(2, count) = 0.25;
            ++count;
        } 

        for(int j = 0; j < n_x/2; ++j){
            type[count] = 1;
            regressors(0, count) = 50; // 50*rnd(gen)
			regressors(1, count) = 0.25;
			regressors(2, count) = 0.5;
            ++count;
        }
    }

    for(int i = n_y/2; i < n_y; ++i){
        for(int j = 0; j < n_x/2; ++j){
            type[count] = 1;
            regressors(0, count) = 50; // 50*rnd(gen)
			regressors(1, count) = 0.25;
			regressors(2, count) = 0.5;
            ++count;
        } 

        for(int j = 0; j < n_x/2; ++j){
            type[count] = 0;
            regressors(0, count) = 50 + 50; //50 + 50*rnd(gen)
			regressors(1, count) = 0.5;
			regressors(2, count) = 0.25;
            ++count;
        }
    }

	sample = xt::zeros<vector<int>>({C,D,T,R});
	nb_observations = xt::zeros<int>({C,D,T,R});
	nb_arrivals = xt::zeros<int>({C,D,T,R});

	int max_obs = 0;
	for(int index = 0; index < nb_weeks*7; ++index){ //each day in sample space
		int day = index % 7;

		if(is_holidays[index].first){
			day = 7 + is_holidays[index].second;
		}
		for(int c = 0; c < C; ++c){ //1
			for(int t = 0; t < T; ++t){ //4
				for(int r = 0; r < R; ++r){ // 100
					double rate = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += theoretical_beta(c,day,t,j)*regressors(j,r);
					}
					// poisson_distribution<int> pd(rate);
					int this_nb_arrival = floor(rate)+1;
					// fmt::print("Sample {} {} {} {} {}: {}\n", c,day,t,r, nb_observations(c,day,t,r),
					// 	this_nb_call);
					sample(c,day,t,r).push_back(this_nb_arrival);
                    ++nb_observations(c,day,t,r);
					nb_arrivals(c,day,t,r) += this_nb_arrival;
				}
			}
		}
	}
    
    obs_days = vector<double>(D, 0);
    obs_before = vector<double>(nb_weeks*7, 0);

    for(int index = 0; index < nb_weeks*7; ++index){
        int day = index % 7;
        if(is_holidays[index].first){
			day = 7 + is_holidays[index].second;
		}
        ++obs_days[day];
        obs_before[index] = obs_days[day];
    }

	g_params.EPS = pow(10,-6);
	sigma = 0.5;
	max_iter = 30;
    weights = vector<double>(groups.size(), 1); 
	l_bounds = xt::zeros<double>(theoretical_beta.shape());
	std::cout << "Initialized Regressors\n";
}

GeneratorRegressor::GeneratorRegressor(GRBEnv& env, std::string calls_path, 
		std::string neighbors_path, std::string info_path): env(env){
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
	// sample_calls = xt::zeros<int>(shape_i);

	// shape_i = {T, nb_holidays_years, R, P};
	// nb_observations_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,nb_holidays_years,G,R,P};
	// nb_calls_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,G,R,P};
	// nb_calls_no_holidays = xt::zeros<int>(shape_i);

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
	// 	sample_calls(t,g,r,p,j) = val;
	// 	if(h == -1){
	// 		nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,j);
	// 	}else{
	// 		is_holidays[j] = make_pair(true, h);
	// 		nb_observations_holidays(t,h,r,p) += 1;
	// 		nb_calls_holidays(t,h,g,r,p) += sample_calls(t,g,r,p,j);
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


void GeneratorRegressor::test(){
	double epsilon = g_params.EPS;
	vector<double> test_weights{0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,
		130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400};


	xt::xarray<double> x = epsilon*xt::ones<double>(theoretical_beta.shape());
	for(int i = 0; i < test_weights.size(); ++i){
		x = epsilon*xt::ones<double>(theoretical_beta.shape());
		weights = vector<double>(groups.size(), test_weights[i]);
		auto f_val = projected_gradient_armijo_feasible(x);
	}
	fmt::print("Armijo feasible alpha=1/weight=1 error: {}\n", average_difference(x));
	auto result = cross_validation(0.2,test_weights);
	fmt::print("Cross validation time = {}\n",result.cpu_time);
	fmt::print("Cross validation weight = {}\n",result.weight);
	x = result.lambda;
	double best_w = result.weight;
	ofstream err_arq(fmt::format("err_reg_a{}_w{}.txt", best_w,best_w), std::ios::out);
	ofstream x_arq(fmt::format("x_reg_a{}_w{}.txt", best_w,best_w), std::ios::out);
	vector<double> difference_l2;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rate = 0;
					double rate_est = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += theoretical_beta(c,d,t,j)*regressors(j,r);
						rate_est += x(c,d,t,j)*regressors(j,r);
					}
					difference_l2.push_back(abs(rate-rate_est) / rate);
				}
			}
		}
	}
	fmt::print("Average difference = {}\n", accumulate(difference_l2.begin(), difference_l2.end(), 0.0) / 
		difference_l2.size());


}

CrossValidationResult GeneratorRegressor::cross_validation(double proportion, 
	std::vector<double>& group_weights){
	auto t0 = std::chrono::high_resolution_clock::now();
	double cpu_time = 0;
	double best_weight = GRB_INFINITY;
	double max_likelihood = GRB_INFINITY;

	int nb_in_block = floor(nb_obs*proportion);
	xt::xarray<int> initial_nb_obs = nb_observations;
	xt::xarray<int> initial_nb_arrivals = nb_arrivals;


	for(int index_weight = 0; index_weight < group_weights.size(); ++index_weight){
		double likelihood = 0;
		weights = vector<double>(groups.size(), group_weights[index_weight]);
		for(int index_cross = 0; index_cross < floor(1/proportion); ++index_cross){
			xt::xarray<int> nb_observations_current = xt::zeros<int>({C,D,T,R});
			xt::xarray<int> nb_calls_current = xt::zeros<int>({C,D,T,R});
			for(int index = index_cross*nb_in_block; index < (index_cross+1)*nb_in_block; ++index){
				for(int c = 0; c < C; ++c){
					for(int d = 0; d < D; ++d){
						for(int t = 0; t < T; ++t){
							for(int r = 0; r < R; ++r){
								int day = index % 7;
								if(is_holidays[index].first){
									day = 7 + is_holidays[index].second;
								}
								++nb_observations_current(c,day,t,r);
								nb_calls_current(c,day,t,r) +=  sample(c,day,t,r)[obs_before[index]];
							}
						}
					}
				}
			}
			xt::xarray<double> x = g_params.EPS*xt::ones<double>({C,D,T,nb_regressors});
			nb_observations = nb_observations_current;
			nb_arrivals = nb_calls_current;
			auto f_val = projected_gradient_armijo_feasible(x);
			xt::xarray<int> nb_observations_remaining = xt::zeros<int>({C,D,T,R});
			xt::xarray<int> nb_calls_remaining = xt::zeros<int>({C,D,T,R});
			for(int index = 0; index < index_cross*nb_in_block; ++index){
				for(int c = 0; c < C; ++c){
					for(int d = 0; d < D; ++d){
						for(int t = 0; t < T; ++t){
							for(int r = 0; r < R; ++r){
								int day = index % 7;
								if(is_holidays[index].first){
									day = 7 + is_holidays[index].second;
								}
								++nb_observations_remaining(c,day,t,r);
								nb_calls_remaining(c,day,t,r) +=  sample(c,day,t,r)[obs_before[index]];
							}
						}
					}
				}
			}
			for(int index = (index_cross+1)*nb_in_block; index < nb_obs; ++index){
				for(int c = 0; c < C; ++c){
					for(int d = 0; d < D; ++d){
						for(int t = 0; t < T; ++t){
							for(int r = 0; r < R; ++r){
								int day = index % 7;
								if(is_holidays[index].first){
									day = 7 + is_holidays[index].second;
								}
								++nb_observations_remaining(c,day,t,r);
								nb_calls_remaining(c,day,t,r) +=  sample(c,day,t,r)[obs_before[index]];
							}
						}
					}
				}
			}

			double f = 0;
			for(int c = 0; c < C; ++c){
				for(int d = 0; d < D; ++d){
					for(int t = 0; t < T; ++t){
						for(int r = 0; r < R; ++r){
							double rates = 0;
							for(int j = 0; j < nb_regressors; ++j){
								rates += x(c,d,t,j)*regressors(j,r);
							}
							if(rates < 0.0000001){
								rates = 0.0000001;
							}
							f += nb_observations_remaining(c,d,t,r)*rates - 
								nb_calls_remaining(c,d,t,r)*log(rates);
						}
					}
				}
			}
			likelihood += f;
		}
		
		if(likelihood < max_likelihood){
			max_likelihood = likelihood;
			best_weight = group_weights[index_weight];
		}
	}

	weights = vector<double>(groups.size(), best_weight);
	nb_observations = xt::zeros<int>({C,D,T,R});
	nb_arrivals = xt::zeros<int>({C,D,T,R});

	for(int index = 0; index < nb_obs; ++index){
		for(int c = 0; c < C; ++c){
			for(int d = 0; d < D; ++d){
				for(int t = 0; t < T; ++t){
					for(int r = 0; r < R; ++r){
						int day = index % 7;
						if(is_holidays[index].first){
							day = 7 + is_holidays[index].second;
						}
						++nb_observations(c,day,t,r);
						nb_arrivals(c,day,t,r) +=  sample(c,day,t,r)[obs_before[index]];
					}
				}
			}
		}
	}
	xt::xarray<double> x = g_params.EPS*xt::ones<double>({C,D,T,nb_regressors});
	auto f_val = projected_gradient_armijo_feasible(x);

	auto dt = std::chrono::high_resolution_clock::now();
	cpu_time = std::chrono::duration_cast<std::chrono::seconds>(dt - t0).count();
	return {cpu_time, best_weight, x};
}



std::vector<double> GeneratorRegressor::projected_gradient_armijo_feasible(
	xt::xarray<double>& x){
	using xt::linalg::dot;
	xt::xarray<double> z = xt::zeros<double>(x.shape());
	double eps = g_params.EPS;

	int k = 0;
	std::vector<double> f_val;
	double b_param = 2;
	double beta_k = b_param;
	f_val.reserve(max_iter);
	int j = 0;
	while(k < max_iter){
		double fold = oracle_objective_model2(x);
		xt::xarray<double> gradient = oracle_gradient_model2(x);
		xt::xarray<double> x_aux = x-beta_k*gradient;
		try{
			z = projection_regressors(x_aux);
		}catch(GRBException& ex){
			fmt::print("{} {}\n",ex.getErrorCode(), ex.getMessage());
			cin.get();
		}
		bool stop = false;
		j = 0;
		xt::xarray<double> diff_aux = x-z;
		double rhs = mat_prod(gradient, diff_aux);
		fmt::print("fold = {} rhs = {}\n", fold, rhs);
		double f = GRB_INFINITY;
		xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
		while(!stop){
			z_aux = x + (1/pow(2,j))*(z-x);
			f = oracle_objective_model2(z_aux);
			fmt::print("\t f = {} test {}\n", f, fold -f -(sigma/pow(2,j))*rhs);
			if(f <= fold-(sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}
		f_val.push_back(f);
		x = z_aux;
		if(k % 1 == 0){
			fmt::print("k = {}, f = {}, j = {}\n", k, f, j);
		}
		beta_k = b_param / pow(2,j);
		++k;
		// cin.get();
	}
	fmt::print("End weight = {}: diff_l2 = {}\n", weights[0], average_difference(x));
	// cin.get();
	return f_val;
}

double GeneratorRegressor::mat_prod(xt::xarray<double>& a, xt::xarray<double>& b){
	double sum = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					sum += a(c,d,t,j)*b(c,d,t,j);
				}
			}
		}
	}
	return sum;
}



xt::xarray<double> GeneratorRegressor::oracle_gradient_model2(xt::xarray<double>& x){
	
	xt::xarray<double> gradient = xt::zeros<double>(x.shape());

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rates = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rates += x(c,d,t,j)*regressors(j,r);
					}
                    if(abs(rates) < 0.0000001){
                        rates = 0.0000001;
                    }
					for(int j = 0; j < nb_regressors; ++j){
						gradient(c,d,t,j) += nb_observations(c,d,t,r)*
							regressors(j,r)-nb_arrivals(c,d,t,r)*regressors(j,r) / 
							rates;
					}
				}
			}
		}
	}
    if(groups.size() > 0){
        for(int c = 0; c < C; ++c){
            for(int d = 0; d < D; ++d){
                for(int t = 0; t < T; ++t){
                    for(auto& e1: groups[which_group[d][t]]){
                        int d1 = e1.first;
                        int t1 = e1.second;
                        if((d1 != d) || (t1 != t)){
                            for(int j = 0; j < nb_regressors; ++j){
                                gradient(c,d,t,j) += (2*weights[which_group[d][t]] / durations[t]) *
                                    ((gradient(c,d,t,j) / 
                                    durations[t]) - (gradient(c,d1,t1,j) / 
                                    durations[t1]));
                            }
                        }
                    }
                }
            }
        }
    }

	return gradient;
}


double GeneratorRegressor::oracle_objective_model2(xt::xarray<double>& x){
	
    double f = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rates = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rates += x(c,d,t,j)*regressors(j,r);
					}
                    if(abs(rates) < 0.0000001){
                        rates = 0.0000001;
                    }
					f += nb_observations(c,d,t,r)*rates - 
						nb_arrivals(c,d,t,r)*log(rates);
				}
			}
		}
	}

	for(int m = 0; m < groups.size(); ++m){ // for m=1:length(Groups)
        auto& group = groups[m];
		for(auto& e1: group){
			int d1 = e1.first;
			int t1 = e1.second;
			for(auto& e2: group){
				if(e1 != e2){
					int d2 = e2.first;
					int t2 = e2.second;
					for(int c = 0; c < C; ++c){
						for(int j = 0; j < nb_regressors; ++j){
							f += (weights[m]/2) * pow(
								(x(c,d1,t1,j)/durations[t1]) - 
								(x(c,d2,t2,j)/durations[t2]),
								2);
						}
					}
				}
			}
		}	
	}
	return f;
}

xt::xarray<double> GeneratorRegressor::projection_regressors(xt::xarray<double>& x){

	xt::xarray<GRBVar> y(x.shape());

	GRBModel model(env);
	stringstream name;

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					name << "y_" << c << "_" << d << "_" << t << "_" << j;
					// double ub = (j == 0) ? 1 : pow(10,5);
					double ub = pow(10,5);
					y(c,d,t,j) = model.addVar(l_bounds(c,d,t,j), ub,
						0,GRB_CONTINUOUS, name.str());
					name.str("");
				}
			}
		}
	}

	GRBQuadExpr obj = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					// fmt::print("x(c{},d{},t{},j{}) = {}\n", c,d,t,j,x(c,d,t,j));
					obj += 0.5*y(c,d,t,j)*y(c,d,t,j) -
						x(c,d,t,j)*y(c,d,t,j);
				}
			}
		}
	}
	try{
		model.setObjective(obj, GRB_MINIMIZE);
	}catch(GRBException& ex){
		cout << ex.getMessage() << "\n";
	}

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					GRBLinExpr con1 = 0;
					for(int j = 0; j < nb_regressors; ++j){
						con1 += y(c,d,t,j)*regressors(j,r);
					}
					name << "con1_" << c << "_" << d << "_" << t << "_" << r;
					model.addConstr(con1, GRB_GREATER_EQUAL, g_params.EPS, name.str());
					name.str("");
					con1 = 0;
				}
			}
		}
	}
	model.update();
	// model.write("test.lp");
	model.set(GRB_IntParam_OutputFlag,0);
	model.set(GRB_IntParam_NumericFocus, 3);
	model.set(GRB_IntParam_DualReductions, 0);

	model.optimize();

	auto status = model.get(GRB_IntAttr_Status);

	xt::xarray<double> y_val = GRB_INFINITY*xt::ones<double>(y.shape());
	if(status == GRB_OPTIMAL){
		for(int c = 0; c < C; ++c){
			for(int d = 0; d < D; ++d){
				for(int t = 0; t < T; ++t){
					for(int j = 0; j < nb_regressors; ++j){
						y_val(c,d,t,j) = y(c,d,t,j).get(GRB_DoubleAttr_X);
					}
				}
			}
		}
	}else{
		fmt::print("Status = {}\n", status);
		model.write("proj_regressors.lp");
		cin.get();
	}
	return y_val;

}


double GeneratorRegressor::average_difference(xt::xarray<double>& x_beta){
	vector<double> difference_l2;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rate = 0;
					double rate_est = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += theoretical_beta(c,d,t,j)*regressors(j,r);
						rate_est += x_beta(c,d,t,j)*regressors(j,r);
					}
					difference_l2.push_back(abs(rate-rate_est) / rate);
				}
			}
		}
	}

	return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) / difference_l2.size();;
}



void GeneratorRegressor::comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps){
	// for(int t = 0; t < T; ++t){
	// 	for(int g = 0; g < G; ++g){
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				z(t,g,r,p) = max(a(t,g,r,p), eps);
	// 			}
	// 		}
	// 	}
	// }
}


bool GeneratorRegressor::is_neighbor(int r, int s){
	return r != s;
}


void GeneratorRegressor::write_params(xt::xarray<double>& x_beta){
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
	// 				out_file << " " << x_beta(t,g,p,j) << "\n";
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
	// 				sum += x_beta(t,g,p,j);
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


// std::vector<double> GeneratorRegressor::projected_gradient_armijo_boundary(
// 	xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, double alpha, 
// 	double sigma, double beta_bar){

// 	std::vector<double> f_val;
// 	// xt::xarray<double> z_beta = xt::zeros<double>(x_beta.shape());
// 	// xt::xarray<double> z_delta = xt::zeros<double>(x_delta.shape());
// 	// int k = 0;
// 	// double eps = 0.1;
// 	// f_val.reserve(max_iter);
// 	// while(k < max_iter){
// 	// 	double fold = oracle_objective_model2(x_beta, x_delta, alpha);
// 	// 	xt::xarray<double> gradient_beta, gradient_delta;
// 	// 	tie(gradient_beta, gradient_delta) = oracle_gradient_model2(x_beta, x_delta, 
// 	// 		alpha);
// 	// 	bool stop = false;
// 	// 	int j = 0;
// 	// 	double f = 0;
// 	// 	do{
// 	// 		xt::xarray<double> x_beta_aux = x_beta - (beta_bar/pow(2,j))*gradient_beta;
// 	// 		xt::xarray<double> x_delta_aux = x_delta - (beta_bar/pow(2,j))*gradient_delta;
// 	// 		xt::xarray<double> z_beta, z_delta;
// 	// 		tie(z_beta, z_delta) = projection_regressors(x_beta_aux, x_delta_aux, eps);
// 	// 		f = oracle_objective_model2(z_beta, z_delta, alpha);
// 	// 		double rhs = fold - sigma*(xt::sum(gradient_beta*(x_beta-z_beta))() + 
// 	// 			xt::sum(gradient_delta*(x_delta-z_delta))());
// 	// 		fmt::print("k = {}, j = {}, f = {}, rhs = {}\n", k, j, f, rhs);
// 	// 		if(f <= rhs){
// 	// 			stop = true;
// 	// 		}else{
// 	// 			++j;
// 	// 		}
// 	// 	}while(!stop);
// 	// 	f_val.push_back(f);
// 	// 	x_beta = z_beta;
// 	// 	x_delta = z_delta;
// 	// 	++k;
// 	// }

// 	return f_val;
// }