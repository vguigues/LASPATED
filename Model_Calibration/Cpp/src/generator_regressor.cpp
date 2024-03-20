#include "../include/generator_regressor.h"

using namespace std;

// Hardcoded test
GeneratorRegressor::GeneratorRegressor(GRBEnv &env) : env(env) {
  x_max = y_max = 10;
  n_x = n_y = 10;
  R = n_x * n_y;
  C = 1;
  T = 4;
  D = 7;

  nb_weeks = 100;
  nb_years = floor(nb_weeks / 52);
  nb_obs = nb_weeks * 7;
  durations = vector<double>(T, 6);
  nb_holidays_years = 8;
  is_holidays = vector<pair<bool, int>>(nb_weeks * 7, make_pair(false, -1));
  // vector<int> days_h;
  // for(int i = 0; i < 8; ++i){
  // 	days_h.push_back(i);
  // }

  // for(int year = 0; year < nb_years; ++year){
  // 	for(int k = 0; k < days_h.size(); ++k){
  // 		is_holidays[year*52*7+days_h[k]] = make_pair(true,k);
  // 	}
  // }

  // for(int k = 0; k < days_h.size(); ++k){
  // 	if(nb_years*52*7 + days_h[k] <= nb_weeks*7){
  // 		is_holidays[nb_years*52*7+days_h[k]] = make_pair(true,k);
  // 	}
  // }
  // for(size_t i = 0; i < nb_weeks*7; ++i){
  // 	auto h = is_holidays[i];
  // 	if(h.first){
  // 		fmt::print("holiday at {}, day {}\n", i, h.second);
  // 	}
  // }
  // cin.get();

  nb_land_types = 2;
  nb_regressors = 1 + nb_land_types;

  theoretical_beta = xt::zeros<double>({C, D, T, nb_regressors});
  regressors = xt::zeros<double>({nb_regressors, R});

  for (int d = 0; d < 7; ++d) {
    theoretical_beta(0, d, 1, 0) = 0.05;
    theoretical_beta(0, d, 3, 0) = 0.05;

    theoretical_beta(0, d, 0, 1) = 6;
    theoretical_beta(0, d, 1, 1) = 18;
    theoretical_beta(0, d, 2, 1) = 6;
    theoretical_beta(0, d, 3, 1) = 18;

    theoretical_beta(0, d, 0, 2) = 3;
    theoretical_beta(0, d, 1, 2) = 6;
    theoretical_beta(0, d, 2, 2) = 3;
    theoretical_beta(0, d, 3, 2) = 6;
  }

  // for(int d = 7; d < 15; ++d){
  // 	theoretical_beta(0,d,1,0) = 0.1;
  // 	theoretical_beta(0,d,3,0) = 0.1;

  // 	theoretical_beta(0,d,0,1) = 12;
  // 	theoretical_beta(0,d,1,1) = 36;
  // 	theoretical_beta(0,d,2,1) = 12;
  // 	theoretical_beta(0,d,3,1) = 36;

  // 	theoretical_beta(0,d,0,2) = 6;
  // 	theoretical_beta(0,d,1,2) = 12;
  // 	theoretical_beta(0,d,2,2) = 6;
  // 	theoretical_beta(0,d,3,2) = 12;
  // }

  // for(int c = 0; c < C; ++c){
  // 	for(int d = 0; d < D; ++d){
  // 		for(int t = 0; t < T; ++t){
  // 			for(int j = 0; j < nb_regressors; ++j){
  // 				fmt::print("{} {} {} {} = {:.2f}\n", c+1,
  // d+1,t+1, j+1, theoretical_beta(c,d,t,j));
  // 			}
  // 		}
  // 	}
  // }
  // cin.get();

  // std::default_random_engine gen;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> rnd_blue(0, 1);
  std::uniform_real_distribution<double> rnd_red(0, 1);

  // u = vector<vector<double>>(2,vector<double>(20,0));

  // for(int k = 0; k < 20; ++k){
  // 	u[0][k] = rnd_blue(gen);
  // }

  // for(int k = 0; k < 20; ++k){
  // 	u[1][k] = rnd_red(gen);
  // }

  u = {{0.2906, 0.5534, 0.9220, 0.3887, 0.3733, 0.3287, 0.0366,
        0.4025, 0.9484, 0.3023, 0.1957, 0.6760, 0.4114, 0.9227,
        0.5340, 0.4614, 0.4626, 0.2645, 0.0999, 0.4704},

       {0.3734, 0.6277, 0.2767, 0.9151, 0.0255, 0.8281, 0.2148,
        0.1626, 0.9430, 0.8117, 0.9916, 0.8755, 0.6834, 0.6933,
        0.3374, 0.7807, 0.2081, 0.0578, 0.2721, 0.3194}};

  type = vector<int>(R, -1);
  vector<bool> is_blue(100, false);
  std::uniform_real_distribution<double> rnd(0, 1);

  int count = 0;
  for (int i = 0; i < n_y / 2; ++i) {
    for (int j = 0; j < n_x / 2; ++j) {
      type[count] = 0;
      is_blue[count] = false;
      // regressors(0, count) = 50 + 50*rnd(gen);
      regressors(0, count) = get_population(count, is_blue);
      // regressors(0,count) = 1;
      regressors(1, count) = 0.5;
      regressors(2, count) = 0.25;
      ++count;
    }

    for (int j = 0; j < n_x / 2; ++j) {
      type[count] = 1;
      is_blue[count] = true;
      // regressors(0, count) = 50*rnd(gen);
      regressors(0, count) = get_population(count, is_blue);
      // regressors(0,count) = 1;
      regressors(1, count) = 0.25;
      regressors(2, count) = 0.5;
      ++count;
    }
  }

  for (int i = n_y / 2; i < n_y; ++i) {
    for (int j = 0; j < n_x / 2; ++j) {
      type[count] = 1;
      is_blue[count] = true;
      // regressors(0, count) = 50*rnd(gen);
      regressors(0, count) = get_population(count, is_blue);
      // regressors(0,count) = 1;
      regressors(1, count) = 0.25;
      regressors(2, count) = 0.5;
      ++count;
    }

    for (int j = 0; j < n_x / 2; ++j) {
      type[count] = 0;
      is_blue[count] = false;
      // regressors(0, count) = 50 + 50*rnd(gen);
      regressors(0, count) = get_population(count, is_blue);
      // regressors(0,count) = 1;
      regressors(1, count) = 0.5;
      regressors(2, count) = 0.25;
      ++count;
    }
  }

  // for(int r = 0; r < R; ++r){
  // 	fmt::print("r{}: ", r);
  // 	for(int j = 0; j < nb_regressors; ++j){
  // 		fmt::print("{} ", regressors(j,r));
  // 	}
  // 	fmt::print("\n");
  // }
  // cin.get();

  sample = xt::zeros<vector<int>>({C, D, T, R});

  nb_observations = xt::zeros<int>({C, D, T, R});
  nb_arrivals = xt::zeros<int>({C, D, T, R});
  int max_obs = 0;
  // for(int index = 0; index < nb_weeks*7; ++index){ //each day in sample space
  // 	int day = index % 7;

  // 	if(is_holidays[index].first){
  // 		day = 7 + is_holidays[index].second;
  // 	}
  // 	for(int c = 0; c < C; ++c){ //1
  // 		for(int t = 0; t < T; ++t){ //4
  // 			for(int r = 0; r < R; ++r){ // 100
  // 				double rate = 0;
  // 				for(int j = 0; j < nb_regressors; ++j){
  // 					rate +=
  // theoretical_beta(c,day,t,j)*regressors(j,r);
  // 				}
  // 				poisson_distribution<int> pd(rate);
  // 				// int this_nb_arrival = floor(rate)+1;

  // 				int this_nb_arrival = pd(gen);
  // 				// fmt::print("rate {} {} {} {} = {:.5f} |
  // thisnbcall = {}\n",
  // 				// 	c,day,t,r,rate, this_nb_arrival);
  // 				// cin.get();
  // 				// fmt::print("Sample {} {} {} {} {}: {}\n",
  // c,day,t,r, nb_observations(c,day,t,r),
  // 				// 	this_nb_call);
  // 				sample(c,day,t,r).push_back(this_nb_arrival);
  // 				// no_reg_sample((index % 7)*t, r,
  // c).push_back(this_nb_arrival);
  //                 ++nb_observations(c,day,t,r);
  // 				nb_arrivals(c,day,t,r) += this_nb_arrival;
  // 				// fmt::print("c{} d{} t{} r{}: obs = {} arr =
  // {}\n",c,day,t,r, nb_observations(c,day,t,r),
  // 				// 	nb_arrivals(c,day,t,r));
  // 				// fmt::print("{} {} {} {}: {}\n",c,day,t,r,
  // sample(c,day,t,r));
  // 				// cin.get();
  // 			}
  // 		}
  // 	}
  // }
  // fmt::print("nb_regressors = {}\n", nb_regressors);
  // cin.get();

  ifstream sample_arq("sample.txt", std::ios::in);
  int sc, sd, st, si, sk, sval;
  while (!sample_arq.eof()) {
    sample_arq >> sc >> sd >> st >> si >> sk >> sval;
    sample(sc - 1, sd - 1, st - 1, si - 1).push_back(sval);
    nb_arrivals(sc - 1, sd - 1, st - 1, si - 1) += sval;
    ++nb_observations(sc - 1, sd - 1, st - 1, si - 1);
  }
  sample_arq.close();

  // for(int c = 0; c < C; ++c){
  // 	for(int d = 0; d < D; ++d){
  // 		for(int t = 0; t < T; ++t){
  // 			for(int r = 0; r < R; ++r){
  // 				fmt::print("{} {} {} {}: ",c,d,t,r);
  // 				for(int k = 0; k < sample(c,d,t,r).size(); ++k){
  // 					fmt::print("{} ", sample(c,d,t,r)[k]);
  // 				}
  // 				fmt::print("\n");
  // 				cin.get();
  // 			}
  // 		}
  // 	}
  // }

  // Read sample
  //  ifstream arq_sample("sample.txt", std::ios::in);
  //  for(int k = 0; k < C*D*T*R*100; ++k){
  //  	int c,d,t,r, j, val;
  //  	arq_sample >> c >> d >> t >> r >> j >> val;
  //  	sample(c,d,t,r) = val;
  //  	nb_arrivals(c,d,t,r) += val;
  //  	nb_observations(c,d,t,r) += 1;
  //  }
  //  arq_sample.close();

  xt::xarray<double> no_reg_rates = xt::zeros<double>({D * T, R, C});
  std::string rates_file;
  if (D == 15) {
    rates_file = "cov_no_reg_rates_holiday.txt";
  } else if (D == 7) {
    rates_file = "cov_no_reg_rates.txt";
  }
  ofstream cov_no_reg_rates(rates_file, std::ios::out);
  for (int d = 0; d < D; ++d) {
    for (int t = 0; t < T; ++t) {
      int m_index_t = d * T + t;
      for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          cov_no_reg_rates << m_index_t << " " << r << " " << c << " " << rate
                           << "\n";
        }
      }
    }
  }
  cov_no_reg_rates.close();

  estimated = xt::zeros<double>({C, D, T, R});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          estimated(c, d, t, r) =
              nb_arrivals(c, d, t, r) / (nb_observations(c, d, t, r));
          // fmt::print("c{} d{} t{} r{}: rate = {}, estimated = {}, blue =
          // {}\n", c,d,t,r, rate, estimated(c,d,t,r), is_blue[r]);
        }
      }
    }
  }
  // cin.get();

  // for()

  g_params.EPS = pow(10, -6);
  sigma = 0.5;
  max_iter = g_params.max_iter;
  weights = vector<double>(groups.size(), 1);
  l_bounds = xt::zeros<double>(theoretical_beta.shape());
  std::cout << "Initialized Regressors\n";
}

GeneratorRegressor::GeneratorRegressor(GRBEnv &env, std::string calls_path,
                                       std::string neighbors_path,
                                       std::string info_path)
    : env(env) {
  auto info_arq = ifstream(info_path, ios::in);
  info_arq >> T >> D >> R >> C >> nb_regressors >> nb_holidays_years;
  slot_duration = 24 / T;
  fmt::print("info: {} {} {} {} {} {}\n", T, D, R, C, nb_regressors,
             nb_holidays_years);
  daily_obs = std::vector<int>(D, 0);
  for (int d = 0; d < D; ++d) {
    info_arq >> daily_obs[d];
  }
  info_arq.close();
  fmt::print("daily_obs = {}\n", daily_obs);
  nb_land_types = nb_regressors - 2;
  nb_regressors = 1 + nb_land_types;
  int max_obs = *max_element(daily_obs.begin(), daily_obs.end());
  int min_obs = *min_element(daily_obs.begin(), daily_obs.end());

  nb_observations = xt::zeros<double>({C, D, T, R});
  sample = xt::zeros<vector<int>>({C, D, T, R});
  nb_arrivals = xt::zeros<double>({C, D, T, R});

  auto calls_arq = ifstream(calls_path, ios::in);
  std::string aux_str;
  // is_holidays = std::vector<std::pair<bool, int>>(max_obs,
  // make_pair(false,-1));
  do {
    std::getline(calls_arq, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    int t, d, r, c, j, h, val;
    ss >> t >> d >> r >> c >> j >> val >> h;
    // fmt::print("calls {} {} {} {} {} {} {}\n",t,d,r,c,j,val,h);
    // cin.get();
    sample(c, d, t, r).push_back(val);
    nb_arrivals(c, d, t, r) += val;
    nb_observations(c, d, t, r) += 1;
  } while (true);
  calls_arq.close();

  durations = vector<double>(T, 0.5);
  estimated = xt::zeros<double>({C, D, T, R});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          nb_observations(c, d, t, r) = daily_obs[d];
          estimated(c, d, t, r) = static_cast<double>(nb_arrivals(c, d, t, r)) /
                                  (nb_observations(c, d, t, r) * durations[t]);
          // fmt::print("c {} d {} t{} r{}, emp = {}\n", c, d, t, r,
          //            estimated(c, d, t, r));
        }
      }
    }
  }
  ofstream est_by_t(fmt::format("est_by_t_{}.txt", R), std::ios::out);
  for (int d = 0; d < D; ++d) {
    for (int t = 0; t < T; ++t) {
      int ind_t = d * T + t;
      double sum = 0;
      for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
          sum += estimated(c, d, t, r);
        }
      }
      // fmt::print("{} ind_t = {}, agg_rate =  {}\n",
      // g_params.generator_folder,
      //            ind_t, sum);
      est_by_t << ind_t << " " << sum << "\n";
    }
  }
  est_by_t.close();
  // cin.get();
  ofstream emp_by_r(fmt::format("emp_regions_{}.txt", R), std::ios::out);

  type = std::vector<int>(R, -1);
  regressors = xt::zeros<double>({nb_regressors, R});
  neighbors = std::vector<vector<int>>(R, std::vector<int>());
  distance = xt::zeros<double>({R, R});
  regions = std::vector<Location>(R, null_location);
  auto neighbors_arq = ifstream(neighbors_path, ios::in);
  double pop1, pop2;
  while (true) {
    int ind, terrain_type, s;
    double lat, longi, dist;
    std::getline(neighbors_arq, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    ss >> ind >> lat >> longi >> terrain_type;
    type[ind] = terrain_type;
    regions[ind] = make_pair(lat, longi);
    for (int j = 0; j < nb_land_types; ++j) {
      ss >> regressors(j, ind);
    }
    ss >> pop1 >> pop2;
    regressors(nb_regressors - 1, ind) = pop1 + pop2;
    while (ss >> s >> dist) {
      distance(ind, s) = dist;
      neighbors[ind].push_back(s);
    }
  }

  // for(int r = 0; r < R; ++r){
  // 	fmt::print("r = {}, neighbors = {}, regressors = {} {} {} {}\n", r,
  // neighbors[r], 		regressors(0,r), regressors(1,r),
  // regressors(2,r), regressors(3,r));
  // }

  // cin.get();

  for (int j = 0; j < nb_regressors; ++j) {
    double sum_j = 0.0;
    for (int r = 0; r < R; ++r) {
      sum_j += regressors(j, r);
      // if (regressors(j, r) == 0.0) {
      //   regressors(j, r) = g_params.EPS;
      // }

      if (regressors(j, r) > 1000) {
        fmt::print("Big reg at {} {} = {}\n", j, r, regressors(j, r));
      }
    }
    fmt::print("j = {}, sum_j = {}\n", j, sum_j);
  }
  // cin.get();
  neighbors_arq.close();

  for (int r = 0; r < R; ++r) {
    double sum_r = 0;
    for (int c = 0; c < C; ++c) {
      double sum_c = 0;
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          int index = d * T + t;
          sum_c += estimated(c, d, t, r);
        }
      }
      sum_r += sum_c;
    }
    // fmt::print("{}, r = {}, rate = {}\n", g_params.generator_folder, r,
    // sum_r);
    emp_by_r << r << " " << regions[r].first << " " << regions[r].second << " "
             << sum_r << "\n";
  }
  emp_by_r.close();
  // cin.get();

  double epsilon = pow(10, -5);
  int t2 = regressors.shape(1);
  for (int i = 0; i < t2; ++i) {
    double sum = 0;
    for (int j = 0; j < regressors.shape(0); ++j) {
      sum += regressors(j, i);
    }
    if (sum < epsilon) {
      regressors(nb_regressors - 1, i) = epsilon;
    }
  }

  g_params.EPS = epsilon;
  sigma = 0.5;
  max_iter = g_params.max_iter;
  weights = vector<double>(groups.size(), 1);
  l_bounds = xt::zeros<double>({C, D, T, nb_regressors});
  theoretical_beta = emp_beta();
  std::cout << "Initialized Regressors Real data\n";
}

GeneratorRegressor::GeneratorRegressor(GRBEnv &env, xt::xarray<int> &N,
                                       xt::xarray<int> &M,
                                       xt::xarray<double> &reg)
    : env(env) {
  if (N.dimension() != 4) { // N should be C,D,T,R
    fmt::print(
        "Error: N has {} dimensions but should be 4. Problem was not set.\n",
        N.dimension());
    exit(1);
  }
  if (M.dimension() != 4) { // M also should be C,D,T,R
    fmt::print(
        "Error: M has {} dimensions but should be 4. Problem was not set.\n",
        M.dimension());
    exit(1);
  }

  if (reg.dimension() != 2) { // R should be nb_regressors,R
    fmt::print("Error: regressor array has {} dimensons but should be 2. "
               "Problem was not set.\n",
               reg.dimension());
    exit(1);
  }

  C = N.shape(0);
  D = N.shape(1);
  T = N.shape(2);
  R = N.shape(3);
  nb_regressors = reg.shape(0);

  nb_observations = N;
  nb_arrivals = M;
  regressors = reg;

  g_params.EPS = pow(10, -5);
  sigma = 0.5;
  max_iter = 30;
  l_bounds = xt::zeros<double>({C, D, T, nb_regressors});
}

void GeneratorRegressor::test() {
  double epsilon = g_params.EPS;
  vector<double> test_weights = g_params.weights_list;
  // vector<double> test_weights = {0,0.2,0.4, 0.6, 0.8,1.0,1.2,1.4,1.6,1.8,2};
  xt::xarray<double> x = epsilon * xt::ones<double>({C, D, T, nb_regressors});
  double min_err = 10e100;
  int min_w = -1;
  fmt::print("Running test method\n");
  // for(int i = 0; i < 1; ++i){
  // 	double w = test_weights[i];
  // 	weights = vector<double>(groups.size(), w);
  // 	x = epsilon*xt::ones<double>({C,D,T,nb_regressors});
  // 	auto f_val = projected_gradient_armijo_feasible(x);
  // 	double err = average_difference(x);
  // 	if(err < min_err){
  // 		min_err = err;
  // 		min_w = i;
  // 	}
  // }

  weights = vector<double>(groups.size(), 0);
  x = epsilon * xt::ones<double>(l_bounds.shape());
  auto f_val = projected_gradient_armijo_feasible(x);

  ofstream x_arq(fmt::format("x_reg.txt"), std::ios::out);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          double rate_est = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, r) * regressors(j, r);
          }

          x_arq << fmt::format("{} {} {} {} {}\n", c, d, t, r, rate);
        }
      }
    }
  }
  // x_arq << fmt::format("{}\n{:.7f}\n", weights[0], average_difference(x));
  x_arq.close();
  fmt::print("Wrote intensities at x_reg.txt\n");
  // fmt::print("weight {}, diff = {:.7f}\n", weights[0],
  // average_difference(x)); auto result = cross_validation(0.2,test_weights);
  // fmt::print("Cross validation time = {}\n",result.cpu_time);
  // fmt::print("Cross validation weight = {}\n",result.weight);
  // x = epsilon*xt::ones<double>(theoretical_beta.shape());
  // double best_w = result.weight;
  // f_val = projected_gradient_armijo_feasible(x);
  // fmt::print("Average difference = {}\n", average_difference(x));
}

void GeneratorRegressor::calibrate() {
  double epsilon = g_params.EPS;
  vector<double> test_weights = {0};
  xt::xarray<double> x = xt::ones<double>({C, D, T, nb_regressors});
  double min_err = 10e100;
  int min_w = -1;

  fmt::print("Running projected gradient\n");
  // ofstream err_arq(fmt::format("err_reg_obs{}_h{}.txt", nb_weeks, D == 15),
  // std::ios::out);
  for (int i = 0; i < test_weights.size(); ++i) {
    double w = test_weights[i];
    weights = vector<double>(groups.size(), w);
    x = 2 * pow(10.0, -5.0) * xt::ones<double>({C, D, T, nb_regressors});
    if (g_params.type_proj_gradient == 1) {
      fmt::print("Running projected gradient 1\n");
      auto f_val = projected_gradient_old(x);
    } else if (g_params.type_proj_gradient == 2) {
      fmt::print("Running projected gradient 2\n");
      auto f_val = projected_gradient_armijo_feasible(x);
    }
  }
  // err_arq.close();
  ofstream x_arq(fmt::format("x_reg.txt"), std::ios::out);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          x_arq << c << " " << d << " " << t << " " << j << " " << x(c, d, t, j)
                << "\n";
        }
      }
    }
  }
  x_arq.close();
  string time_rates =
      fmt::format("err_no_reg/rates_reg_real{}_g336_obs105_n1.txt", R);
  ofstream arq_rates(time_rates, std::ios::out);
  fmt::print("C\tT\tRate_Result\tRate_Emp\n");
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = (d * T) + t;
        double sum = 0;
        double sum_theo = 0;
        double sum_emp = 0;
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          double theoretical_rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
            theoretical_rate += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          sum += rate / durations[t];
          sum_theo += theoretical_rate;
          sum_emp += estimated(c, d, t, r);
        }
        arq_rates << c << " " << index_t << " " << sum << " " << sum_emp
                  << "\n";
        cout << c << " " << index_t << " " << sum << " " << sum_emp << "\n";
      }
    }
  }
  arq_rates.close();

  xt::xarray<double> beta_rect = xt::zeros<double>({C, D, T, nb_regressors});
  ifstream x_reg76("x_reg_76.txt", ios::in);
  int c, d, t, j;
  double val;
  while (true) {
    x_reg76 >> c >> d >> t >> j >> val;
    beta_rect(c, d, t, j) = val;
    // fmt::print("{} {} {} {} {}\n", c, d, t, j, val);
    if (x_reg76.eof()) {
      break;
    }
  }
  x_reg76.close();

  fmt::print("f(beta) = {}, f(x) = {}, f(rect) = {}\n",
             oracle_objective_model2(theoretical_beta),
             oracle_objective_model2(x), oracle_objective_model2(beta_rect));
  fmt::print("feasible beta? = {}\n", is_feasible(theoretical_beta));
  // fmt::print("Wrote intensities at x_reg.txt\n");
  // fmt::print("Wrote time rates at {}\n", time_rates);

  string region_rates =
      fmt::format("region_rates_real{}_g336_obs105_n1.txt", R);
  ofstream arq_regions(region_rates, std::ios::out);
  for (int r = 0; r < R; ++r) {
    vector<double> sums(C, 0.0);
    vector<double> sum_ests(C, 0.0);
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }
          sums[c] += rate;
          sum_ests[c] += estimated(c, d, t, r);
        }
      }
    }
    arq_regions << r << " ";
    for (int c = 0; c < C; ++c) {
      arq_regions << sums[c] << " ";
    }

    for (int c = 0; c < C; ++c) {
      arq_regions << sum_ests[c] << " ";
    }

    arq_regions << "\n";
  }
  arq_regions.close();
  // fmt::print("Wrote region rates at {}\n", time_rates);
}

std::vector<double>
GeneratorRegressor::projected_gradient_old(xt::xarray<double> &x) {
  int k = 0;
  std::vector<double> f_val;
  double b_param = 2.0;
  double beta_k = b_param;

  double accuracy = 0.01;
  double eps = g_params.EPS;
  double upper_lambda = 1e3;
  double upper_bound = GRB_INFINITY;
  int j = 0;

  x = projection_regressors(x, upper_lambda);
  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());
  fmt::print("MAX ITER = {}\n", max_iter);
  while (k < max_iter) {
    double fold = oracle_objective_model2(x);
    xt::xarray<double> gradient = oracle_gradient_model2(x);
    x_aux = x - beta_k * gradient;
    z = projection_regressors(x_aux, upper_lambda);
    // print_vars(z,"z");
    bool stop = false;
    j = 0;
    diff_aux = x - z;
    double rhs = mat_prod(gradient, diff_aux);
    double f = GRB_INFINITY;
    fmt::print("k = {} fold = {} rhs = {} |grad(x)| = {}, x(0,3,43,3) = {}, "
               "z(0,3,43,3) = {}\n",
               k, fold, rhs, norm_l1(gradient), x(0, 3, 43, 3), z(0, 3, 43, 3));
    while (!stop) {
      z_aux = x + (1 / pow(2, j)) * (z - x); // z_aux expected to tend to x
      f = oracle_objective_model2(z_aux);
      double test_cond = fold - (sigma / pow(2, j)) * rhs - f;
      // fmt::print("\tj = {}, f = {} test_cond = {}, 1/2j = {}, my_z = {}, my_x
      // = {}, my_zaux = {}\n",j,f, test_cond, 	(1/pow(2,j)), z(0,3,34,3),
      // x(0,3,34,3), z_aux(0,3,34,3));
      if (test_cond > -1e-3) {
        stop = true;
      } else {
        ++j;
      }
    }
    f_val.push_back(f);
    x = z_aux;
    ++k;
    beta_k = b_param / pow(2, j);
    // cin.get();
    // fold = oracle_objective_model2(x);
    // gradient = oracle_gradient_model2(x);
    // upper_bound = min(f, upper_bound);
    // double lb = get_lower_bound(x, gradient, eps, upper_lambda);
    // double lower_bound = f + lb;
    // if(upper_bound - lower_bound < accuracy){
    // 	fmt::print("k = {}, gap = {:.3f}, f = {:.3f}, lb = {:.3f}, ub =
    // {:.3f}\n", k, upper_bound - lower_bound, f, lower_bound, upper_bound);
    // 	check_solution(x,gradient, eps, upper_lambda);
    // 	break;
    // }
    // fmt::print("k = {}, gap = {:.3f}, f = {:.3f}, lb = {:.3f}, ub =
    // {:.3f}\n", k, upper_bound - lower_bound, f, lower_bound, upper_bound);
    // check_solution(x,gradient, eps, upper_lambda);
  }
  cin.get();
  return f_val;
}

std::vector<double>
GeneratorRegressor::projected_gradient_armijo_feasible(xt::xarray<double> &x) {
  int k = 0;
  std::vector<double> f_val;
  double b_param = 2.0;
  double beta_k = b_param;
  f_val.reserve(max_iter);

  double accuracy = 1;
  double upper_bound = GRB_INFINITY;

  double eps = g_params.EPS; // lower bound decision variables
  double upper_lambda = 1e3; // upper_bound
  // max_iter = 1000;

  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());

  ifstream x_reg76("x_reg_76.txt", ios::in);
  int c, d, t, j;
  double val;
  while (true) {
    x_reg76 >> c >> d >> t >> j >> val;
    x(c, d, t, j) = val;
    // fmt::print("{} {} {} {} {} {}\n", c, d, t, j, val, x(c, d, t, j));
    if (x_reg76.eof()) {
      break;
    }
  }

  x_reg76.close();

  if (!is_feasible(x)) {
    fmt::print("Betas from rectangular are not feasible\n");
    cin.get();
  }
  fmt::print("Read rectangular betas. f(rect) = {}\n",
             oracle_objective_model2(x));
  cin.get();
  if (!is_feasible(x)) {
    fmt::print("Projecting initial solution\n");
    cin.get();
    x = projection_regressors(x, upper_lambda);
  }
  // fmt::print("F(x_0) = {}, F(beta) = {}\n", oracle_objective_model2(x),
  // oracle_objective_model2(theoretical_beta)); cin.get();
  double fold = oracle_objective_model2(x);

  xt::xarray<double> gradient = oracle_gradient_model2(x);
  while (k < max_iter) {
    double pre_fold = fold;
    x_aux = x - beta_k * gradient;
    z = projection_regressors(x_aux, upper_lambda);
    // for (int c = 0; c < C; ++c) {
    //   for (int d = 0; d < D; ++d) {
    //     for (int t = 0; t < T; ++t) {
    //       for (int j = 0; j < nb_regressors; ++j) {
    //         fmt::print("z({},{},{},{}) = {}\n", c, d, t, j, z(c, d, t, j));
    //       }
    //     }
    //   }
    // }
    // cin.get();
    if (!is_feasible(z)) {
      fmt::print("Infeasible z\n");
      cin.get();
    }
    diff_aux = x - z;
    double rhs = mat_prod(gradient, diff_aux);
    fmt::print("rhs = {:.15f}\n", rhs);
    // cin.get();
    // if(rhs < -0.001){
    // 	fmt::print("NEGATIVE rhs = {}\n", rhs);
    // 	break;
    // }

    double f = oracle_objective_model2(z);
    if (rhs > 0.0 && f > fold - sigma * rhs) {
      bool stop = false;
      z_aux = xt::zeros<double>(x.shape());
      double this_pow = 1.0;
      int count = 0;
      double best_pow = -1.0;
      double best_val = f;
      while (!stop) {
        z_aux = x + (1.0 / this_pow) * (z - x);
        if (!is_feasible(z_aux)) {
          fmt::print("Infeasible z_aux\n");
          cin.get();
        }
        f = oracle_objective_model2(z_aux);
        if (f < best_val) {
          best_val = f;
          best_pow = this_pow;
        }
        xt::xarray<double> x_minus_zaux = x - z_aux;
        double test_cond = fold - (sigma / this_pow) * rhs - f;
        double sum_diff = 0;
        for (int c = 0; c < C; ++c) {
          for (int d = 0; d < D; ++d) {
            for (int t = 0; t < T; ++t) {
              for (int j = 0; j < nb_regressors; ++j) {
                double diff = abs(x(c, d, t, j) - z_aux(c, d, t, j));
                sum_diff += diff;
              }
            }
          }
        }
        fmt::print("\tj = {}, f = {}, 1/2j = {}, norm(x-z_aux) = {:.3f}, "
                   "sum_diff = {}\n",
                   count, f, (1.0 / this_pow), norm_l1(x_minus_zaux), sum_diff);
        if (f <= fold - (sigma / this_pow) * rhs) {
          stop = true;
        } else {
          this_pow *= 2;
        }
        ++count;
      }
      // vector<double> sum_js;
      // for (int j = 0; j < nb_regressors; ++j) {
      //   double sum_j = 0.0;
      //   for (int r = 0; r < R; ++r) {
      //     sum_j += regressors(j, r);
      //   }
      //   sum_js.push_back(sum_j);
      //   fmt::print("j = {}, sum = {}\n", j, sum_j);
      // }
      // double sum_diff = 0;
      // double greatest_diff = -1;
      // double f_zaux = oracle_objective_model2(z_aux);
      // fmt::print("f_zaux  = {}\n", f_zaux);
      // double sum_influence = 0.0;
      // double max_influence = -GRB_INFINITY;
      // for (int c = 0; c < C; ++c) {
      //   for (int d = 0; d < D; ++d) {
      //     for (int t = 0; t < T; ++t) {
      //       for (int j = 0; j < nb_regressors; ++j) {
      //         double diff = abs(x(c, d, t, j) - z_aux(c, d, t, j));
      //         sum_diff += diff;
      //         if (diff > greatest_diff) {
      //           greatest_diff = diff;
      //           // fmt::print("new greatest_diff at {} {} {} {}, val = {}\n",
      //           c,
      //           // d,
      //           //            t, j, x(c, d, t, j) - z_aux(c, d, t, j));
      //         }
      //         if (abs(diff) > 10.0) {
      //           fmt::print("c{} d{} t{} j{}: diff = {} | (x-z)*sum_j = {}\n",
      //           c,
      //                      d, t, j, diff, diff * sum_js[j]);
      //           cin.get();
      //         }
      //         double old_val = z_aux(c, d, t, j);
      //         z_aux(c, d, t, j) = 0.0;
      //         double curr_f = oracle_objective_model2(z_aux);
      //         double influence = abs(f_zaux - curr_f);
      //         sum_influence += influence;
      //         if (influence > max_influence) {
      //           max_influence = influence;
      //         }
      //         fmt::print("c{} d{} t{} j{}, diff = {},  influence = {}\n", c,
      //         d,
      //                    t, j, old_val, influence);
      //         z_aux(c, d, t, j) = old_val;
      //       }
      //     }
      //   }
      // }
      // fmt::print("sum_diff = {}, avg_influence = {}, max_influence = {}\n",
      //            sum_diff, sum_influence / (C * D * T * nb_regressors),
      //            max_influence);
      // cin.get();
      if (best_pow < 0.0) {
        // x = z_aux;
        f = fold;
        beta_k *= 2.0;
      } else {
        f = best_val;
        beta_k *= 2.0 / best_pow;
        x = x + (1.0 / best_pow) * (z - x);
      }
      f_val.push_back(f);
      // beta_k = beta_k*2 / this_pow;
      // beta = b_param / this_pow;
    } else {
      if (rhs > 0.0) {
        x = z;
        beta_k *= 2.0;
        // beta = b_param / this_pow;
      } else {
        f = fold;
        beta_k *= 2.0;
      }
    }

    if (!is_feasible(x)) {
      fmt::print("Infeasible x after iteration\n");
      cin.get();
    }

    fold = f;
    gradient = oracle_gradient_model2(x);

    upper_bound = min(f, upper_bound);
    double lb = get_lower_bound(x, gradient, upper_lambda);
    double lower_bound = f + lb;
    double gap = abs((upper_bound - lower_bound) / upper_bound);
    fmt::print(
        "k = {} f0 = {:.5f} rhs = {:.5f}, b_k = {}, f = {:.5f}, gap = {}\n", k,
        pre_fold, rhs, beta_k, fold, gap);
    if (gap < 0.01) {
      // fmt::print("k = {}, gap = {}, f = {}, lb = {}, ub = {} rhs = {}\n", k,
      //            gap, f, lower_bound, upper_bound, rhs);
      break;
    }
    // fmt::print("k = {}, gap = {:.3f}, f = {:.3f}, lb = {:.3f}, ub = {:.3f}
    // rhs = {}\n", k, 	gap, f, lower_bound, upper_bound, rhs);
    ++k;
    // check_solution(x,gradient, eps, upper_lambda);
    // cin.get();
  }
  cin.get();
  // vector<double> conv_aux;
  // for(int i = 0; i < 100; ++i){
  // 	double theta = i/100.0;
  // 	xt::xarray<double> diff = xt::zeros<double>(x.shape());
  // 	for(int c = 0; c < C; ++c){
  // 		for(int d = 0; d  < D; ++d){
  // 			for(int t = 0; t  < T; ++t){
  // 				for(int j = 0; j  < nb_regressors; ++j){
  // 					diff(c,d,t,j) = x(c,d,t,j) +
  // theta*(theoretical_beta(c,d,t,j)- x(c,d,t,j));
  // 					// fmt::print("beta({}, {}, {}, {}) =
  // result {} | theoretical {}\n", c,d,t,j, x(c,d,t,j),
  // 					// 	theoretical_beta(c,d,t,j));
  // 				}
  // 			}
  // 		}
  // 	}
  // 	fmt::print("theta = {}, f(diff) = {}\n", theta,
  // oracle_objective_model2(diff));
  // }
  // fmt::print("Avg_difference: {}\n", average_difference(x,
  // theoretical_beta)); cin.get();
  return f_val;
}

void GeneratorRegressor::check_solution(xt::xarray<double> &x,
                                        xt::xarray<double> &grad,
                                        double lower_lambda,
                                        double upper_lambda) {
  double avg_grad = 0.0;
  double avg_x = 0.0;
  double avg_rate = 0.0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          avg_grad += fabs(grad(c, d, t, j));
          avg_x += x(c, d, t, j);
        }
      }
    }
  }
  double avg_reg = 0.0;
  for (int j = 0; j < nb_regressors; ++j) {
    for (int r = 0; r < R; ++r) {
      avg_reg += regressors(j, r);
    }
  }

  bool feasible = true;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }
          // fmt::print("{} {} {} {}, NB_obs = {}, NB_arr = {}, Rate =
          // {:.4f}\n",c,d,t,r, 	nb_observations(c,d,t,r),
          // nb_arrivals(c,d,t,r), rate);
          avg_rate += rate;
          if (rate < g_params.EPS) {
            feasible = false;
            fmt::print("Rate unfeasible {:.4f} REGS = {:.4f} {:.4f} {:.4f} X = "
                       "{:.4f} {:.4f} {:.4f}\n",
                       rate, regressors(0, r), regressors(1, r),
                       regressors(2, r), x(c, d, t, 0), x(c, d, t, 1),
                       x(c, d, t, 2));
            goto end_rate_eval;
          }
        }
      }
    }
  }
// cin.get();
end_rate_eval:
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          if (x(c, d, t, j) < lower_lambda || x(c, d, t, j) > upper_lambda) {
            feasible = false;
            fmt::print("bound unfeasible\n");
            goto end_bound_eval;
          }
        }
      }
    }
  }
end_bound_eval:

  avg_grad /= C * D * T * nb_regressors;
  avg_x /= C * D * T * nb_regressors;
  avg_rate /= C * D * T * R;
  avg_reg /= nb_regressors * R;
  fmt::print("Avg x = {:.4f}, Avg grad = {:.4f}, Avg Rate = {:.4f}, Avg Reg = "
             "{:.4f}\n",
             avg_x, avg_grad, avg_rate, avg_reg);
}

double GeneratorRegressor::mat_prod(xt::xarray<double> &a,
                                    xt::xarray<double> &b) {
  double sum = 0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          sum += a(c, d, t, j) * b(c, d, t, j);
        }
      }
    }
  }
  return sum;
}

xt::xarray<double>
GeneratorRegressor::oracle_gradient_model2(xt::xarray<double> &x) {

  xt::xarray<double> gradient = xt::zeros<double>(x.shape());

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rates = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rates += x(c, d, t, j) * regressors(j, r);
          }
          double grad_coeff = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            grad_coeff += nb_observations(c, d, t, r) * regressors(j, r) -
                          nb_arrivals(c, d, t, r) * regressors(j, r) / rates;
            gradient(c, d, t, j) +=
                nb_observations(c, d, t, r) * regressors(j, r) -
                nb_arrivals(c, d, t, r) * regressors(j, r) / rates;
            // fmt::print("\tRegressor {} = {}\n", j, regressors(j, r));
          }
          // fmt::print("rates {} {} {} {} = {}, grad_coeff = {}, obs = {}, "
          //            "arrivals = {}\n",
          //            c, d, t, r, rates, grad_coeff, nb_observations(c, d, t,
          //            r), nb_arrivals(c, d, t, r));
        }
      }
    }
  }
  // fmt::print("C = {}, D = {}, T = {}, j = {}\n", C, D, T, nb_regressors);
  // cin.get();
  // for (int c = 0; c < C; ++c) {
  //   for (int d = 0; d < D; ++d) {
  //     for (int t = 0; t < T; ++t) {
  //       for (int j = 0; j < nb_regressors; ++j) {
  //         fmt::print("grad({}, {}, {}, {}) = {}\n", c, d, t, j,
  //                    gradient(c, d, t, j));
  //       }
  //     }
  //   }
  // }
  // cin.get();
  return gradient;
}

double GeneratorRegressor::oracle_objective_model2(xt::xarray<double> &x) {
  double f = 0;

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rates = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rates += x(c, d, t, j) * regressors(j, r);
          }
          // if(rates < g_params.EPS/10){
          // 	fmt::print("ENTERED c{} d{} t{} r{} = {} | REGS = {} {} {}\n",
          // c,d,t,r, rates, 		regressors(0,r), regressors(1,r),
          // regressors(2,r)); 	cin.get();
          //     rates = g_params.EPS/10;
          // }

          f += nb_observations(c, d, t, r) * rates -
               nb_arrivals(c, d, t, r) * log(rates);
        }
      }
    }
  }

  // for(int m = 0; m < groups.size(); ++m){ // for m=1:length(Groups)
  //     auto& group = groups[m];
  // 	for(auto& e1: group){
  // 		int d1 = e1.first;
  // 		int t1 = e1.second;
  // 		for(auto& e2: group){
  // 			if(e1 != e2){
  // 				int d2 = e2.first;
  // 				int t2 = e2.second;
  // 				for(int c = 0; c < C; ++c){
  // 					for(int j = 0; j < nb_regressors; ++j){
  // 						f += (weights[m]/2) * pow(
  // 							(x(c,d1,t1,j)/durations[t1])
  // - (x(c,d2,t2,j)/durations[t2]),
  // 2);
  // 					}
  // 				}
  // 			}
  // 		}
  // 	}
  // }
  return f;
}

xt::xarray<double>
GeneratorRegressor::projection_regressors(xt::xarray<double> &x,
                                          double upper_beta) {

  xt::xarray<GRBVar> y({C, D, T, nb_regressors});

  GRBModel model(env);
  stringstream name;

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          name << "y_" << c << "_" << d << "_" << t << "_" << j;
          double ub = (j == 0) ? 1 : upper_beta;
          // double ub = (j == 0) ? 1 : pow(10,5);
          // double ub = pow(10,5);
          y(c, d, t, j) = model.addVar(0, ub, 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }
  }

  GRBQuadExpr obj = 0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          // fmt::print("x(c{},d{},t{},j{}) = {}\n", c,d,t,j,x(c,d,t,j));
          obj += 0.5 * y(c, d, t, j) * y(c, d, t, j) -
                 x(c, d, t, j) * y(c, d, t, j);
        }
      }
    }
  }
  try {
    model.setObjective(obj, GRB_MINIMIZE);
  } catch (GRBException &ex) {
    cout << ex.getMessage() << "\n";
  }

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          GRBLinExpr con1 = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            con1 += y(c, d, t, j) * regressors(j, r);
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
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_NumericFocus, 3);
  model.set(GRB_IntParam_DualReductions, 0);
  model.optimize();

  auto status = model.get(GRB_IntAttr_Status);

  xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());
  if (status == GRB_OPTIMAL) {
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            y_val(c, d, t, j) = y(c, d, t, j).get(GRB_DoubleAttr_X);
            // fmt::print("c{} d{} t{} j{} = {:.8f}\n",c,d,t,j, y_val(c,d,t,j));
          }
        }
      }
    }
    // cin.get();
  } else {
    fmt::print("Status = {}\n", status);
    model.write("proj_regressors.lp");
    cin.get();
  }
  return y_val;
}

double GeneratorRegressor::average_difference(xt::xarray<double> &x,
                                              xt::xarray<double> &beta) {
  vector<double> difference_l2;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          double rate_est = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
            rate_est += beta(c, d, t, j) * regressors(j, r);
          }
          difference_l2.push_back(abs(rate - rate_est) / rate);
        }
      }
    }
  }

  return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) /
         difference_l2.size();
}

void GeneratorRegressor::comp_wise_max(xt::xarray<double> &z,
                                       xt::xarray<double> &a, double eps) {
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

double GeneratorRegressor::get_lower_bound(xt::xarray<double> &x,
                                           xt::xarray<double> &grad,
                                           double upper_beta) {
  GRBModel model(env);
  std::stringstream name;
  xt::xarray<GRBVar> y({C, D, T, nb_regressors});

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          name << "y_" << c << "_" << d << "_" << t << "_" << j;
          double ub = (j == 0) ? 1 : upper_beta;
          // double ub = pow(10,5);
          y(c, d, t, j) = model.addVar(0, ub, 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }
  }

  GRBLinExpr obj = 0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          // fmt::print("x(c{},d{},t{},j{}) = {:.7f}, grad = {:.3f}\n", c,d,t,j,
          // 	x(c,d,t,j),grad(c,d,t,j));
          obj += grad(c, d, t, j) * (y(c, d, t, j) - x(c, d, t, j));
        }
      }
    }
  }

  model.setObjective(obj, GRB_MINIMIZE);

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          GRBLinExpr con1 = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            con1 += y(c, d, t, j) * regressors(j, r);
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
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_NumericFocus, 3);
  model.set(GRB_IntParam_DualReductions, 0);

  model.optimize();

  auto status = model.get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
    // model.write("lower_bound_regressors.lp");
    // cin.get();
    return model.get(GRB_DoubleAttr_ObjVal);
  } else {
    fmt::print("Status = {}\n", status);
    model.write("lower_bound_regressors.lp");
    cin.get();
    return GRB_INFINITY;
  }
}

bool GeneratorRegressor::is_neighbor(int r, int s) { return r != s; }

void GeneratorRegressor::print_vars(xt::xarray<double> &x, std::string prefix) {
  fmt::print("================== Printing {} =========================\n",
             prefix);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          if (abs(x(c, d, t, j)) > g_params.EPS || t == 3) {
            fmt::print("{} {} {} {} {} = {}\n", prefix, c, d, t, j,
                       x(c, d, t, j));
          }
        }
      }
    }
  }
  fmt::print("================== End print {} =========================\n",
             prefix);
}

xt::xarray<double> laspated_reg(xt::xarray<int> &N, xt::xarray<int> &M,
                                xt::xarray<double> &reg,
                                xt::xarray<double> &x) {
  GRBEnv env;
  GeneratorRegressor gen(env, N, M, reg);
  if (x.dimension() != 4) {
    fmt::print("Error: x has {} dimensions but must be 4.\n", x.dimension());
    exit(1);
  }

  if (x.shape(0) != gen.C || x.shape(1) != gen.D || x.shape(2) != gen.T ||
      x.shape(3) != gen.nb_regressors) {
    std::string shape_x = fmt::format("({},{},{},{})", x.shape(0), x.shape(1),
                                      x.shape(2), x.shape(3));
    fmt::print("Error: x has shape {}, but expected is ({},{},{},{}).\n",
               shape_x, gen.C, gen.D, gen.T, gen.nb_regressors);
    exit(1);
  }

  auto lambda = x;
  auto f_val = gen.projected_gradient_armijo_feasible(lambda);
  return lambda;
}
double GeneratorRegressor::inner_integral(int i, int k) {
  double pi_app = 3.14159;
  double ceil_val_i = ceil((k / 5.0) * (i - 1));
  double condition_i = (5.0 / k) * ceil_val_i;
  double coeff_i = 0;

  if (condition_i <= i) {
    coeff_i =
        (5 / (pi_app * k)) * pow(-1, ceil_val_i) *
            (cos(pi_app * ceil_val_i) - cos(pi_app * (k / 5) * (i - 1))) +
        (10 / (pi_app * k)) * (floor((k / 5) * i) - ceil_val_i) -
        (5 / (pi_app * k)) * pow(-1, floor((k / 5) * (i))) *
            (cos(pi_app * (k / 5) * i) - cos(pi_app * floor((k / 5) * i)));
  } else {
    coeff_i = (5 / (pi_app * k)) * pow(-1, floor((k / 5) * (i))) *
              (cos(pi_app * (k / 5) * i) - cos(pi_app * (k / 5) * (i)));
  }
  return coeff_i;
}

double GeneratorRegressor::get_population(int r, std::vector<bool> &is_blue) {
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rnd(0, 1);
  double pi_app = 3.14159;
  // int i = (r+1) % 10;
  // int j = (r+1) / 10;
  int i = (r + 1) % 10;
  int j = 1 + floor((r + 1) / 10.0);
  if (i == 0) {
    i = 10;
    j = (r + 1) / 10;
  }
  // fmt::print("r = {}, i = {}, j = {}, blue = {}\n", r+1, i, j, is_blue[r]);
  double val = 0.0;
  for (int k = 1; k <= 10; ++k) {
    double ceil_val_i = ceil((k / 5.0) * (i - 1)); // ceil_val_i == int3
    double floor_val_i = floor((k / 5.0) * i);     // floor_val_i == int4
    double condition_i = (5.0 / k) * ceil_val_i;
    double coeff_i = 0;
    if (condition_i <= i) {
      double pw1 = pow(-1.0, ceil_val_i);
      coeff_i =
          (5.0 / (pi_app * k)) * pw1 *
              (cos(pi_app * ceil_val_i) - cos(pi_app * (k / 5.0) * (i - 1))) +
          (10.0 / (pi_app * k)) * (floor_val_i - ceil_val_i) -
          (5.0 / (pi_app * k)) * pow(-1.0, floor_val_i) *
              (cos(pi_app * k * (i / 5.0)) - cos(pi_app * floor_val_i));
    } else {
      coeff_i =
          (5.0 / (pi_app * k)) * pow(-1.0, floor_val_i) *
          (cos(pi_app * k * ((i - 1) / 5.0)) - cos(pi_app * k * (i / 5.0)));
    }

    double ceil_val_j = ceil((k / 5.0) * (j - 1));
    double floor_val_j = floor((k / 5.0) * j);
    double condition_j = (5.0 / k) * ceil_val_j;
    double coeff_j = 0;
    if (condition_j <= j) {
      double pw1 = pow(-1.0, ceil_val_j);
      coeff_j =
          (5.0 / (pi_app * k)) * pw1 *
              (cos(pi_app * ceil_val_j) - cos(pi_app * (k / 5.0) * (j - 1))) +
          (10.0 / (pi_app * k)) * (floor_val_j - ceil_val_j) -
          (5.0 / (pi_app * k)) * pow(-1.0, floor_val_j) *
              (cos(pi_app * (k / 5.0) * j) - cos(pi_app * floor_val_j));
    } else {
      coeff_j =
          (5.0 / (pi_app * k)) * pow(-1.0, floor_val_j) *
          (cos(pi_app * (k / 5.0) * (j - 1)) - cos(pi_app * (k / 5.0) * (j)));
    }
    if (is_blue[r]) {
      val += (1.0 / pow(2.0, k)) *
             (u[0][(2 * k) - 2] * coeff_i + u[0][2 * k - 1] * coeff_j);
    } else {
      val += (1.0 / pow(2.0, k)) * ((u[1][(2 * k) - 2] + 1) * coeff_i +
                                    (u[1][2 * k - 1] + 1) * coeff_j);
    }
    // fmt::print("\tint1 = {}, int2 = {}\n", coeff_i, coeff_j);
  }
  // cin.get();
  return val;
}

xt::xarray<double> GeneratorRegressor::emp_beta() {
  GRBModel model(env);

  xt::xarray<double> beta = xt::zeros<double>({C, D, T, nb_regressors});
  xt::xarray<GRBVar> x = xt::xarray<GRBVar>({C, D, T, nb_regressors});

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          x(c, d, t, j) =
              model.addVar(0.0, 10e5, 0, GRB_CONTINUOUS,
                           fmt::format("beta_{}_{}_{}_{}", c, d, t, j));
        }
      }
    }
  }

  GRBQuadExpr obj = 0;

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          GRBLinExpr rate = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }

          obj +=
              (rate - estimated(c, d, t, r)) * (rate - estimated(c, d, t, r));
        }
      }
    }
  }
  model.setObjective(obj, GRB_MINIMIZE);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          GRBLinExpr rate = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }

          model.addConstr(rate, GRB_GREATER_EQUAL, g_params.EPS,
                          fmt::format("rate_{}_{}_{}_{}", c, d, t, r));
        }
      }
    }
  }
  model.set(GRB_IntParam_OutputFlag, 0);
  model.optimize();

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          beta(c, d, t, j) = x(c, d, t, j).get(GRB_DoubleAttr_X);
        }
      }
    }
  }

  return beta;
}

double GeneratorRegressor::norm_l1(xt::xarray<double> &x) {
  double val = 0.0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          val += abs(x(c, d, t, j));
        }
      }
    }
  }

  return val;
}

bool GeneratorRegressor::is_feasible(xt::xarray<double> &x) {
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }
          // if (c == 0 && d == 0 && t == 0 && r == 62) {
          //   fmt::print("Rate {} {} {} {} = {}\n", c, d, t, r, rate);
          //   cin.get();
          // }
          if (rate < g_params.EPS) {
            fmt::print("Rate {} {} {} {} = {} is unfeasible\n", c, d, t, r,
                       rate);
            return false;
          }
        }
      }
    }
  }

  return true;
}

void GeneratorRegressor::write_params(xt::xarray<double> &x_beta) {
  // ofstream
  // out_file(fmt::format("{}/xRegressorT{}G{}I{}P{}K{}J{}_alpha{}.txt",
  // 	g_params.generator_folder,T,G,R,P, nb_holidays_years,
  // nb_land_types,alpha), ios::out);

  // for(int t = 0; t < T; ++t){
  // 	for(int k = 0; k < nb_holidays_years; ++k){
  // 		for(int r = 0; r < R; ++r){
  // 			for(int p = 0; p < P; ++p){
  // 				out_file << t << " " << k << " " << r << " " <<
  // p; 				out_file << " " << x_delta(t,k,r,p) <<
  // "\n";
  // 			}
  // 		}
  // 	}
  // }

  // for(int t = 0; t < T; ++t){
  // 	for(int g = 0; g < G; ++g){
  // 		for(int p = 0; p < P; ++p){
  // 			for(int j = 0; j < nb_land_types; ++j){
  // 				out_file << t << " " << g << " " << p << " " <<
  // j; 				out_file << " " << x_beta(t,g,p,j) <<
  // "\n";
  // 			}
  // 		}
  // 	}
  // }

  // out_file << "END";
  // out_file.close();

  // ofstream plot(fmt::format("{}/params_plot_reg_{}.txt",
  // g_params.generator_folder, alpha), 	ios::out); for(int k = 0; k <
  // nb_holidays_years; ++k){ 	plot << k << " "; 	for(int t = 0; t < T;
  // ++t){ 		double sum = 0; 		for(int r = 0; r < R;
  // ++r){ 			for(int p = 0; p < P;
  // ++p){ 				sum += x_delta(t,k,r,p);
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

// 	which_group = vector<vector<int>>(D,
// 	vector<int>(T, 0));

// for(int d = 0; d < 7; ++d){
// 	which_group[d][0] = 0;
// 	which_group[d][1] = 1;
// 	which_group[d][2] = 0;
// 	which_group[d][3] = 1;
// }

// // for(int d = 7; d < 15; ++d){
// // 	which_group[d][0] = 2;
// // 	which_group[d][1] = 3;
// // 	which_group[d][2] = 2;
// // 	which_group[d][3] = 3;
// // }

// vector<pair<int,int>> aux_group;
// for(int i = 0; i < 7; ++i){
// 	aux_group.push_back(make_pair(i,0));
// }
// groups.push_back(aux_group); //groups[0]
// aux_group.clear();
// for(int i = 0; i < 7; ++i){
// 	aux_group.push_back(make_pair(i,1));
// }
// groups.push_back(aux_group); //groups[1]
// aux_group.clear();
// for(int i = 0; i < 7; ++i){
// 	aux_group.push_back(make_pair(i,2));
// }
// groups[0].insert(groups[0].end(), aux_group.begin(),
// 	aux_group.end());
// aux_group.clear();
// for(int i = 0; i < 7; ++i){
// 	aux_group.push_back(make_pair(i,3));
// }
// groups[1].insert(groups[1].end(), aux_group.begin(),
// 	aux_group.end());
// aux_group.clear();

// for(int i = 7; i < 15; ++i){
// 	aux_group.push_back(make_pair(i,0));
// }
// groups.push_back(aux_group); //groups[2]
// aux_group.clear();
// for(int i = 7; i < 15; ++i){
// 	aux_group.push_back(make_pair(i,1));
// }
// groups.push_back(aux_group); //groups[3]
// aux_group.clear();
// for(int i = 7; i < 15; ++i){
// 	aux_group.push_back(make_pair(i,2));
// }
// groups[2].insert(groups[2].end(), aux_group.begin(),
// 	aux_group.end());
// aux_group.clear();
// for(int i = 7; i < 15; ++i){
// 	aux_group.push_back(make_pair(i,3));
// }
// groups[3].insert(groups[3].end(), aux_group.begin(),
// 	aux_group.end());
// aux_group.clear();

// fmt::print("Groups\n");
// for(int i = 0; i < groups.size(); ++i){
// 	fmt::print("{}: ", i+1);
// 	for(auto elem: groups[i]){
// 		fmt::print("({}, {}) ", elem.first+1, elem.second+1);
// 	}
// 	fmt::print("\n");
// }
// fmt::print("Which Group:\n");
// for(int d = 0; d < D; ++d){
// 	fmt::print("{}: ",d+1);
// 	for(int t = 0; t < which_group[d].size(); ++t){
// 		fmt::print("{} ", which_group[d][t]+1);
// 	}
// 	fmt::print("\n");
// }
// cin.get();

// int count = 0;
// for(int i = 0; i < n_y / 2; ++i){
//     for(int j = 0; j < n_x/2; ++j){
//         type[count] = 0;
// 		is_blue[count] = false;
//         // regressors(0, count) = 50 + 50*rnd(gen);
// 		regressors(0,count) = get_population(count, is_blue);
// 		regressors(1, count) = 0.5;
// 		regressors(2, count) = 0.25;
//         ++count;
//     }

//     for(int j = 0; j < n_x/2; ++j){
//         type[count] = 1;
// 		is_blue[count] = true;
//         // regressors(0, count) = 50*rnd(gen);
// 		regressors(0,count) = get_population(count, is_blue);
// 		regressors(1, count) = 0.25;
// 		regressors(2, count) = 0.5;
//         ++count;
//     }
// }

// for(int i = n_y/2; i < n_y; ++i){
//     for(int j = 0; j < n_x/2; ++j){
//         type[count] = 1;
// 		is_blue[count] = true;
//         // regressors(0, count) = 50*rnd(gen);
// 		regressors(0,count) = get_population(count, is_blue);
// 		regressors(1, count) = 0.25;
// 		regressors(2, count) = 0.5;
//         ++count;
//     }

//     for(int j = 0; j < n_x/2; ++j){
//         type[count] = 0;
// 		is_blue[count] = false;
//         // regressors(0, count) = 50 + 50*rnd(gen);
// 		regressors(0,count) = get_population(count, is_blue);
// 		regressors(1, count) = 0.5;
// 		regressors(2, count) = 0.25;
//         ++count;
//     }
// }