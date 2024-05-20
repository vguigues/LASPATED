#include "../include/generator_no_regressor.h"

using namespace std;

// Hardcoded test
GeneratorNoRegressor::GeneratorNoRegressor() {
  x_max = y_max = 10;
  n_x = n_y = 10;
  R = n_x * n_y;
  C = 1;
  T = 4 * 7;
  // T = 4 * 15;
  int nb_groups = 2;
  nb_weeks = 1;
  full_neighbors = false;  // if true, allows diagonal neighbors
  neighbor_factor = 1;     // If true, apply neighborhood regularization
  read_covariates = false; // If true, read lambdas from covariates file
  use_simulation = false;  // If true, lambdas are generated via simulation
  constant_lambdas = true; // If true, lambda(s,t) is constant

  which_group = vector<int>(T, -1);
  groups = vector<vector<int>>(nb_groups, vector<int>());
  for (int t = 0; t < T; ++t) {
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

  type_region = vector<int>(R, -1);

  nb_years = 1;
  nb_observations_total = nb_weeks * nb_years;
  int max_obs = nb_observations_total;
  durations = vector<double>(T, 1);
  // std::default_random_engine gen(600);
  std::random_device rd;
  std::mt19937 gen(rd());

  is_red = vector<bool>(100, false);

  for (int r = 0; r < R; ++r) {
    int x = r % 10;
    int y = r / 10;
    double cx = x + 0.5;
    double cy = y + 0.5;
    int mid_x = n_x / 2;
    int mid_y = n_y / 2;

    if ((x < mid_x && y < mid_y) || (x >= mid_x && y >= mid_y)) {
      is_red[r] = true;
      type_region[r] = 0;
    } else if ((x >= mid_x && y < mid_y) || (x < mid_x && y >= mid_y)) {
      is_red[r] = false;
      type_region[r] = 1;
    } else {
      fmt::print("ERROR: Impossible region {}, ({},{})\n", r, x, y);
      exit(1);
    }
    // fmt::print("Region {}: ({}, {}) blue = {}, type = {}\n", r+1, cx, cy,
    // !is_red[r], type_region[r]); if(is_red[r] && type_region[r] != 0){
    // 	fmt::print("Region {} is red but type == 1\n", r);
    // 	cin.get();
    // }
  }
  // cin.get();

  sample = xt::zeros<int>({T, R, C, static_cast<ulong>(max_obs)});
  nb_arrivals = xt::zeros<int>({C, R, T});
  uniform_real_distribution<double> rand_coord(0, 1);
  nb_observations = nb_observations_total * xt::ones<int>({C, R, T});
  theoretical_lambda = xt::zeros<double>({T, R, C});
  if (read_covariates) {
    std::string rates_file;
    if (T == 4 * 15) {
      rates_file = "cov_no_reg_rates_holiday.txt";
    } else {
      rates_file = "cov_no_reg_rates.txt";
    }
    ifstream arq_cov(rates_file, ios::in);
    xt::xarray<double> lambda_cov = xt::zeros<double>({T, R, C});
    for (int k = 0; k < T * R * C; ++k) {
      int t, r, c;
      double val;
      arq_cov >> t >> r >> c >> val;
      // fmt::print("lam_trc {} {} {} {}\n",t,r,c,val);
      lambda_cov(t, r, c) = val;
      theoretical_lambda(t, r, c) = lambda_cov(t, r, c);
    }

    arq_cov.close();
    for (int t = 0; t < T; ++t) {
      for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
          for (int k = 0; k < nb_observations_total; ++k) {
            poisson_distribution<int> pd_cov(lambda_cov(t, r, c));
            int this_nb_arrival = pd_cov(gen);
            sample(t, r, c, k) = this_nb_arrival;
            nb_arrivals(c, r, t) += this_nb_arrival;
          }
        }
      }
    }
  } else {
    for (int t = 0; t < T; ++t) {
      for (int r = 0; r < R; ++r) {
        int x = r % 10;
        int y = r / 10;
        double cx = x + 0.5;
        double cy = y + 0.5;
        for (int c = 0; c < C; ++c) {
          if ((!is_red[r] && which_group[t] % 2 == 0) ||
              (is_red[r] && which_group[t] % 2 != 0)) {
            theoretical_lambda(t, r, c) = (constant_lambdas) ? 0.1 : (cx + cy);
          } else if ((!is_red[r] && which_group[t] % 2 != 0) ||
                     (is_red[r] && which_group[t] % 2 == 0)) {
            theoretical_lambda(t, r, c) =
                (constant_lambdas) ? 0.5 : 5 * (cx + cy);
          } else {
            fmt::print("ERROR: Unpredictable lambda case at t = {} and r = {} "
                       "({},{})\n",
                       t, r, cx, cy);
            exit(1);
          }
          for (int n = 0; n < nb_observations_total; ++n) {
            poisson_distribution<int> pd(theoretical_lambda(t, r, c) *
                                         durations[t]);
            int this_nb_arrival = pd(gen);
            sample(t, r, c, n) = this_nb_arrival;
            nb_arrivals(c, r, t) += this_nb_arrival;
          }
        }
      }
    }
  }

  estimated = xt::zeros<double>({C, R, T});
  vector<double> aux;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      int x = r % 10;
      int y = r / 10;
      double cx = x + 0.5;
      double cy = y + 0.5;
      for (int t = 0; t < T; ++t) {
        estimated(c, r, t) =
            nb_arrivals(c, r, t) /
            static_cast<double>(nb_observations(c, r, t) * durations[t]);
        aux.push_back(abs(estimated(c, r, t) - theoretical_lambda(t, r, c)));
        // fmt::print("R = {} ({},{}), blue = {}, t = {}, est = {}, teorico =
        // {}\n",r + 1, cx,cy, !is_red[r], t, estimated(c,r,t),
        // 	theoretical_lambda(t,r,c));
      }
      // cin.get();
    }
  }
  double diff_est_theoretical =
      accumulate(aux.begin(), aux.end(), 0.0) / aux.size();
  aux.clear();

  // fmt::print("Avg diff est and theoretical {}\n", diff_est_theoretical);
  // cin.get();

  string region_est_file;
  if (read_covariates) {
    region_est_file = "est_by_region_cov.txt";
  } else {
    region_est_file = "est_by_region.txt";
  }
  ofstream avg_by_region(region_est_file, ios::out);
  for (int r = 0; r < R; ++r) {
    double sum = 0.0;
    int count = 0;
    for (int c = 0; c < C; ++c) {
      for (int t = 0; t < T; ++t) {
        sum += estimated(c, r, t);
        ++count;
      }
    }
    avg_by_region << sum / count << "\n";
  }
  avg_by_region.close();

  xt::xarray<double> est = xt::zeros<double>({2, 4});
  for (int t = 0; t < 4; ++t) {
    int nb_type0 = 0;
    int nb_type1 = 0;
    for (int r = 0; r < R; ++r) {
      if (type_region[r] == 0) {
        est(0, t) += nb_arrivals(0, r, t);
        nb_type0 += nb_observations_total;
      } else {
        est(1, t) += nb_arrivals(0, r, t);
        nb_type1 += nb_observations_total;
      }
    }
    // fmt::print("t = {}, est1 = {}, est2 = {}\n",t+1, est(0,t), est(1,t));
    est(0, t) = est(0, t) / (nb_type0 * durations[t]);
    est(1, t) = est(1, t) / (nb_type1 * durations[t]);
    // fmt::print("t = {}, est1 = {}, est2 = {}\n",t+1, est(0,t), est(1,t));
  }

  neighbors = vector<vector<int>>(R, vector<int>());
  distance = GRB_INFINITY * xt::ones<double>({R, R});
  for (int r = 0; r < R; ++r) {
    int Xi = (r + 1) % n_x;
    int Yi = -1;
    if (Xi == 0) {
      Yi = (r + 1) / n_x;
    } else {
      Yi = 1 + (r + 1 - Xi) / n_x;
    }
    if (Xi == 0) {
      Xi = n_x;
    }

    if (Yi - 1 >= 1) {
      neighbors[r].push_back((Yi - 2) * n_x + Xi - 1);
    }

    if (Yi + 1 <= n_y) {
      neighbors[r].push_back(Yi * n_x + Xi - 1);
    }

    if (Xi - 1 >= 1) {
      neighbors[r].push_back((Yi - 1) * n_x + Xi - 1 - 1);
    }

    if (Xi + 1 <= n_x) {
      neighbors[r].push_back((Yi - 1) * n_x + Xi + 1 - 1);
    }

    if (full_neighbors) {
      if ((Yi + 1 <= n_y) && (Xi - 1 >= 1)) {
        neighbors[r].push_back(Yi * n_x + Xi - 1 - 1);
      }

      if ((Yi - 1 >= 1) && (Xi - 1) >= 1) {
        neighbors[r].push_back((Yi - 2) * n_x + Xi - 1 - 1);
      }

      if ((Yi + 1 <= n_y) && (Xi + 1) <= n_x) {
        neighbors[r].push_back(Yi * n_x + Xi + 1 - 1);
      }

      if ((Yi - 1 >= 1) && (Xi + 1) <= n_x) {
        neighbors[r].push_back((Yi - 2) * n_x + Xi + 1 - 1);
      }
    }
  }

  vector<double> absc(n_x, 0);
  for (int i = 0; i < n_x; ++i) {
    absc[i] = x_max / (2 * n_x) + (x_max / n_x) * i;
  }
  vector<double> ord(n_x, 0);
  for (int j = 0; j < n_y; ++j) {
    ord[j] = y_max / (2 * n_y) + (y_max / n_y) * j;
  }

  for (int r = 0; r < R; ++r) {
    int x_r = r % n_x;
    int y_r = r / n_y;

    for (int s = 0; s < R; ++s) {
      int x_s = s % n_x;
      int y_s = s / n_y;

      // distance(r,s) = sqrt(pow(x_r - x_s, 2) + pow(y_r - y_s, 2));
      // if(x_r == x_s || y_r == y_s){
      // 	distance(r,s) = 1;
      // }else{
      // 	distance(r,s) = pow(2,0.5);
      // }
      // distance(r,s) = sqrt(pow(absc[x_r] - absc[x_s], 2) +
      // pow(ord[y_r]-ord[y_s], 2));
      distance(r, s) = 1;
    }
    // fmt::print("Neighbors {} (type_region = {}): ", r+1, type_region[r]+1);
    // for(auto s: neighbors[r]){
    // 	fmt::print("{} ", s+1);
    // }
    // fmt::print("\n");
  }
  // cin.get();

  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        cout << c << " " << r << " " << t
             << ": obs = " << nb_observations(c, r, t)
             << ", arr = " << nb_arrivals(c, r, t)
             << ", est = " << estimated(c, r, t) << "\n";
      }
    }
  }
  cin.get();

  l_bounds = xt::zeros<double>({C, R, T});
  g_params.EPS = pow(10, -3);
  alpha = 1;
  weight = 1;
  weights = vector<double>(nb_groups, 1);
  sigma = 0.5;
  beta_bar = 1;
  max_iter = 30;
  std::cout << "Initialized No Regressor: "
            << fmt::format("nb_weeks = {}, nb_groups = {}, neighbor = {}, "
                           "constant_lam = {}\n",
                           nb_weeks, groups.size(), neighbor_factor,
                           constant_lambdas);
}

GeneratorNoRegressor::GeneratorNoRegressor(std::string calls_path,
                                           std::string neighbors_path,
                                           std::string info_path) {
  full_neighbors = true;
  neighbor_factor = 1;
  read_covariates = false;
  use_simulation = false;
  constant_lambdas = false;
  auto info_arq = ifstream(info_path, ios::in);
  info_arq >> T >> D >> R >> C >> nb_regressors >> nb_holidays_years;
  slot_duration = 24 / T;
  fmt::print("info: {} {} {} {} {} {}\n", T, D, R, C, nb_regressors,
             nb_holidays_years);
  daily_obs = std::vector<int>(D, 0);
  for (int d = 0; d < D; ++d) {
    info_arq >> daily_obs[d];
  }
  fmt::print("daily obs: {}\n", daily_obs);
  info_arq.close();
  nb_land_types = nb_regressors - 1;
  // TODO: remove comments
  //  nb_land_types = nb_regressors - 2;
  //  nb_regressors = 1 + nb_land_types;
  unsigned long max_obs = *max_element(daily_obs.begin(), daily_obs.end());
  unsigned long min_obs = *min_element(daily_obs.begin(), daily_obs.end());
  nb_weeks = max_obs;
  xt::xarray<int> nb_observations_file = xt::zeros<int>({C, D, T, R});
  xt::xarray<int> sample_file = xt::zeros<int>({C, D, T, R, max_obs});
  xt::xarray<int> nb_arrivals_file = xt::zeros<int>({C, D, T, R});

  auto calls_arq = ifstream(calls_path, ios::in);
  std::string aux_str;
  // is_holidays = std::vector<std::pair<bool, int>>(max_obs,
  // make_pair(false,-1));
  int count_line = 0;
  do {
    std::getline(calls_arq, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    int t, d, r, c, j, h, val;
    // TODO: remove comment
    ss >> t >> d >> r >> c >> j >> val; // >> h;
    // fmt::print("{} {} {} {} {} {} {}\n", t,d,r,c,j,val,h);
    sample_file(c, d, t, r, j) = val;
    nb_arrivals_file(c, d, t, r) += val;
    // nb_observations_file(c,d,t,r) += 1;
    ++count_line;
    // if(count_line % 45 == 0){
    // 	cin.get();
    // }
  } while (true);
  fmt::print("calls.dat {} lines read\n", count_line);
  // cin.get();
  calls_arq.close();

  type_region = std::vector<int>(R, -1);
  regressors = xt::zeros<double>({nb_regressors, R});
  neighbors = std::vector<vector<int>>(R, std::vector<int>());
  distance = xt::zeros<double>({R, R});
  regions = std::vector<Location>(R, null_location);
  auto neighbors_arq = ifstream(neighbors_path, ios::in);
  double pop1, pop2;
  count_line = 0;
  while (true) {
    int ind, terrain_type, s;
    double lat, longi, dist;
    std::getline(neighbors_arq, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    // TODO: Remove comments
    ss >> ind >> lat >> longi; // >> terrain_type;
    type_region[ind] = 0;
    regions[ind] = make_pair(lat, longi);
    // TODO: nb_landtypes
    for (int j = 0; j < nb_regressors; ++j) {
      ss >> regressors(j, ind);
    }
    // TODO: Remove comments
    // ss >> pop1 >> pop2;
    // regressors(nb_regressors - 2, ind) = pop1;
    // regressors(nb_regressors - 1, ind) = pop1 + pop2;
    while (ss >> s >> dist) {
      distance(ind, s) = dist;
      neighbors[ind].push_back(s);
    }
    ++count_line;
    // if(count_line % 45 == 0){
    // 	cin.get();
    // }
  }
  fmt::print("neighbors.dat {} lines read\n", count_line);
  neighbors_arq.close();

  sample = xt::zeros<int>({7 * T, R, C, min_obs});
  nb_observations = xt::zeros<int>({C, R, 7 * T});
  nb_arrivals = xt::zeros<int>({C, R, 7 * T});
  durations = vector<double>(7 * T, 0.5);
  estimated = xt::zeros<double>({C, R, 7 * T});
  theoretical_lambda = xt::zeros<double>({7 * T, R, C});
  for (int r = 0; r < R; ++r) {
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          int index = d * T + t;
          nb_observations(c, r, index) = daily_obs[d];
          nb_arrivals(c, r, index) = nb_arrivals_file(c, d, t, r);
          for (int j = 0; j < min_obs; ++j) {
            sample(index, r, c, j) = sample_file(c, d, t, r, j);
          }
          estimated(c, r, index) =
              static_cast<double>(nb_arrivals(c, r, index)) /
              (nb_observations(c, r, index) * durations[t]);
          // fmt::print("c{} r{} t{}, emp = {}\n", c, r, d * T + t,
          //            estimated(c, d, t, r));
        }
      }
    }
  }

  // for (int c = 0; c < C; ++c) {
  //   for (int d = 0; d < D; ++d) {
  //     for (int t = 0; t < T; ++t) {
  //       for (int r = 0; r < R; ++r) {
  //         fmt::print("c{} r{} t{}, emp = {}\n", c, r, d * T + t,
  //                    estimated(c, r, d * T + t));
  //       }
  //     }
  //   }
  // }

  // fmt::print("END EMP\n");
  // cin.get();

  // for (int r = 0; r < R; ++r) {
  //   fmt::print("r{} (t{}): ", r, type_region[r]);
  //   for (auto s : neighbors[r]) {
  //     fmt::print("({},{:.1f}) ", s, distance(r, s));
  //   }
  //   fmt::print("\n");
  // }
  // fmt::print("END Neighbors\n");
  // cin.get();

  T = 7 * T;
  ulong DT = T;

  which_group = vector<int>(T, 0);
  groups = vector<vector<int>>(T, vector<int>());
  for (int t = 0; t < T; ++t) {
    groups[t].push_back(t);
    which_group[t] = t;
    // fmt::print("which_group {} = {}\n", t + 1, which_group[t] + 1);
  }
  // fmt::print("END Groups\n");
  // cin.get();

  // for (int c = 0; c < C; ++c) {
  //   for (int r = 0; r < R; ++r) {
  //     for (int t = 0; t < T; ++t) {
  //       fmt::print("c{}, r{} , t{}, arr = {}\n", c, r, t, nb_arrivals(c, r,
  //       t));
  //     }
  //   }
  // }

  // fmt::print("End sample/arr/obs {} {} {}\n", C, R, T);
  // cin.get();
  for (int r = 0; r < R; ++r) {
  }
  sigma = 0.5;
  beta_bar = 1;
  max_iter = 30;
  weight = 0.1;
  alpha = 0.1 * xt::ones<double>({R, R});
  durations = vector<double>(T, 0.5);
  std::cout << "Initialized No Regressor Real Data\n";
}

GeneratorNoRegressor::GeneratorNoRegressor(
    xt::xarray<int> &N, xt::xarray<int> &M, std::vector<double> &a_durations,
    std::vector<std::vector<int>> &a_groups, std::vector<double> &a_weights,
    xt::xarray<double> &alphas, xt::xarray<double> &a_distance,
    std::vector<int> &a_type, std::vector<std::vector<int>> &a_neighbors) {
  if (N.dimension() != 3) { // N should be C,R,T
    fmt::print(
        "Error: N has {} dimensions but should be 3. Problem was not set.\n",
        N.dimension());
    exit(1);
  }
  if (M.dimension() != 3) { // M also should be C,R,T
    fmt::print(
        "Error: M has {} dimensions but should be 3. Problem was not set.\n",
        M.dimension());
    exit(1);
  }

  full_neighbors = true;
  neighbor_factor = 1;
  read_covariates = false;
  use_simulation = false;

  C = N.shape(0);
  R = N.shape(1);
  T = N.shape(2);

  nb_observations = N;
  nb_arrivals = M;
  durations = a_durations;
  alpha = alphas;
  weights = a_weights;
  distance = a_distance;
  neighbors = a_neighbors;
  type_region = a_type;
  groups = a_groups;
  int nb_groups = groups.size();

  which_group = vector<int>(T, 0);
  for (int g = 0; g < nb_groups; ++g) {
    for (auto i : groups[g]) {
      which_group[i] = g;
    }
  }

  g_params.EPS = pow(10, -3);
  sigma = 0.5;
  beta_bar = 1;
  max_iter = 30;
}

xt::xarray<double> laspated_no_reg(
    xt::xarray<int> &N, xt::xarray<int> &M, std::vector<double> &a_durations,
    std::vector<std::vector<int>> &a_groups, std::vector<double> &a_weights,
    xt::xarray<double> &alphas, xt::xarray<double> &a_distance,
    std::vector<int> &a_type, std::vector<std::vector<int>> &a_neighbors,
    xt::xarray<double> &x) {

  GeneratorNoRegressor gen(N, M, a_durations, a_groups, a_weights, alphas,
                           a_distance, a_type, a_neighbors);

  if (x.dimension() != 3) {
    fmt::print("Error: x has {} dimensions but must be 3.\n", x.dimension());
    exit(1);
  }

  if (x.shape(0) != gen.C || x.shape(1) != gen.R || x.shape(2) != gen.T) {
    std::string shape_x =
        fmt::format("({},{},{})", x.shape(0), x.shape(1), x.shape(2));
    fmt::print("Error: x has shape {}, but expected is ({},{},{}).\n", shape_x,
               gen.C, gen.R, gen.T);
    exit(1);
  }

  auto lambda = x;
  if (g_params.method == "calibration") {
    auto f_val = gen.projected_gradient_armijo_feasible(lambda);
  } else if (g_params.method == "cross_validation") {
    auto result = gen.cross_validation(
        g_params.cv_proportion, g_params.weights_list, g_params.weights_list);
    lambda = result.lambda;
    fmt::print("Cross validation best weight = {}\n", result.weight);
  }
  return lambda;
}

void GeneratorNoRegressor::test() {
  double epsilon = g_params.EPS;
  fmt::print("ENTERED test\n");
  if (g_params.generator_folder != "") {
    // vector<double> test_weights = {0,   0.2, 0.4, 0.6, 0.8, 1.0,
    //                                1.2, 1.4, 1.6, 1.8, 2};
    vector<double> test_weights = {0,      0.0001, 0.0002, 0.0003,
                                   0.0005, 0.0007, 0.0009, 0.001};
    for (size_t i = 0; i < test_weights.size(); ++i) {
      test_weights[i] *= 8;
    }
    vector<double> test_alphas = test_weights;
    xt::xarray<double> x = epsilon * xt::ones<double>({C, R, T});
    auto result = cross_validation(0.2, test_alphas, test_weights);
    double best_w = result.weight;
    x = epsilon * xt::ones<double>({C, R, T});
    weights = vector<double>(groups.size(), best_w);
    alpha = best_w * 1 * xt::ones<double>({R, R});

    vector<double> f_val;
    if (g_params.type_proj_gradient == 2) {
      f_val = projected_gradient_armijo_feasible(x);
    } else {
      f_val = projected_gradient_old(x);
    }

    // for(int r = 0; r < R; ++r){
    // 	fmt::print("r = {}, regressors = {:.3f} {:.3f} {:.3f}\n", r,
    // regressors(0,r), 		regressors(1,r), regressors(2,r));
    // }

    string rates_no_reg_filename =
        fmt::format("err_no_reg/rates_no_reg_real{}_g336_obs105_n1.txt", R);
    ofstream rates_no_reg_real(rates_no_reg_filename, std::ios::out);
    for (int c = 0; c < C; ++c) {
      for (int t = 0; t < T; ++t) {
        double sum = 0;
        double sum_est = 0;
        for (int r = 0; r < R; ++r) {
          sum += x(c, r, t);
          sum_est += estimated(c, r, t);
          // fmt::print("r = {}, rate = {}, est = {}\n", r,  x(c, r, t),
          // estimated(c,r,t));
        }
        // fmt::print("Sum = {}\n", sum);
        rates_no_reg_real << fmt::format("{} {} {} {}\n", c, t, sum / 0.5,
                                         sum_est);
        fmt::print("{} {} {} {}\n", c, t, sum / 0.5, sum_est);
      }
    }
    rates_no_reg_real.close();
    fmt::print("Wrote time rates at {}\n", rates_no_reg_filename);
    string region_rates =
        fmt::format("region_rates_no_reg_real{}_g336_obs105_n1.txt", R);
    ofstream arq_regions(region_rates, std::ios::out);
    for (int r = 0; r < R; ++r) {
      vector<double> sums(C, 0.0);
      vector<double> sum_ests(C, 0.0);
      for (int c = 0; c < C; ++c) {
        for (int t = 0; t < T; ++t) {
          sums[c] += x(c, r, t);
          sum_ests[c] += estimated(c, r, t);
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
    fmt::print("Wrote region rates at {}\n", rates_no_reg_filename);
    fmt::print("best_w = {}\n", best_w);
    exit(0);
  }

  vector<double> test_weights;
  test_weights =
      (constant_lambdas)
          ? vector<double>({0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2})
          : vector<double>(
                {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1});
  if (nb_weeks == 20) {
    for (size_t i = 0; i < test_weights.size(); ++i) {
      test_weights[i] /= 10.0;
    }
  } else if (nb_weeks == 100) {
    for (size_t i = 0; i < test_weights.size(); ++i) {
      test_weights[i] /= 50.0;
    }
  } else if (nb_weeks == 1000) {
    for (size_t i = 0; i < test_weights.size(); ++i) {
      test_weights[i] /= 500.0;
    }
  }
  vector<double> test_alphas = test_weights;
  xt::xarray<double> x = xt::ones<double>({C, R, T});
  string ins_type = (read_covariates) ? "_cov" : "";
  string type_neighborhood = (full_neighbors) ? "n1" : "n2";
  string type_alpha = (neighbor_factor == 0) ? "_a0" : "";
  string real = (g_params.generator_folder == "") ? "" : "_real";
  string err_arq_name =
      fmt::format("err_no_reg/err{}_no_reg{}_g{}_obs{}_{}{}.txt", ins_type,
                  real, groups.size(), nb_weeks, type_neighborhood, type_alpha);
  ofstream err_arq(err_arq_name, std::ios::out);

  std::random_device rd;
  std::mt19937 gen(rd());

  int num_trials = 1;
  vector<vector<double>> err_by_replications(
      num_trials, vector<double>(test_weights.size(), 0.0));

  vector<pair<double, double>> cv_by_replications(num_trials,
                                                  make_pair(GRB_INFINITY, -1));
  for (int k = 0; k < num_trials; ++k) {
    fmt::print("trial {}\n", k);
    // for (int t = 0; t < T; ++t) {
    //   for (int r = 0; r < R; ++r) {
    //     int x = r % 10;
    //     int y = r / 10;
    //     double cx = x + 0.5;
    //     double cy = y + 0.5;
    //     for (int c = 0; c < C; ++c) {
    //       if ((!is_red[r] && which_group[t] % 2 == 0) ||
    //           (is_red[r] && which_group[t] % 2 != 0)) {
    //         theoretical_lambda(t, r, c) = (constant_lambdas) ? 0.1 : (cx +
    //         cy);
    //       } else if ((!is_red[r] && which_group[t] % 2 != 0) ||
    //                  (is_red[r] && which_group[t] % 2 == 0)) {
    //         theoretical_lambda(t, r, c) =
    //             (constant_lambdas) ? 0.5 : 5 * (cx + cy);
    //       } else {
    //         fmt::print("ERROR: Unpredictable lambda case at t = {} and r = {}
    //         "
    //                    "({},{})\n",
    //                    t, r, cx, cy);
    //         exit(1);
    //       }
    //       nb_arrivals(c, r, t) = 0;
    //       for (int n = 0; n < nb_observations_total; ++n) {
    //         poisson_distribution<int> pd(theoretical_lambda(t, r, c) *
    //                                      durations[t]);
    //         int this_nb_arrival = pd(gen);
    //         sample(t, r, c, n) = this_nb_arrival;
    //         nb_arrivals(c, r, t) += this_nb_arrival;
    //       }
    //     }
    //   }
    // }
    for (int i = 0; i < test_weights.size(); ++i) {
      double avg_err = 0;
      x = epsilon * xt::ones<double>({C, R, T});
      alpha = test_alphas[i] * neighbor_factor * xt::ones<double>({R, R});
      // alpha = 0;
      weights = vector<double>(groups.size(), test_weights[i]);
      // weights = vector<double>(groups.size(), 0);
      auto f_val = projected_gradient_armijo_feasible(x);
      double err = average_difference(x);
      err_by_replications[k][i] = err;
      fmt::print("\tw = {}, err = {}\n", test_weights[i], err);
    }

    test_alphas = test_weights;
    x = epsilon * xt::ones<double>({C, R, T});
    auto result = cross_validation(0.2, test_alphas, test_weights);
    double best_w = result.weight;
    x = epsilon * xt::ones<double>({C, R, T});
    weights = vector<double>(groups.size(), best_w);
    auto f_val = projected_gradient_armijo_feasible(x);
    cv_by_replications[k] = make_pair(average_difference(x), best_w);
  }

  for (int i = 0; i < test_weights.size(); ++i) {
    double avg = 0.0;
    for (int k = 0; k < num_trials; ++k) {
      avg += err_by_replications[k][i];
    }
    // avg /= (num_trials - 1);
    avg /= (num_trials);
    fmt::print("w = {}, avg = {}\n", test_weights[i], avg);
    // double std_dev = 0.0;
    // for (int k = 0; k < num_trials; ++k) {
    //   std_dev += pow(err_by_replications[k][i] - avg, 2);
    // }
    // std_dev /= (num_trials - 1);

    // fmt::print("Weight {}: mean {:.6f} dev {:.6f} std_err {:.6f}\n",
    //            test_weights[i], avg, std_dev, std_dev / sqrt(num_trials));

    err_arq << avg << "\n";
  }

  err_arq.close();
  double avg = 0.0;
  for (int k = 0; k < num_trials; ++k) {
    avg += cv_by_replications[k].first;
  }
  avg /= (num_trials - 1);
  double std_dev = 0.0;
  for (int k = 0; k < num_trials; ++k) {
    std_dev += pow(avg - cv_by_replications[k].first, 2);
  }
  std_dev /= (num_trials - 1);
  fmt::print("CV: mean {:.6f} dev {:.1f} std_err {:.6f}\n", avg, std_dev,
             std_dev / sqrt(num_trials));
}

void GeneratorNoRegressor::calibrate() {
  double epsilon = g_params.EPS;
  vector<double> test_weights = g_params.weights_list;
  vector<double> alphas = test_weights;
  xt::xarray<double> x = xt::ones<double>({C, R, T});
  fmt::print("Running projected gradient for weights {}\n", test_weights);
  string ins_type = (read_covariates) ? "_cov" : "";
  if (read_covariates) {
    ins_type = "_cov";
  } else if (!use_simulation) {
    ins_type = "_lam";
  }
  string type_neighborhood = (full_neighbors) ? "n1" : "n2";
  string type_alpha = (neighbor_factor == 0) ? "_a0" : "";
  string holidays = (T == 4 * 15) ? "_h1" : "";
  string err_name =
      fmt::format("err{}_no_reg_g{}_obs{}_{}{}{}.txt", ins_type, groups.size(),
                  nb_weeks, type_neighborhood, type_alpha, holidays);
  ofstream err_arq(err_name, std::ios::out);
  ofstream week_arq(fmt::format("week_no_reg.txt"), std::ios::out);
  double min_err = 10e100;
  int min_w = -1;
  for (int i = 0; i < test_weights.size(); ++i) {
    // x = epsilon*xt::ones<double>({C,R,T});
    x = 0.001 * xt::ones<double>({C, R, T});
    alpha = alphas[i] * neighbor_factor * xt::ones<double>({R, R});
    // alpha = xt::zeros<double>({R,R});
    // alpha = 0;
    weights = vector<double>(groups.size(), test_weights[i]);
    // weights = vector<double>(groups.size(), 0);
    auto f_val = projected_gradient_armijo_feasible(x);
    // fmt::print("alpha = {}, weight = {}, diff = {}\n",alpha, weights[0],
    // average_difference(x));

    if (g_params.generator_folder == "") {
      double err = average_difference(x);
      if (err < min_err) {
        min_err = err;
        min_w = i;
      }
      err_arq << fmt::format("{:.7f}", err) << "\n";
      fmt::print("w({}) = {}, err = {}\n", i, test_weights[i], err);
    }
  }
  alpha = alphas[min_w];
  if (g_params.generator_folder == "") {
    fmt::print("Best weight = {} {}\n", alpha, min_w);
    fmt::print("Wrote weight results at {}\n", err_name);
  }

  weights = vector<double>(groups.size(), test_weights[min_w]);
  x = epsilon * xt::ones<double>({C, R, T});
  auto f_val = projected_gradient_armijo_feasible(x);

  ofstream x_arq(fmt::format("x_no_reg.txt"), std::ios::out);
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        x_arq << fmt::format("{} {} {} = {}\n", c, r, t, x(c, r, t));
      }
    }
  }
  x_arq.close();
  fmt::print("Wrote intensities at x_no_reg.txt\n");
}

void GeneratorNoRegressor::write_cv_results(CrossValidationResult &cv_result) {
  string ins_type = (read_covariates) ? "_cov" : "";
  string type_neighborhood = (full_neighbors) ? "n1" : "n2";
  string type_alpha = (neighbor_factor == 0) ? "_a0" : "";
  string holidays = (T == 4 * 15) ? "_h1" : "";
  string err_name =
      fmt::format("cv_{}_no_reg_g{}_obs{}_{}{}{}.txt", ins_type, groups.size(),
                  nb_weeks, type_neighborhood, type_alpha, holidays);
  ofstream arq(err_name, std::ios::out);
  arq << "Best weight = " << cv_result.weight
      << ", cpu time (s) = " << cv_result.cpu_time << "\n";
  xt::xarray<double> x =
      g_params.EPS * xt::ones<double>(cv_result.lambda.shape());
  auto f_val = projected_gradient_armijo_feasible(x);
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        arq << c << " " << r << " " << t << " " << x(c, r, t) << "\n";
      }
    }
  }
  arq.close();
  fmt::print("Wrote intensities at cv_x_no_reg.txt\n");
  ofstream weekly_total("weekly_total_no_reg.txt", std::ios::out);
  for (int t = 0; t < T; ++t) {
    double sum = 0;
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        sum += x(c, r, t);
      }
    }
    weekly_total << sum << "\n";
  }
  weekly_total.close();

  ofstream weekly0("weekly_0_no_reg.txt", std::ios::out);
  for (int t = 0; t < T; ++t) {
    double sum = 0;
    for (int r = 0; r < R; ++r) {
      sum += x(0, r, t);
    }
    weekly0 << sum << "\n";
  }
  weekly0.close();

  ofstream weekly1("weekly_1_no_reg.txt", std::ios::out);
  for (int t = 0; t < T; ++t) {
    double sum = 0;
    for (int r = 0; r < R; ++r) {
      sum += x(1, r, t);
    }
    weekly1 << sum << "\n";
  }
  weekly1.close();

  ofstream weekly2("weekly_2_no_reg.txt", std::ios::out);
  for (int t = 0; t < T; ++t) {
    double sum = 0;
    for (int r = 0; r < R; ++r) {
      sum += x(1, r, t);
    }
    weekly2 << sum << "\n";
  }
  weekly2.close();
}

std::vector<double>
GeneratorNoRegressor::projected_gradient_old(xt::xarray<double> &x) {
  int k = 0;
  std::vector<double> f_val;
  double b_param = 2.0;
  double beta_k = b_param;

  double accuracy = 0.01;
  double eps = g_params.EPS;
  double upper_lambda = 1e3;
  double upper_bound = GRB_INFINITY;
  max_iter = g_params.max_iter;
  int j = 0;
  comp_wise_max(x, eps);
  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());
  // fmt::print("MAX ITER = {}\n", max_iter);

  while (k < max_iter) {
    double fold = oracle_objective_model_new(x);
    xt::xarray<double> gradient = oracle_gradient_model_new(x);
    z = x - beta_k * gradient;
    comp_wise_max(z, eps);
    bool stop = false;
    j = 0;
    diff_aux = x - z;
    double rhs = mat_prod(gradient, diff_aux);
    double f = GRB_INFINITY;
    fmt::print("k = {}, fold = {}, rhs = {}\n", k, fold, rhs);
    while (!stop) {
      z_aux = x + (1 / pow(2, j)) * (z - x); // z_aux expected to tend to x
      f = oracle_objective_model_new(z_aux);
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
    // printf("k = %d, fold = %.8f, f = %.8f, rhs = %f, j = %d\n", k, fold,
    //        f / pow(10, 6), rhs, j);

    f_val.push_back(f);
    x = z_aux;
    ++k;
    beta_k = b_param / pow(2, j);
  }
  cin.get();
  return f_val;
}

std::vector<double> GeneratorNoRegressor::projected_gradient_armijo_feasible(
    xt::xarray<double> &x) {
  int k = 0;
  std::vector<double> f_val;
  double b_param = 2.0;
  double beta_k = b_param;
  f_val.reserve(max_iter);

  double accuracy = 0.01;
  double upper_bound = GRB_INFINITY;
  double upper_lambda = 100;
  double lower_bound = -GRB_INFINITY;

  double eps = g_params.EPS; // lower bound decision variables
  max_iter = g_params.max_iter;
  comp_wise_max(x, eps);

  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  int j = 0;
  double fold = oracle_objective_model_new(x);
  xt::xarray<double> gradient = oracle_gradient_model_new(x);
  while (k < max_iter) {
    z = x - beta_k * gradient;
    comp_wise_max(z, eps);
    diff_aux = x - z;
    double rhs = mat_prod(gradient, diff_aux);
    double f = oracle_objective_model_new(z);
    // fmt::print("k = {}, fold = {}, f(z) = {}, rhs = {}\n", k, fold, f, rhs);
    if (rhs > 0.0 && f > fold - (sigma)*rhs) {
      bool stop = false;
      z_aux = xt::zeros<double>(x.shape());
      double this_pow = 1.0;
      int count = 0;
      double best_pow = -1.0;
      double best_val = f;
      while (!stop) {
        z_aux = x + (1.0 / this_pow) * (z - x);
        f = oracle_objective_model_new(z_aux);
        // fmt::print("\tj = {}, f = {}\n", j, f);
        if (f < best_val) {
          best_val = f;
          best_pow = this_pow;
        }
        if (f <= fold - (sigma / this_pow) * rhs) {
          stop = true;
        } else {
          this_pow *= 2;
        }
        ++j;
      }
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
      // print_var(x,fmt::format("final x iter {}", k));
      // fmt::print("k = {}\n", k+1);
      // beta_k = b_param / pow(2, j); // current
      // beta_k = beta_k*2 / this_pow; //suggestion
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
    fold = f;
    gradient = oracle_gradient_model_new(x);
    upper_bound = min(upper_bound, f);
    lower_bound =
        max(lower_bound, f + get_lower_bound(x, gradient, eps, upper_lambda));
    if (upper_bound - lower_bound < accuracy) {
      break;
    }
    // fmt::print("Final this_pow = {}, beta_k = {}\n", this_pow, beta_k);
    ++k;
    // cin.get();
  }
  // cin.get();
  return f_val;
}

double GeneratorNoRegressor::mat_prod(xt::xarray<double> &a,
                                      xt::xarray<double> &b) {
  double sum = 0;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        sum += a(c, r, t) * b(c, r, t);
      }
    }
  }
  return sum;
}

double GeneratorNoRegressor::get_lower_bound(xt::xarray<double> &x,
                                             xt::xarray<double> &grad,
                                             double lower_lambda,
                                             double upper_lambda) {
  double change = 0.0;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        if (grad(c, r, t) < 0.0) {
          change += grad(c, r, t) * (upper_lambda - x(c, r, t));
        } else if (grad(c, r, t) > 0.0) {
          change += grad(c, r, t) * (lower_lambda - x(c, r, t));
        }
      }
    }
  }
  return change;
}

xt::xarray<double>
GeneratorNoRegressor::oracle_gradient_model_new(xt::xarray<double> &x) {
  xt::xarray<double> gradient = xt::zeros<double>(x.shape());
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double current_lambda = x(c, r, t);
        double grad_component = nb_observations(c, r, t) * durations[t] -
                                (nb_arrivals(c, r, t) / current_lambda);
        double pre_n = grad_component;
        for (int s : neighbors[r]) {
          if (type_region[r] == type_region[s]) {
            grad_component += 2 * alpha(r, s) * nb_observations(c, r, t) *
                              nb_observations(c, s, t) *
                              (x(c, r, t) - x(c, s, t)) / (distance(r, s));
          }
        }
        double post_n = grad_component - pre_n;
        gradient(c, r, t) = grad_component;
      }
    }
  }
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double sum_groups = 0;
        auto &group = groups[which_group[t]];
        for (int j = 0; j < group.size(); ++j) {
          int tp = group[j];
          if (tp != t) {
            gradient(c, r, t) +=
                2 * weights[which_group[t]] * nb_observations(c, r, t) *
                nb_observations(c, r, tp) * (x(c, r, t) - x(c, r, tp));
          }
        }
      }
    }
  }

  return gradient;
}

double GeneratorNoRegressor::oracle_objective_model_new(xt::xarray<double> &x) {
  double f = 0;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double current_lambda = x(c, r, t);
        f += nb_observations(c, r, t) * current_lambda * durations[t] -
             nb_arrivals(c, r, t) * log(current_lambda * durations[t]);
        for (int s : neighbors[r]) {
          if (type_region[r] == type_region[s]) {
            f += (0.5 * alpha(r, s)) * nb_observations(c, r, t) *
                 nb_observations(c, s, t) * pow(x(c, r, t) - x(c, s, t), 2) /
                 distance(r, s);
          }
          // f += (0.5*alpha(r,s))*pow(x(c,r,t)- x(c,s,t), 2) / distance(r,s);
        }
      }
    }
  }
  // sum_{c \in C}sum{r \in R}sum_{t \in T}\sum_{t' \neq t \in which_group[t]}
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int grindex = 0; grindex < groups.size(); ++grindex) {
        auto &group = groups[grindex];
        for (int j = 0; j < group.size(); ++j) {
          int t = group[j];
          for (int itp = 0; itp < group.size(); ++itp) {
            int tp = group[itp];
            if (tp != t) {
              f += (0.5 * weights[grindex]) * nb_observations(c, r, t) *
                   nb_observations(c, r, tp) *
                   (pow(x(c, r, t) - x(c, r, tp), 2));
            }
          }
        }
      }
    }
  }

  return f;
}

xt::xarray<double>
GeneratorNoRegressor::oracle_gradient_model(xt::xarray<double> &x) {

  xt::xarray<double> gradient = xt::zeros<double>(x.shape());

  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double current_lambda = x(c, r, t);
        double grad_component = nb_observations(c, r, t) * durations[t] -
                                (nb_arrivals(c, r, t) / current_lambda);
        // double prev_comp = grad_component;
        // double sum_neighbors = 0;
        for (int s : neighbors[r]) {
          if (!full_neighbors || (type_region[r] == type_region[s])) {
            grad_component +=
                2 * alpha(r, s) * (x(c, r, t) - x(c, s, t)) / (distance(r, s));
          }
        }
        gradient(c, r, t) = grad_component;
      }
    }
  }

  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double sum_groups = 0;
        auto &group = groups[which_group[t]];
        for (int j = 0; j < group.size(); ++j) {
          int tp = group[j];
          if (tp != t) {
            gradient(c, r, t) +=
                2 * weights[which_group[t]] * (x(c, r, t) - x(c, r, tp));
          }
        }
      }
    }
  }

  return gradient;
}

double GeneratorNoRegressor::oracle_objective_model(xt::xarray<double> &x) {
  double f = 0;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        double current_lambda = x(c, r, t);
        f += nb_observations(c, r, t) * current_lambda * durations[t] -
             nb_arrivals(c, r, t) * log(current_lambda * durations[t]);
        for (int s : neighbors[r]) {
          if (!full_neighbors || (type_region[r] == type_region[s])) {
            f += (0.5 * alpha(r, s)) * pow(x(c, r, t) - x(c, s, t), 2) /
                 distance(r, s);
          }
          // f += (0.5*alpha(r,s))*pow(x(c,r,t)- x(c,s,t), 2) / distance(r,s);
        }
      }
    }
  }
  // sum_{c \in C}sum{r \in R}sum_{t \in T}\sum_{t' \neq t \in which_group[t]}
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int grindex = 0; grindex < groups.size(); ++grindex) {
        auto &group = groups[grindex];
        for (int j = 0; j < group.size(); ++j) {
          int t = group[j];
          for (int itp = 0; itp < group.size(); ++itp) {
            int tp = group[itp];
            if (tp != t) {
              f +=
                  (0.5 * weights[grindex]) * (pow(x(c, r, t) - x(c, r, tp), 2));
            }
          }
        }
      }
    }
  }

  return f;
}

double GeneratorNoRegressor::average_difference(xt::xarray<double> &x) {
  vector<double> difference_l2;
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        difference_l2.push_back(abs(theoretical_lambda(t, r, c) - x(c, r, t)) /
                                theoretical_lambda(t, r, c));
      }
    }
  }
  // fmt::print("Error: {}\n", difference_l2);
  return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) /
         difference_l2.size();
}

void GeneratorNoRegressor::comp_wise_max(xt::xarray<double> &z, double eps) {
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        z(c, r, t) = max(z(c, r, t), eps);
      }
    }
  }
}

bool GeneratorNoRegressor::is_neighbor(int r, int s) { return r != s; }

void GeneratorNoRegressor::print_var(xt::xarray<double> &x,
                                     std::string prefix) {
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        if (abs(x(c, r, t)) > g_params.EPS) {
          fmt::print("nonzero: c = {}, r = {}, t = {} (type_region = {}), val  "
                     "= {}, theoretical = {}\n",
                     c + 1, r + 1, t + 1, type_region[r] + 1, x(c, r, t),
                     theoretical_lambda(t, r, c));
        }
      }
    }
  }
  // fmt::print("{} c = 0, r = 84, t = 27, val = {}\n", prefix, x(0,84,27));
  // fmt::print("{} c = 0, r = 85, t = 27, val = {}\n", prefix, x(0,85,27));
  cin.get();
}

CrossValidationResult GeneratorNoRegressor::cross_validation(
    double proportion, vector<double> &alphas, vector<double> &group_weights) {
  auto t0 = std::chrono::high_resolution_clock::now();
  int nb_observations_total = sample.shape(3);
  int nb_groups = groups.size();
  double min_loss = GRB_INFINITY;
  double cpu_time = 0;
  int nb_in_block = floor(nb_observations_total * proportion);
  xt::xarray<int> initial_nb_obs = nb_observations;
  xt::xarray<int> initial_nb_arrivals = nb_arrivals;
  // fmt::print("cross validation weights.size = {}\n", group_weights.size());
  double best_alpha = GRB_INFINITY;
  double best_weight = GRB_INFINITY;
  // fmt::print("Running cross validation with proportion = {} and weights =
  // {}\n", proportion, group_weights);
  for (int index_alpha = 0; index_alpha < alphas.size(); ++index_alpha) {
    double likelihood = 0;
    alpha = alphas[index_alpha] * neighbor_factor * xt::ones<double>({R, R});
    weights = vector<double>(groups.size(), group_weights[index_alpha]);
    // fmt::print("Testing weight = {}\n", group_weights[index_alpha]);
    for (int index_cross = 0; index_cross < floor(1 / proportion);
         ++index_cross) {
      xt::xarray<int> nb_observations_current = xt::zeros<int>({C, R, T});
      xt::xarray<int> nb_calls_current = xt::zeros<int>({C, R, T});
      for (int index = index_cross * nb_in_block;
           index < (index_cross + 1) * nb_in_block; ++index) {
        for (int c = 0; c < C; ++c) {
          for (int r = 0; r < R; ++r) {
            for (int t = 0; t < T; ++t) {
              ++nb_observations_current(c, r, t);
              nb_calls_current(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      xt::xarray<double> x = g_params.EPS * xt::ones<double>({C, R, T});
      nb_observations = nb_observations_current;
      nb_arrivals = nb_calls_current;
      vector<double> f_val;
      if (g_params.type_proj_gradient == 2) {
        f_val = projected_gradient_armijo_feasible(x);
      } else if (g_params.type_proj_gradient == 1) {
        f_val = projected_gradient_old(x);
      } else {
        printf("Invalid type_proj_gradient = %d\n",
               g_params.type_proj_gradient);
        exit(1);
      }
      xt::xarray<int> nb_calls_remaining = xt::zeros<int>({C, R, T});
      for (int index = 0; index < index_cross * nb_in_block; ++index) {
        for (int c = 0; c < C; ++c) {
          for (int r = 0; r < R; ++r) {
            for (int t = 0; t < T; ++t) {
              nb_calls_remaining(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      for (int index = (index_cross + 1) * nb_in_block;
           index < nb_observations_total; ++index) {
        for (int c = 0; c < C; ++c) {
          for (int r = 0; r < R; ++r) {
            for (int t = 0; t < T; ++t) {
              nb_calls_remaining(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      double f = 0;
      for (int c = 0; c < C; ++c) {
        for (int r = 0; r < R; ++r) {
          for (int t = 0; t < T; ++t) {
            double current_lambda = x(c, r, t);
            f += (nb_observations_total - nb_in_block) * current_lambda *
                     durations[t] -
                 nb_calls_remaining(c, r, t) * log(current_lambda);
          }
        }
      }
      // fmt::print("\tf test set {} = {}, first f_val = {:.1f}, last f_val =
      // {:.1f}\n",index_cross, f, f_val[0], f_val.back());
      likelihood += f;
    }
    likelihood = likelihood / floor(1 / proportion);
    // fmt::print("Likelihood current_alpha {} = {}\n", alphas[index_alpha],
    // likelihood);
    if (likelihood < min_loss) {
      min_loss = likelihood;
      best_alpha = alphas[index_alpha];
      best_weight = group_weights[index_alpha];
    }
  }
  alpha = best_alpha;
  weights = vector<double>(groups.size(), best_weight);
  xt::xarray<int> nb_observations_current = xt::zeros<int>({C, R, T});
  xt::xarray<int> nb_calls_current = xt::zeros<int>({C, R, T});
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        for (int index = 0; index < nb_observations_total; ++index) {
          nb_calls_current(c, r, t) += sample(t, r, c, index);
          ++nb_observations_current(c, r, t);
        }
      }
    }
  }

  nb_observations = nb_observations_current;
  nb_arrivals = nb_calls_current;
  xt::xarray<double> x = g_params.EPS * xt::ones<double>({C, R, T});
  vector<double> f_val;
  if (g_params.type_proj_gradient == 2) {
    f_val = projected_gradient_armijo_feasible(x);
  } else if (g_params.type_proj_gradient == 1) {
    f_val = projected_gradient_old(x);
  } else {
    fmt::print("Invalid type_proj_gradient = {} at line {}\n",
               g_params.type_proj_gradient, __LINE__);
    exit(1);
  }
  auto dt = std::chrono::high_resolution_clock::now();
  cpu_time = std::chrono::duration_cast<std::chrono::seconds>(dt - t0).count();
  nb_observations = initial_nb_obs;
  nb_arrivals = initial_nb_arrivals;
  // fmt::print("best_weight = {}\n", best_weight);
  return {cpu_time, best_weight, x};
}

// // theoretical_lambda = xt::zeros<double>({T,R,C});
// xt::xarray<double> integrated_lambda = xt::zeros<double>({T,R});

// for(int r = 1; r <= R; ++r){
// 	int i = r % 10;
// 	int j = floor((r-1) / 10) + 1;

// 	if(i == 0){
// 		i = 10;
// 	}

// 	for(int t = 0; t < T; ++t){
// 		if(t % 2 != 0){
// 			if(!is_red[r-1]){
// 				if(use_simulation){
// 					theoretical_lambda(t,r-1, 0) =
// 5*(i+j-1);
// 				}
// 				integrated_lambda(t,r-1)
// = 2.5*(pow(i,2)-pow(i-1,2)) + 2.5*(pow(j,2) - pow(j-1,2));
// }else{ 				if(use_simulation){
// theoretical_lambda(t,r-1, 0) = i+j-1;
// 				}
// 				integrated_lambda(t,r-1) =
// 0.5*(pow(i,2)-pow(i-1,2)) + 0.5*(pow(j,2) - pow(j-1,2));
// 			}
// 		}else{
// 			if(!is_red[r-1]){
// 				if(use_simulation){
// 					theoretical_lambda(t,r-1, 0) = i+j-1;
// 				}
// 				integrated_lambda(t,r-1) =
// 0.5*(pow(i,2)-pow(i-1,2)) + 0.5*(pow(j,2) - pow(j-1,2));
// }else{ 				if(use_simulation){
// theoretical_lambda(t,r-1, 0) = 5*(i+j-1);
// 				}
// 				integrated_lambda(t,r-1)
// = 2.5*(pow(i,2)-pow(i-1,2)) + 2.5*(pow(j,2) - pow(j-1,2));
// 			}
// 		}
// 	}
// }