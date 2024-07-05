#include "laspated.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void read_weights(std::ifstream& time_groups_file, int nb_groups,
                  std::vector<double>& weights) {
  weights = std::vector<double>(nb_groups, 0);
  for (int g = 0; g < nb_groups; ++g) {
    time_groups_file >> weights[g];
    if (weights[g] < 0) {
      printf("ERROR: parsing time groups info. Invalid weight %f\n",
             weights[g]);
      exit(1);
    }
  }
}

void read_groups(laspated::AppParameters& app_params, ulong total_T,
                 std::vector<std::vector<int>>& groups,
                 std::vector<int>& which_group) {
  std::ifstream time_groups(app_params.time_groups_file);
  if (!time_groups) {
    printf("ERROR: could not open file %s.\n",
           app_params.time_groups_file.c_str());
    exit(1);
  }

  int nb_groups;
  time_groups >> nb_groups;

  if (nb_groups < 0 || nb_groups > total_T) {
    printf("ERROR: Invalid number of groups: %d. Must be between 0 and %ld",
           nb_groups, total_T);
  }
  groups = std::vector<std::vector<int>>(nb_groups, std::vector<int>());
  which_group = std::vector<int>(total_T, -1);
  for (int t = 0; t < total_T; ++t) {
    int index_group;
    time_groups >> index_group;
    if (index_group >= nb_groups) {
      printf("ERROR: Invalid index_group %d. Must be less than %d.\n",
             index_group, nb_groups);
    }
    groups[index_group].push_back(t);
    which_group[t] = index_group;
  }
  time_groups.close();
}

void read_groups(laspated::AppParameters& app_params, ulong total_T,
                 std::vector<std::vector<int>>& groups,
                 std::vector<int>& which_group, std::vector<double>& weights) {
  std::ifstream time_groups(app_params.time_groups_file);
  if (!time_groups) {
    printf("ERROR: could not open file %s.\n",
           app_params.time_groups_file.c_str());
    exit(1);
  }

  int nb_groups;
  time_groups >> nb_groups;

  if (nb_groups < 0 || nb_groups > total_T) {
    printf("ERROR: Invalid number of groups: %d. Must be between 0 and %ld",
           nb_groups, total_T);
  }
  groups = std::vector<std::vector<int>>(nb_groups, std::vector<int>());
  which_group = std::vector<int>(total_T, -1);
  for (int t = 0; t < total_T; ++t) {
    int index_group;
    time_groups >> index_group;
    if (index_group >= nb_groups) {
      printf("ERROR: Invalid index_group %d. Must be less than %d.\n",
             index_group, nb_groups);
      exit(1);
    }
    groups[index_group].push_back(t);
    which_group[t] = index_group;
  }
  read_weights(time_groups, nb_groups, weights);
  time_groups.close();
}

void read_durations(laspated::AppParameters& app_params,
                    std::vector<double>& durations, const ulong T) {
  durations = std::vector<double>(T, 1);
  std::ifstream times_file(app_params.durations_file, std::ios::in);
  ulong file_T;
  times_file >> file_T;
  if (file_T != T) {
    printf(
        "ERROR: number of times in durations_file (%ld) is different from T "
        "(%ld).\n",
        file_T, T);
    exit(1);
  }
  for (ulong t = 0; t < T; ++t) {
    times_file >> durations[t];
  }
  times_file.close();
}

void read_alpha_matrix(laspated::AppParameters& app_params, ulong R,
                       xt::xarray<double>& alphas) {
  std::ifstream alpha_regions(app_params.alpha_regions_file);
  if (!alpha_regions) {
    printf("ERROR: Unable to open file %s.\n",
           app_params.alpha_regions_file.c_str());
  }
  for (int r = 0; r < R; ++r) {
    for (int s = 0; s < R; ++s) {
      alpha_regions >> alphas(r, s);
      if (alphas(r, s) < 0) {
        printf(
            "ERROR: invalid alpha(%d, %d) = %f. Must be greater or equal to "
            "zero.\n",
            r, s, alphas(r, s));
        exit(1);
      }
    }
  }
  alpha_regions.close();
}

void laspated_no_reg(laspated::AppParameters& app_params) {
  using namespace std;
  stringstream info_file_name;
  stringstream arrivals_file_name;
  stringstream neighbors_file_name;
  info_file_name << app_params.info_file;
  arrivals_file_name << app_params.arrivals_file;
  neighbors_file_name << app_params.neighbors_file;
  cout << "Running laspated for the Regularized Model\n";
  cout << "info_file_name = " << info_file_name.str() << "\n"
       << "arrivals_file_name = " << arrivals_file_name.str() << "\n"
       << "neighbors_file_name = " << neighbors_file_name.str() << "\n";
  ulong C, D, T, R, nb_regressors, nb_holidays_years;
  auto info_file = ifstream(info_file_name.str(), ios::in);
  info_file >> T >> D >> R >> C >> nb_regressors >> nb_holidays_years;
  std::vector<int> daily_obs(D, 0);

  // printf("%ld %ld %ld %ld %ld %ld\ndaily_obs = ", T, D, R, C, nb_regressors,
  //        nb_holidays_years);
  for (int d = 0; d < D; ++d) {
    info_file >> daily_obs[d];
    // printf("%d ", daily_obs[d]);
  }
  // printf("\n");
  info_file.close();
  printf("Read info_file\n");
  xt::xarray<int> nb_observations = xt::zeros<int>({C, D, T, R});
  ulong nb_obs = *min_element(daily_obs.begin(), daily_obs.end());
  xt::xarray<vector<int>> sample = xt::zeros<vector<int>>({
      C,
      D,
      T,
      R,
  });
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, D, T, R});

  auto arrivals_file = ifstream(arrivals_file_name.str(), ios::in);
  string aux_str;

  do {
    getline(arrivals_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    istringstream ss(aux_str);
    int t, d, r, c, j, val;
    ss >> t >> d >> r >> c >> j >> val;

    // printf("%d %d %d %d %d %d\n", t, d, r, c, j, val);
    sample(c, d, t, r).push_back(val);
    nb_observations(c, d, t, r) += 1;
    nb_arrivals(c, d, t, r) += val;
  } while (true);

  arrivals_file.close();
  printf("Read arrivals_file\n");

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          nb_observations(c, d, t, r) = daily_obs[d];
        }
      }
    }
  }

  auto type_region = std::vector<int>(R, -1);
  xt::xarray<double> regressors = xt::zeros<double>({nb_regressors, R});
  auto neighbors = std::vector<vector<int>>(R, std::vector<int>());
  xt::xarray<double> distance = xt::zeros<double>({R, R});

  auto neighbors_file = ifstream(neighbors_file_name.str(), ios::in);
  while (true) {
    int ind, terrain_type, s, land_type;
    double lat, longi, dist;
    std::getline(neighbors_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    // std::cout << aux_str << "\n";
    std::istringstream ss(aux_str);
    ss >> ind >> lat >> longi >> land_type;
    type_region[ind] = land_type;
    for (int j = 0; j < nb_regressors; ++j) {
      ss >> regressors(j, ind);
    }
    while (ss >> s >> dist) {
      distance(ind, s) = dist;
      neighbors[ind].push_back(s);
    }
  }
  neighbors_file.close();
  printf("Read neighbors_file\n");
  xt::xarray<int> nb_observations_no_cov = xt::zeros<int>({C, R, D * T});
  xt::xarray<int> nb_arrivals_no_cov = xt::zeros<int>({C, R, D * T});
  xt::xarray<int> sample_no_cov = xt::zeros<int>({D * T, R, C, nb_obs});

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          nb_observations_no_cov(c, r, index_t) = nb_observations(c, d, t, r);
          nb_arrivals_no_cov(c, r, index_t) = nb_arrivals(c, d, t, r);
          for (int k = 0; k < sample(c, d, t, r).size() && k < nb_obs; ++k) {
            sample_no_cov(index_t, r, c, k) = sample(c, d, t, r)[k];
          }
        }
      }
    }
  }

  ifstream groups_file(app_params.time_groups_file);
  if (!groups_file) {
    printf("ERROR: Could not open time_groups_file %s\n",
           app_params.time_groups_file.c_str());
    exit(1);
  }
  vector<int> which_group(D * T, -1);
  vector<vector<int>> groups;
  vector<double> weights;
  std::vector<double> durations_no_cov;
  read_durations(app_params, durations_no_cov, D * T);
  printf("Read durations\n");
  if (app_params.model_type == "no_reg") {
    read_groups(app_params, D * T, groups, which_group, weights);
  }
  printf("Read groups_file\n");
  xt::xarray<double> alphas = xt::zeros<double>({R, R});
  if (app_params.model_type == "no_reg" && app_params.method == "calibration") {
    read_alpha_matrix(app_params, R, alphas);
  }
  printf("Read alpha_file\n");

  using laspated::Param;
  using laspated::RegularizedModel;
  Param param(app_params);
  xt::xarray<double> lambda = xt::zeros<double>(nb_observations_no_cov.shape());
  xt::xarray<double> lambda0 = param.EPS * xt::ones<double>(lambda.shape());
  if (app_params.method == "calibration") {
    RegularizedModel m(nb_observations_no_cov, nb_arrivals_no_cov,
                       durations_no_cov, groups, weights, alphas, distance,
                       type_region, neighbors, param);
    printf("Running Projected Gradient %s\n", app_params.algorithm.c_str());
    if (app_params.algorithm == "feasible") {
      lambda = laspated::projected_gradient_armijo_feasible<RegularizedModel>(
          m, param, lambda0);
    } else if (app_params.algorithm == "boundary") {
      lambda = laspated::projected_gradient_armijo_boundary<RegularizedModel>(
          m, param, lambda0);
    } else {
      printf("ERROR: Unknown algorithm %s\n", app_params.algorithm.c_str());
      exit(1);
    }
  } else if (app_params.method == "cross_validation") {
    ifstream weights_file(app_params.cv_weights_file);
    if (!weights_file) {
      printf("ERROR: Could not open cv_weights_file %s\n",
             app_params.cv_weights_file.c_str());
      exit(1);
    }
    std::vector<double> test_weights;
    while (!weights_file.eof()) {
      double weight;
      weights_file >> weight;
      if (weight < 0.0) {
        printf(
            "ERROR: Parsing cv_weights_file %s: invalid weight %f. Weights "
            "must be positive.\n",
            app_params.cv_weights_file.c_str(), weight);
        exit(1);
      }
      test_weights.push_back(weight);
    }
    weights_file.close();
    using laspated::cross_validation;
    std::vector<double> weights(D * T, 0);
    RegularizedModel m1(nb_observations_no_cov, nb_arrivals_no_cov,
                        durations_no_cov, groups, weights, alphas, distance,
                        type_region, neighbors, param);
    printf("Running Cross Validation with Projected Gradient %s\n",
           app_params.algorithm.c_str());
    auto cv_result = cross_validation(param, m1, sample_no_cov, test_weights);
    lambda = cv_result.lambda;
    printf("Finished cross validation in %f seconds. Best weight = %f\n",
           cv_result.wall_time, cv_result.weight);
  } else {
    printf("ERROR: Unknown method: %s\n", app_params.method.c_str());
    exit(1);
  }
  ofstream output(app_params.output_file);
  if (!output) {
    printf("ERROR: Could not open outputfile %s.\n",
           app_params.output_file.c_str());
    exit(1);
  }
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        output << c << " " << r << " " << t << " " << lambda(c, r, t) << "\n";
      }
    }
  }
  printf("Intensities saved at %s\n", app_params.output_file.c_str());
}

#if USE_GUROBI == 1
void laspated_reg(laspated::AppParameters& app_params) {
  using namespace std;
  stringstream info_file_name;
  stringstream arrivals_file_name;
  stringstream neighbors_file_name;
  info_file_name << app_params.info_file;
  arrivals_file_name << app_params.arrivals_file;
  neighbors_file_name << app_params.neighbors_file;

  cout << info_file_name.str() << "\n"
       << arrivals_file_name.str() << "\n"
       << neighbors_file_name.str() << "\n";
  ulong C, D, T, R, nb_regressors, nb_holidays_years;
  auto info_file = ifstream(info_file_name.str(), ios::in);
  info_file >> T >> D >> R >> C >> nb_regressors >> nb_holidays_years;
  std::vector<int> daily_obs(D, 0);

  printf("%ld %ld %ld %ld %ld %ld\ndaily_obs = ", T, D, R, C, nb_regressors,
         nb_holidays_years);
  for (int d = 0; d < D; ++d) {
    info_file >> daily_obs[d];
    printf("%d ", daily_obs[d]);
  }
  printf("\n");
  info_file.close();

  xt::xarray<int> nb_observations = xt::zeros<int>({C, D, T, R});
  ulong nb_obs = *min_element(daily_obs.begin(), daily_obs.end());
  xt::xarray<vector<int>> sample = xt::zeros<vector<int>>({
      C,
      D,
      T,
      R,
  });
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, D, T, R});

  auto arrivals_file = ifstream(arrivals_file_name.str(), ios::in);
  string aux_str;

  do {
    getline(arrivals_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    istringstream ss(aux_str);
    int t, d, r, c, j, val;
    ss >> t >> d >> r >> c >> j >> val;
    // fmt::print("calls {} {} {} {} {} {} {}\n",t,d,r,c,j,val,h);
    // cin.get();
    sample(c, d, t, r).push_back(val);
    nb_observations(c, d, t, r) += 1;
    nb_arrivals(c, d, t, r) += val;
  } while (true);

  arrivals_file.close();
  auto durations = vector<double>(T, 1);
  read_durations(app_params, durations, T);
  xt::xarray<double> empirical_rates = xt::zeros<double>({C, R, D * T});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          nb_observations(c, d, t, r) = daily_obs[d];
          empirical_rates(c, r, d * T + t) =
              static_cast<double>(nb_arrivals(c, d, t, r)) /
              (nb_observations(c, d, t, r));
          empirical_rates(c, r, d * T + t) /= durations[t];
          // printf("c%d r%d t%ld, emp = %f, obs = %d arr = %d\n", c, r, d * T
          // + t,
          //        empirical_rates(c, r, d * T + t), nb_observations(c, d, t,
          //        r), nb_arrivals(c, d, t, r));
        }
      }
    }
  }
  // cin.get();
  auto type_region = std::vector<int>(R, -1);
  xt::xarray<double> regressors = xt::zeros<double>({nb_regressors, R});
  auto neighbors = std::vector<vector<int>>(R, std::vector<int>());
  xt::xarray<double> distance = xt::zeros<double>({R, R});

  auto neighbors_file = ifstream(neighbors_file_name.str(), ios::in);
  while (true) {
    int ind, terrain_type, s, land_type;
    double lat, longi, dist;
    std::getline(neighbors_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    // std::cout << aux_str << "\n";
    std::istringstream ss(aux_str);
    ss >> ind >> lat >> longi >> land_type;
    type_region[ind] = land_type;
    for (int j = 0; j < nb_regressors; ++j) {
      ss >> regressors(j, ind);
    }
    while (ss >> s >> dist) {
      distance(ind, s) = dist;
      neighbors[ind].push_back(s);
    }
  }
  neighbors_file.close();

  laspated::Param param(app_params);
  param.upper_lambda = 1e3;
  param.EPS = 1e-5;
  param.max_iter = 30;

  laspated::CovariatesModel m(nb_observations, nb_arrivals, regressors, param);

  xt::xarray<double> beta0 =
      2 * pow(10.0, -3) * xt::ones<double>({C, D, T, nb_regressors});
  xt::xarray<double> beta = xt::zeros<double>(beta0.shape());

  printf("Running projected gradient %s\n", app_params.algorithm.c_str());
  if (app_params.algorithm == "boundary") {
    beta =
        laspated::projected_gradient_armijo_boundary<laspated::CovariatesModel>(
            m, param, beta0);
  } else if (app_params.algorithm == "feasible") {
    beta =
        laspated::projected_gradient_armijo_feasible<laspated::CovariatesModel>(
            m, param, beta0);
  } else {
    printf("Error: unknown algorithm %s.\n", app_params.algorithm.c_str());
    exit(1);
  }

  ofstream output(app_params.output_file);
  if (!output) {
    printf("ERROR: Could not open outputfile %s.\n",
           app_params.output_file.c_str());
    exit(1);
  }

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int j = 0; j < nb_regressors; ++j) {
          output << c << " " << d << " " << t << " " << j << " "
                 << beta(c, d, t, j) << "\n";
        }
      }
    }
  }
  printf("Intensities saved at %s\n", app_params.output_file.c_str());
}
#endif

void run_model(laspated::AppParameters& app_params) {
  if (app_params.model_type == "no_reg") {
    // Validation
    if (app_params.method == "cross_validation" &&
        app_params.algorithm == "boundary") {
      printf(
          "Warning: Cross validation only supports projected gradient along "
          "the boundary.\n");
    }
    laspated_no_reg(app_params);
  } else {
#if USE_GUROBI == 1
    if (app_params.method == "cross_validation") {
      printf("Error: Cross validation not implemented for model_type reg.\n");
      exit(1);
    }
    laspated_reg(app_params);
#else
    printf(
        "ERROR: This executable was not compiled with Gurobi support. Please "
        "see INSTALL.md\n");
    exit(1);
#endif
  }
}

int main(int argc, char* argv[]) {
  po::variables_map vm;
  laspated::AppParameters app_params = laspated::load_options(argc, argv, vm);
  run_model(app_params);
  return 0;
}
