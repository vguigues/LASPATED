#include "missing_data.hpp"

namespace po = boost::program_options;
ResultModel1 run_model1(xt::xarray<int>& nb_observations,
                        xt::xarray<int>& nb_arrivals,
                        xt::xarray<int>& nb_missing_arrivals,
                        xt::xarray<double>& durations) {
  ulong C = nb_observations.shape(0);
  ulong D = nb_observations.shape(1);
  ulong T = nb_observations.shape(2);
  ulong R = nb_observations.shape(3);

  double without_location = 0;
  double with_location = 0;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        without_location += nb_missing_arrivals(c, d, t);
        for (int r = 0; r < R; ++r) {
          with_location += nb_arrivals(c, d, t, r);
        }
      }
    }
  }

  double prob = without_location / (without_location + with_location);
  xt::xarray<int> no_location = xt::zeros<int>({C, D, T});
  xt::xarray<int> yes_location = xt::zeros<int>({C, D, T});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        no_location(c, d, t) += nb_missing_arrivals(c, d, t);
        for (int r = 0; r < R; ++r) {
          yes_location(c, d, t) += nb_arrivals(c, d, t, r);
        }
      }
    }
  }

  xt::xarray<double> probs = xt::zeros<double>({C, D, T});
  xt::xarray<double> lambdas_missing = xt::zeros<double>({C, D, T, R});
  xt::xarray<double> lambdas_no_missing = xt::zeros<double>({C, D, T, R});

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        probs(c, d, t) = no_location(c, d, t) /
                         (no_location(c, d, t) + yes_location(c, d, t));
      }
    }
  }

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          lambdas_no_missing(c, d, t, r) =
              nb_arrivals(c, d, t, r) /
              (nb_observations(c, d, t, r) * durations(d, t));
          lambdas_missing(c, d, t, r) =
              nb_arrivals(c, d, t, r) /
              ((1 - probs(c, d, t)) * nb_observations(c, d, t, r) *
               durations(d, t));
        }
      }
    }
  }

  return ResultModel1{prob, probs, lambdas_missing, lambdas_no_missing};
}

void run_models(AppParams& app_params) {
  std::ifstream info_file(app_params.info_file, std::ios::in);
  ulong C, D, T, R, nb_regressors, nb_holidays_years;
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
  ulong nb_obs = *max_element(daily_obs.begin(), daily_obs.end());
  xt::xarray<int> sample = xt::zeros<int>({C, D, T, R, nb_obs});
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, D, T, R});
  std::string aux_str;
  std::ifstream arrivals_file(app_params.arrivals_file, std::ios::in);
  do {
    getline(arrivals_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    int t, d, r, c, j, val;
    ss >> t >> d >> r >> c >> j >> val;

    // printf("%d %d %d %d %d %d\n", t, d, r, c, j, val);
    sample(c, d, t, r, j) = val;
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
  auto neighbors = std::vector<std::vector<int>>(R, std::vector<int>());
  xt::xarray<double> distance = xt::zeros<double>({R, R});

  std::ifstream neighbors_file(app_params.neighbors_file, std::ios::in);
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

  xt::xarray<int> sample_missing_arrivals = xt::zeros<int>({C, D, T, nb_obs});
  xt::xarray<int> nb_missing_arrivals = xt::zeros<int>({C, D, T});

  std::ifstream missing_file(app_params.missing_file, std::ios::in);
  do {
    getline(missing_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    double d, r;
    int t, c, j, val;
    ss >> t >> d >> r >> c >> j >> val;

    // printf("%d %d %d %d %d %d\n", t, d, r, c, j, val);
    sample_missing_arrivals(c, static_cast<int>(d), t, j) = val;
    nb_missing_arrivals(c, static_cast<int>(d), t) += val;
  } while (true);
  missing_file.close();
  printf("Read missing_arrivals file\n");

  xt::xarray<double> durations = 0.5 * xt::ones<double>({D, T});
  auto result =
      run_model1(nb_observations, nb_arrivals, nb_missing_arrivals, durations);

  printf("Run model1: P = %f\n", result.prob);

  std::ofstream model1_lambda_file("results/lambda_model1.txt", std::ios::out);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int ind_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          model1_lambda_file << c << " " << r << " " << ind_t << " "
                             << result.lambda_no_missing(c, d, t, r) << " "
                             << result.lambda_missing(c, d, t, r) << "\n";
        }
      }
    }
  }
  model1_lambda_file.close();

  int nb_groups = 8;
  std::vector<std::vector<std::pair<int, int>>> groups(
      nb_groups, std::vector<std::pair<int, int>>());

  // Weekday morning
  for (int d = 0; d < 5; ++d) {
    for (int t = 12; t < 20; ++t) {
      groups[0].push_back(std::make_pair(d, t));
    }
  }

  // Weekday afternoon
  for (int d = 0; d < 5; ++d) {
    for (int t = 20; t < 36; ++t) {
      groups[1].push_back(std::make_pair(d, t));
    }
  }

  // Weekday evening
  for (int d = 0; d < 5; ++d) {
    for (int t = 36; t < 44; ++t) {
      groups[2].push_back(std::make_pair(d, t));
    }
  }

  // Weekday Night
  for (int d = 0; d < 4; ++d) {
    for (int t = 44; t < 48; ++t) {
      groups[3].push_back(std::make_pair(d, t));
    }
  }

  // Sunday Night
  for (int t = 44; t < 48; ++t) {
    groups[3].push_back(std::make_pair(6, t));
  }

  // Weekday Early Morning
  for (int d = 0; d < 5; ++d) {
    for (int t = 0; t < 12; ++t) {
      groups[3].push_back(std::make_pair(d, t));
    }
  }

  // Friday Night
  for (int t = 44; t < 48; ++t) {
    groups[4].push_back(std::make_pair(4, t));
  }

  // Saturday Night
  for (int t = 44; t < 48; ++t) {
    groups[4].push_back(std::make_pair(5, t));
  }

  // Saturday Early Morning
  for (int t = 0; t < 12; ++t) {
    groups[4].push_back(std::make_pair(5, t));
  }

  // Sunday Early Morning
  for (int t = 0; t < 12; ++t) {
    groups[4].push_back(std::make_pair(6, t));
  }

  // Saturday Morning
  for (int t = 12; t < 20; ++t) {
    groups[5].push_back(std::make_pair(5, t));
  }
  // Sunday Morning
  for (int t = 12; t < 20; ++t) {
    groups[5].push_back(std::make_pair(6, t));
  }

  // Saturday Afternoon
  for (int t = 20; t < 36; ++t) {
    groups[6].push_back(std::make_pair(5, t));
  }

  // Sunday Afternoon
  for (int t = 20; t < 36; ++t) {
    groups[6].push_back(std::make_pair(6, t));
  }
  // Evening
  for (int t = 36; t < 44; ++t) {
    groups[7].push_back(std::make_pair(5, t));
  }
  // Evening
  for (int t = 36; t < 44; ++t) {
    groups[7].push_back(std::make_pair(6, t));
  }

  std::vector<double> test_weights{0, 0.001, 0.005, 0.01, 0.03};
  laspated::Param param;
  param.EPS = app_params.EPS;
  param.sigma = app_params.sigma;
  param.lower_lambda = app_params.lower_lambda;
  param.max_iter = app_params.max_iter;
  param.beta_bar = app_params.beta_bar;

  std::stringstream filename;

  for (double w : test_weights) {
    std::vector<double> weights(nb_groups, w);
    xt::xarray<double> alphas = w * xt::ones<double>({R, R});
    MissingLambdaModel m1(nb_observations, nb_arrivals, nb_missing_arrivals,
                          alphas, weights, groups, neighbors, durations, param);
    xt::xarray<double> x0 = param.EPS * xt::ones<double>({C, D, T, R});

    xt::xarray<double> lambda =
        laspated::projected_gradient_armijo_feasible<MissingLambdaModel>(
            m1, param, x0);
    filename.str("");
    filename << "results/lambda_model2_w" << w << ".txt";
    std::ofstream lambda_reg_file(filename.str(), std::ios::out);
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          int ind_t = d * T + t;
          for (int r = 0; r < R; ++r) {
            lambda_reg_file << c << " " << r << " " << ind_t << " "
                            << lambda(c, d, t, r) << "\n";
          }
        }
      }
    }
    lambda_reg_file.close();
  }
  printf("Finished Model Regularized\n");
  ulong S = 30;
  xt::xarray<int> mn_samples = xt::zeros<int>({C, D, T, nb_obs, S, R});

  std::ifstream mn_samples_file("Rect10x10/mn_samples.dat", std::ios::in);
  int count = 0;
  while (true) {
    std::getline(mn_samples_file, aux_str);
    if (aux_str == "END") {
      break;
    }
    // std::cout << aux_str << "\n";
    std::istringstream ss(aux_str);
    int c, d, t, n, s, r, val;
    ss >> c >> d >> t >> n >> s >> r >> val;
    mn_samples(c, d, t, n, s, r) = val;
    ++count;
    if (count % (39836160) == 0) {
      printf("Reading mn_samples\n");
    }
  }
  mn_samples_file.close();
  printf("Read mn_samples.txt");
  MissingModel2 m2(nb_observations, nb_arrivals, nb_missing_arrivals, sample,
                   sample_missing_arrivals, mn_samples, durations, param);

  xt::xarray<double> x0 = param.EPS * xt::ones<double>({C, D, T, R});

  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int ind_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          x0(c, d, t, r) = result.lambda_missing(c, d, t, r);
        }
      }
    }
  }

  xt::xarray<double> lambda =
      laspated::projected_gradient_armijo_feasible<MissingModel2>(m2, param,
                                                                  x0);

  std::ofstream lambda_model2_file("results/lambda_model3.txt", std::ios::out);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int ind_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          lambda_model2_file << c << " " << r << " " << ind_t << " "
                             << result.lambda_missing(c, d, t, r) << " "
                             << lambda(c, d, t, r) << "\n";
        }
      }
    }
  }
  lambda_model2_file.close();
}

AppParams load_options(int argc, char* argv[]) {
  std::string config_file;
  po::variables_map vm;
  // Declare a group of options that will be
  // allowed only on command line
  po::options_description generic("Generic Options");
  generic.add_options()("help,h", "Display this help message.")(
      "file,f", po::value<std::string>()->default_value(""),
      "Path to configuration file.");

  po::options_description config("Configuration");
  config.add_options()(
      "EPS,E", po::value<double>()->default_value(1e-5),
      "Epsilon for feasibility and convergence checks. Default = 1e-5")(
      "sigma,s", po::value<double>()->default_value(0.5),
      "Sigma parameter of armijo step. Default = 0.5")(
      "max_iter,I", po::value<int>()->default_value(30),
      "Max number of iterations used in stopping criterion. Default = 30")(
      "lower_lambda, L", po::value<double>()->default_value(1e-6),
      "lower bound on decision variables for both models. Default = 1e-6")(
      "beta_bar, B", po::value<double>()->default_value(2.0),
      "Initial step size for projected gradient. Default = 2.0")(
      "info_file,i", po::value<std::string>()->default_value(""),
      "Path to file with general information about the model.")(
      "arrivals_file,a", po::value<std::string>()->default_value(""),
      "Path to file with arrivals data. Default = ''")(
      "neighbors_file,n", po::value<std::string>()->default_value(""),
      "Path to file with neighbors data. Default = ''")(
      "missing_file,n", po::value<std::string>()->default_value(""),
      "Path to file with missing arrivals data. Default = ''");

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  po::options_description config_file_options;
  config_file_options.add(config);

  po::options_description visible("Allowed Options");
  visible.add(generic).add(config);

  store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  notify(vm);

  config_file = vm["file"].as<std::string>();
  std::ifstream ifs(config_file);
  if (config_file != "" && ifs) {
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  } else if (config_file != "" && !ifs) {
    printf("Could not open config file: %s\n", config_file.c_str());
  }
  if (vm.count("help")) {
    std::cout << visible << "\n";
    exit(0);
  }

  AppParams app_params;
  app_params.EPS = vm["EPS"].as<double>();
  app_params.sigma = vm["sigma"].as<double>();
  app_params.max_iter = vm["max_iter"].as<int>();
  app_params.lower_lambda = vm["lower_lambda"].as<double>();
  app_params.beta_bar = vm["beta_bar"].as<double>();
  app_params.info_file = vm["info_file"].as<std::string>();
  app_params.arrivals_file = vm["arrivals_file"].as<std::string>();
  app_params.neighbors_file = vm["neighbors_file"].as<std::string>();
  app_params.missing_file = vm["missing_file"].as<std::string>();
  return app_params;
}

int main(int argc, char* argv[]) {
  AppParams app_params = load_options(argc, argv);

  run_models(app_params);
  return 0;
}
