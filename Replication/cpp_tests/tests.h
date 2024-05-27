#include <fstream>
#include <sstream>

#include "laspated.h"

typedef struct {
  double mean_emp;
  double min_error;
  int min_w;
  double mean_cv;
  double w_cv;
  std::vector<double> error_by_weight;
  ulong C;
  ulong R;
  ulong T;
  xt::xarray<double> theoretical_lambda;
  xt::xarray<double> empirical_lambda;
  xt::xarray<double> estimated_lambda;
} Result1;

Result1 test1(int nb_weeks, int nb_groups, int neighbor_factor,
              bool constant_lambdas, std::vector<double> &test_weights,
              bool do_cross_validation) {
  using namespace std;

  int x_max = 10;
  int y_max = 10;
  int n_x = 10;
  int n_y = 10;
  ulong C = 1;
  ulong R = n_x * n_y;
  ulong T = 4 * 7;

  vector<vector<int>> groups = vector<vector<int>>(nb_groups, vector<int>());
  vector<int> which_group(T, -1);
  for (int t = 0; t < T; ++t) {
    groups[t % nb_groups].push_back(t);
    which_group[t] = t % nb_groups;
  }

  // std::random_device rd;
  // std::mt19937 gen(rd());

  std::default_random_engine gen;

  vector<bool> is_red(n_x * n_y, false);
  vector<int> type_region = vector<int>(R, -1);
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
      cout << "ERROR: Impossible region " << r << "(" << x << "," << y << ")\n";
      exit(1);
    }
  }

  int nb_observations_total = nb_weeks;
  xt::xarray<int> sample =
      xt::zeros<int>({T, R, C, static_cast<ulong>(nb_observations_total)});

  vector<double> durations(T, 1);
  xt::xarray<int> nb_observations =
      nb_observations_total * xt::ones<int>({C, R, T});
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, R, T});
  xt::xarray<double> empirical_lambda = xt::zeros<double>({C, R, T});
  xt::xarray<double> theoretical_lambda = xt::zeros<double>({C, R, T});
  for (int t = 0; t < T; ++t) {
    for (int r = 0; r < R; ++r) {
      int x = r % 10;
      int y = r / 10;
      double cx = x + 0.5;
      double cy = y + 0.5;
      for (int c = 0; c < C; ++c) {
        if ((!is_red[r] && t % 2 == 0) || (is_red[r] && t % 2 != 0)) {
          theoretical_lambda(c, r, t) = (constant_lambdas) ? 0.1 : (cx + cy);
        } else if ((!is_red[r] && t % 2 != 0) || (is_red[r] && t % 2 == 0)) {
          theoretical_lambda(c, r, t) =
              (constant_lambdas) ? 0.5 : 5 * (cx + cy);
        } else {
          std::cout << "ERROR: Unpredictable lambda case at t = ";
          std::cout << t << " and r = " << r << "(" << cx << "," << cy << ")\n";
          exit(1);
        }
        // printf("rate r%d t%d: %f\n", r, t, theoretical_lambda(c, r, t));
        for (int n = 0; n < nb_observations_total; ++n) {
          poisson_distribution<int> pd(theoretical_lambda(c, r, t) *
                                       durations[t]);
          int this_nb_arrival = pd(gen);
          sample(t, r, c, n) = this_nb_arrival;
          nb_arrivals(c, r, t) += this_nb_arrival;
        }
        empirical_lambda(c, r, t) =
            nb_arrivals(c, r, t) /
            static_cast<double>(nb_observations(c, r, t) * durations[t]);
      }
    }
    // cin.get();
  }

  vector<vector<int>> neighbors = vector<vector<int>>(R, vector<int>());
  xt::xarray<double> distance = GRB_INFINITY * xt::ones<double>({R, R});
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

    //   if ((Yi + 1 <= n_y) && (Xi - 1 >= 1)) {
    //     neighbors[r].push_back(Yi * n_x + Xi - 1 - 1);
    //   }

    //   if ((Yi - 1 >= 1) && (Xi - 1) >= 1) {
    //     neighbors[r].push_back((Yi - 2) * n_x + Xi - 1 - 1);
    //   }

    //   if ((Yi + 1 <= n_y) && (Xi + 1) <= n_x) {
    //     neighbors[r].push_back(Yi * n_x + Xi + 1 - 1);
    //   }

    //   if ((Yi - 1 >= 1) && (Xi + 1) <= n_x) {
    //     neighbors[r].push_back((Yi - 2) * n_x + Xi + 1 - 1);
    //   }
  }

  for (int r = 0; r < R; ++r) {
    for (int s : neighbors[r]) {
      distance(r, s) = 1;
    }
  }

  using laspated::Param;
  using laspated::projected_gradient_armijo_feasible;
  using laspated::RegularizedModel;

  Param param;
  double min_err = 1e100;
  int min_w = -1;
  double mean_emp = 0;
  vector<double> err_by_weight;
  xt::xarray<double> min_lambda;

  param.EPS = 0.001;
  param.sigma = 0.5;
  param.beta_bar = 2;
  param.max_iter = 100;
  param.upper_lambda = 1e3;
  param.lower_lambda = param.EPS;
  xt::xarray<double> x0 = param.EPS * xt::ones<double>({C, R, T});

  for (size_t i = 0; i < test_weights.size(); ++i) {
    double w = test_weights[i];
    type_region = vector<int>(R, 0);
    if (constant_lambdas) {
      param.max_iter = 30;
      x0 = param.EPS * xt::ones<double>({C, R, T});
    }

    xt::xarray<double> alphas = w * neighbor_factor * xt::ones<double>({R, R});

    vector<double> weights = vector<double>(groups.size(), w);
    RegularizedModel model(nb_observations, nb_arrivals, durations, groups,
                           weights, alphas, distance, type_region, neighbors,
                           param);

    xt::xarray<double> x =
        projected_gradient_armijo_feasible<RegularizedModel>(model, param, x0);
    double err = model.average_rate_difference(theoretical_lambda, x);

    cout << "\tw = " << w << ", err = " << err << "\n";
    // cin.get();
    err_by_weight.push_back(err);
    if (i == 0) {
      mean_emp = err;
    }
    if (err < min_err) {
      min_err = err;
      min_w = i;
      min_lambda = x;
    }
    x0 = x;
  }
  double mean_cv = -1;
  double w_cv = -1;
  if (do_cross_validation) {
    // alphas and weights will be set by cross_validation
    xt::xarray<double> alphas = xt::zeros<double>({R, R});
    vector<double> weights = vector<double>(groups.size(), 0);
    type_region = vector<int>(R, 0);

    param.EPS = 0.001;
    param.cv_proportion = 0.2;
    param.max_iter = 30;
    param.sigma = 0.5;
    param.beta_bar = 1;
    param.upper_lambda = 1e3;
    param.lower_lambda = param.EPS;
    std::vector<double> test_alphas = test_weights;

    RegularizedModel model(nb_observations, nb_arrivals, durations, groups,
                           weights, alphas, distance, type_region, neighbors,
                           param);
    auto result = cross_validation(param, model, sample, test_weights);
    mean_cv = model.average_rate_difference(theoretical_lambda, result.lambda);
    w_cv = result.weight;
  }
  if (constant_lambdas && (nb_weeks == 1 || nb_weeks == 10) &&
      neighbor_factor == 1 && nb_groups == 2) {
    param.max_iter = 100;
    xt::xarray<double> alphas = test_weights[min_w] * xt::ones<double>({R, R});
    vector<double> weights = vector<double>(groups.size(), test_weights[min_w]);
    RegularizedModel model(nb_observations, nb_arrivals, durations, groups,
                           weights, alphas, distance, type_region, neighbors,
                           param);
    min_lambda = projected_gradient_armijo_feasible<RegularizedModel>(
        model, param, min_lambda);
  }

  return Result1{
      mean_emp,         min_err,   min_w, mean_cv, w_cv,
      err_by_weight,    C,         R,     T,       theoretical_lambda,
      empirical_lambda, min_lambda};
}

typedef struct {
  double err_emp;
  double err_cov;
  double min_err1;
  double min_err2;
  double min_err3;
} Result2;

double get_population(int r, std::vector<bool> &is_blue,
                      std::vector<std::vector<double>> &u) {
  double pi_app = 3.14159;
  // int i = (r+1) % 10;
  // int j = (r+1) / 10;
  int i = (r + 1) % 10;
  int j = 1 + floor((r + 1) / 10.0);
  if (i == 0) {
    i = 10;
    j = (r + 1) / 10;
  }
  double val = 0.0;
  for (int k = 1; k <= 10; ++k) {
    double ceil_val_i = ceil((k / 5.0) * (i - 1));  // ceil_val_i == int3
    double floor_val_i = floor((k / 5.0) * i);      // floor_val_i == int4
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
  }
  return val;
}

Result2 test2(int nb_years, bool use_holidays) {
  using namespace std;
  int x_max = 10;
  int y_max = 10;
  int n_x = 10;
  int n_y = 10;
  int nb_holidays_years = 8;
  ulong R = n_x * n_y;
  ulong C = 1;
  ulong T = 4;
  ulong D = 7 + nb_holidays_years * (use_holidays);

  int nb_weeks = 52;  // fixed number of weeks.
  ulong nb_obs = nb_weeks * 7;
  vector<double> durations(T, 1);
  vector<pair<bool, int>> is_holidays(nb_years * nb_weeks * 7,
                                      make_pair(false, -1));
  vector<int> days_h;
  if (use_holidays) {
    for (int i = 0; i < nb_holidays_years; ++i) {
      days_h.push_back(i);
    }

    // first nb_holidays_years days are holidays.
    for (int year = 0; year < nb_years; ++year) {
      for (int k = 0; k < days_h.size(); ++k) {
        is_holidays[year * nb_weeks * 7 + days_h[k]] = make_pair(true, k);
      }
    }
  }

  int nb_land_types = 2;
  ulong nb_regressors = 1 + nb_land_types;

  xt::xarray<double> theoretical_beta =
      xt::zeros<double>({C, D, T, nb_regressors});
  xt::xarray<double> regressors = xt::zeros<double>({nb_regressors, R});

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
  if (use_holidays) {
    for (int d = 7; d < D; ++d) {
      theoretical_beta(0, d, 1, 0) = 0.1;
      theoretical_beta(0, d, 3, 0) = 0.1;

      theoretical_beta(0, d, 0, 1) = 12;
      theoretical_beta(0, d, 1, 1) = 36;
      theoretical_beta(0, d, 2, 1) = 12;
      theoretical_beta(0, d, 3, 1) = 36;

      theoretical_beta(0, d, 0, 2) = 6;
      theoretical_beta(0, d, 1, 2) = 12;
      theoretical_beta(0, d, 2, 2) = 6;
      theoretical_beta(0, d, 3, 2) = 12;
    }
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> rnd_blue(0, 1);
  std::uniform_real_distribution<double> rnd_red(0, 1);

  vector<vector<double>> u(2, vector<double>(20, 0));
  for (int k = 0; k < 20; ++k) {
    u[0][k] = rnd_blue(gen);
    // u[0][k] = 1;
  }

  for (int k = 0; k < 20; ++k) {
    u[1][k] = rnd_red(gen);
    // u[1][k] = 1;
  }

  // u = {{0.3734, 0.6277, 0.2767, 0.9151, 0.0255, 0.8281, 0.2148,
  //       0.1626, 0.9430, 0.8117, 0.9916, 0.8755, 0.6834, 0.6933,
  //       0.3374, 0.7807, 0.2081, 0.0578, 0.2721, 0.3194},
  //      {0.2906, 0.5534, 0.9220, 0.3887, 0.3733, 0.3287, 0.0366,
  //       0.4025, 0.9484, 0.3023, 0.1957, 0.6760, 0.4114, 0.9227,
  //       0.5340, 0.4614, 0.4626, 0.2645, 0.0999, 0.4704}};

  vector<int> type_region(R, -1);
  vector<bool> is_blue(100, false);

  int count = 0;
  for (int i = 0; i < n_y / 2; ++i) {
    for (int j = 0; j < n_x / 2; ++j) {
      type_region[count] = 0;
      is_blue[count] = false;
      // regressors(0, count) = 50 + 50*rnd(gen);
      // regressors(0, count) = get_population(count, is_blue, u);
      regressors(0, count) = 1;
      regressors(1, count) = 0.5;
      regressors(2, count) = 0.25;
      ++count;
    }

    for (int j = 0; j < n_x / 2; ++j) {
      type_region[count] = 1;
      is_blue[count] = true;
      // regressors(0, count) = 50*rnd(gen);
      // regressors(0, count) = get_population(count, is_blue, u);
      regressors(0, count) = 1;
      regressors(1, count) = 0.25;
      regressors(2, count) = 0.5;
      ++count;
    }
  }

  for (int i = n_y / 2; i < n_y; ++i) {
    for (int j = 0; j < n_x / 2; ++j) {
      type_region[count] = 1;
      is_blue[count] = true;
      // regressors(0, count) = 50*rnd(gen);
      // regressors(0, count) = get_population(count, is_blue, u);
      regressors(0, count) = 1;
      regressors(1, count) = 0.25;
      regressors(2, count) = 0.5;
      ++count;
    }

    for (int j = 0; j < n_x / 2; ++j) {
      type_region[count] = 0;
      is_blue[count] = false;
      // regressors(0, count) = 50 + 50 * rnd(gen);
      // regressors(0, count) = get_population(count, is_blue, u);
      regressors(0, count) = 1;
      regressors(1, count) = 0.5;
      regressors(2, count) = 0.25;
      ++count;
    }
  }

  xt::xarray<vector<int>> sample = xt::zeros<vector<int>>({C, D, T, R});
  xt::xarray<int> nb_observations = xt::zeros<int>({C, D, T, R});
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, D, T, R});
  for (int index = 0; index < nb_years * nb_weeks * 7; ++index) {
    int day = index % 7;
    if (is_holidays[index].first) {
      day = 7 + is_holidays[index].second;
    }

    for (int c = 0; c < C; ++c) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += theoretical_beta(c, day, t, j) * regressors(j, r);
          }
          // printf("index = %d, day = %d, c = %d, r = %d, t = %d: rate = %f\n",
          //        index + 1, day + 1, c + 1, r + 1, t + 1, rate);
          poisson_distribution<int> pd(rate / durations[t]);
          int this_nb_arrival = pd(gen);
          sample(c, day, t, r).push_back(this_nb_arrival);
          ++nb_observations(c, day, t, r);
          nb_arrivals(c, day, t, r) += this_nb_arrival;
        }
      }
    }
  }
  // printf("nb_obs = %d\n", nb_observations(0, 7, 0, 0));
  // cout << "END HOLIDAYS\n";
  // // cin.get();
  xt::xarray<double> empirical_rate = xt::zeros<double>({C, D, T, R});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          empirical_rate(c, d, t, r) =
              static_cast<double>(nb_arrivals(c, d, t, r)) /
              (nb_observations(c, d, t, r) * durations[t]);

          double rate = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          // cout << "c" << c << ", d" << d << ", t" << t << ", r" << r << ": "
          //      << "emp = " << empirical_rate(c, d, t, r)
          //      << ", theo = " << rate / durations[t] << "\n";
        }
      }
    }
  }

  vector<double> difference_l2;
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        for (int r = 0; r < R; ++r) {
          double rate2 = empirical_rate(c, d, t, r);
          double rate1 = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate1 += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          difference_l2.push_back(abs(rate1 - rate2) / rate1);
        }
      }
    }
  }
  double err_emp = accumulate(difference_l2.begin(), difference_l2.end(), 0.0) /
                   difference_l2.size();

  // printf("Error emp = %f\n", err_emp);
  // cout << "END BETA/EMP\n";
  // cin.get();

  // TODO: Check table 2
  using laspated::CovariatesModel;
  using laspated::Param;
  using laspated::projected_gradient_armijo_feasible;
  using laspated::RegularizedModel;
  // Cov: Run with covariates
  Param param;
  param.EPS = 1e-6;
  param.lower_lambda = 1e-5;
  param.upper_lambda = 1e3;
  param.max_iter = (nb_years == 1) ? 100 : nb_years * 20;
  CovariatesModel model(nb_observations, nb_arrivals, regressors, param);
  xt::xarray<double> x0 =
      2 * pow(10.0, -3) * xt::ones<double>({C, D, T, nb_regressors});
  auto x =
      projected_gradient_armijo_feasible<CovariatesModel>(model, param, x0);

  double err_cov = model.average_rate_difference(theoretical_beta, x);
  // printf("Err cov = %.6f\n", err_cov);
  // cin.get();

  // Setup for no covariates experiments
  xt::xarray<double> theoretical_lambda = xt::zeros<double>({C, R, D * T});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          double rate = 0.0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += theoretical_beta(c, d, t, j) * regressors(j, r);
          }
          theoretical_lambda(c, r, index_t) = rate / durations[t];
        }
      }
    }
  }

  vector<double> test_weights{0,   0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
                              1.6, 1.8, 2,   2.2, 2.4, 2.6, 2.8, 3.0};
  for (size_t i = 0; i < test_weights.size(); ++i) {
    test_weights[i] /= nb_weeks;
  }

  int nb_groups = (use_holidays) ? 4 : 2;
  vector<vector<int>> groups(nb_groups, vector<int>());
  vector<int> which_group(D * T, -1);
  if (use_holidays) {
    for (int d = 0; d < 7; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        groups[index_t % 2].push_back(index_t);
        which_group[index_t] = index_t % 2;
      }
    }
    for (int d = 7; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        groups[(index_t % 2) + 2].push_back(index_t);
        which_group[index_t] = (index_t % 2) + 2;
      }
    }
  } else {
    for (int d = 0; d < 7; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        groups[index_t % 2].push_back(index_t);
        which_group[index_t] = index_t % 2;
      }
    }
  }

  xt::xarray<int> nb_obs_no_cov1 = xt::zeros<int>({C, R, D * T});
  xt::xarray<int> nb_arrivals_no_cov1 = xt::zeros<int>({C, R, D * T});
  vector<double> durations1(D * T, -1);
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        durations1[index_t] = durations[t];
        for (int r = 0; r < R; ++r) {
          nb_obs_no_cov1(c, r, index_t) = nb_observations(c, d, t, r);
          nb_arrivals_no_cov1(c, r, index_t) = nb_arrivals(c, d, t, r);
        }
      }
    }
  }

  vector<vector<int>> neighbors(R, vector<int>());
  xt::xarray<double> distance = GRB_INFINITY * xt::ones<double>({R, R});
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
  for (int r = 0; r < R; ++r) {
    for (int s = 0; s < R; ++s) {
      distance(r, s) = 1;
    }
  }
  // cin.get();
  double min_err1 = 1e100;
  int min_w1 = 0;
  xt::xarray<double> x0_no_cov = 2 * 0.001 * xt::ones<double>({C, R, D * T});
  for (size_t i = 0; i < test_weights.size(); ++i) {
    vector<double> weights(nb_groups, test_weights[i]);
    xt::xarray<double> alphas = xt::zeros<double>({R, R});
    RegularizedModel m1(nb_obs_no_cov1, nb_arrivals_no_cov1, durations1, groups,
                        weights, alphas, distance, type_region, neighbors,
                        param);
    xt::xarray<double> x = projected_gradient_armijo_feasible<RegularizedModel>(
        m1, param, x0_no_cov);

    double err = m1.average_rate_difference(theoretical_lambda, x);
    if (err < min_err1) {
      min_err1 = err;
      min_w1 = i;
      x0_no_cov = x;
    }
  }
  // NoCov2: regularized model with usual space regularization
  double min_err2 = 1e100;
  int min_w2 = 0;
  x0_no_cov = 2 * 0.001 * xt::ones<double>({C, R, D * T});
  for (size_t i = 0; i < test_weights.size(); ++i) {
    vector<double> weights(nb_groups, test_weights[i]);
    xt::xarray<double> alphas = test_weights[i] * xt::ones<double>({R, R});
    RegularizedModel m1(nb_obs_no_cov1, nb_arrivals_no_cov1, durations1, groups,
                        weights, alphas, distance, type_region, neighbors,
                        param);
    xt::xarray<double> x = projected_gradient_armijo_feasible<RegularizedModel>(
        m1, param, x0_no_cov);

    double err = m1.average_rate_difference(theoretical_lambda, x);
    if (err < min_err2) {
      min_err2 = err;
      min_w2 = i;
      x0_no_cov = x;
    }
  }

  // NoCov3: regularized model with space regularization without colors
  for (int r = 0; r < R; ++r) {
    type_region[r] = 0;
  }
  double min_err3 = 1e100;
  int min_w3 = 0;
  x0_no_cov = 2 * 0.001 * xt::ones<double>({C, R, D * T});
  for (size_t i = 0; i < test_weights.size(); ++i) {
    vector<double> weights(nb_groups, test_weights[i]);
    xt::xarray<double> alphas = test_weights[i] * xt::ones<double>({R, R});
    RegularizedModel m1(nb_obs_no_cov1, nb_arrivals_no_cov1, durations1, groups,
                        weights, alphas, distance, type_region, neighbors,
                        param);
    xt::xarray<double> x = projected_gradient_armijo_feasible<RegularizedModel>(
        m1, param, x0_no_cov);

    double err = m1.average_rate_difference(theoretical_lambda, x);
    if (err < min_err3) {
      min_err3 = err;
      min_w3 = i;
      x0_no_cov = x;
    }
  }

  return Result2{err_emp, err_cov, min_err1, min_err2, min_err3};
}

typedef struct {
  ulong C, R, T;
  xt::xarray<double> empirical_rates;
  xt::xarray<double> regularized_rates;
  xt::xarray<double> covariates_rates;
  std::vector<double> durations;
} Result3;

Result3 test3(std::string &base_path) {
  using namespace std;
  stringstream info_file_name;
  stringstream arrivals_file_name;
  stringstream neighbors_file_name;
  info_file_name << base_path << "/info.dat";
  arrivals_file_name << base_path << "/arrivals.dat";
  neighbors_file_name << base_path << "/neighbors.dat";

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
  auto durations = vector<double>(T, 0.5);
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
          // printf("c%d r%d t%ld, emp = %f, obs = %d arr = %d\n", c, r, d * T +
          // t,
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

  // Setup for no covariates experiments
  xt::xarray<int> nb_observations_no_cov = xt::zeros<int>({C, R, D * T});
  xt::xarray<int> nb_arrivals_no_cov = xt::zeros<int>({C, R, D * T});
  xt::xarray<int> sample_no_cov = xt::zeros<int>({D * T, R, C, nb_obs});
  vector<double> durations_no_cov(D * T, 0.5);
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

  vector<int> which_group(D * T, 0);
  vector<vector<int>> groups(D * T, vector<int>());
  for (int t = 0; t < D * T; ++t) {
    groups[t].push_back(t);
    which_group[t] = t;
    // fmt::print("which_group {} = {}\n", t+1, which_group[t]+1);
  }

  vector<double> test_weights{0, 0.0001, 0.0003, 0.0004, 0.0008, 0.001};
  for (size_t i = 0; i < test_weights.size(); ++i) {
    test_weights[i] /= nb_obs;
  }

  vector<double> test_alphas = test_weights;
  using laspated::cross_validation;
  // using laspated::cross_validation2;
  using laspated::Param;
  using laspated::RegularizedModel;
  Param param;
  param.upper_lambda = 1e3;
  param.EPS = 1e-5;
  param.max_iter = 30;
  int t2 = regressors.shape(1);
  for (int i = 0; i < t2; ++i) {
    double sum = 0;
    for (int j = 0; j < regressors.shape(0); ++j) {
      sum += regressors(j, i);
    }
    if (sum < param.EPS) {
      regressors(nb_regressors - 1, i) = param.EPS;
    }
  }

  // Real weights and alphas will be set by cross_validation
  vector<double> weights(groups.size(), 0);
  xt::xarray<double> alphas = xt::zeros<double>({R, R});

  RegularizedModel m1(nb_observations_no_cov, nb_arrivals_no_cov,
                      durations_no_cov, groups, weights, alphas, distance,
                      type_region, neighbors, param);

  cout << "Running Cross Validation\n";
  auto cv_result = cross_validation(param, m1, sample_no_cov, test_weights);
  xt::xarray<double> regularized_rates = cv_result.lambda;
  // cout << "CV best Weight = " << cv_result.weight << "\n";
  // cout << "OBJ CV = " << m1.f(regularized_rates) << "\n";
  // Covariates test
  using laspated::CovariatesModel;

  CovariatesModel m2(nb_observations, nb_arrivals, regressors, param);
  xt::xarray<double> x0 =
      2 * pow(10.0, -3) * xt::ones<double>({C, D, T, nb_regressors});
  cout << "Running Model with Covariates\n";
  xt::xarray<double> x =
      laspated::projected_gradient_armijo_feasible<CovariatesModel>(m2, param,
                                                                    x0);

  // cout << "OBJ REG = " << m2.f(x) << "\n";
  xt::xarray<double> covariates_rates = xt::zeros<double>({C, R, D * T});
  for (int c = 0; c < C; ++c) {
    for (int d = 0; d < D; ++d) {
      for (int t = 0; t < T; ++t) {
        int index_t = d * T + t;
        for (int r = 0; r < R; ++r) {
          double rate = 0;
          for (int j = 0; j < nb_regressors; ++j) {
            rate += x(c, d, t, j) * regressors(j, r);
          }
          covariates_rates(c, r, index_t) = rate / durations[t];
        }
      }
    }
  }

  return Result3{C,
                 R,
                 D * T,
                 empirical_rates,
                 regularized_rates,
                 covariates_rates,
                 durations_no_cov};
}

Result1 test1_deterministic(std::string &filename, int nb_weeks, int nb_groups,
                            int neighbor_factor, bool constant_lambdas,
                            std::vector<double> &test_weights,
                            bool do_cross_validation) {
  using namespace std;

  int x_max = 10;
  int y_max = 10;
  int n_x = 10;
  int n_y = 10;
  ulong C = 1;
  ulong R = n_x * n_y;
  ulong T = 4 * 7;

  vector<vector<int>> groups = vector<vector<int>>(nb_groups, vector<int>());
  vector<int> which_group(T, -1);
  for (int t = 0; t < T; ++t) {
    groups[t % nb_groups].push_back(t);
    which_group[t] = t % nb_groups;
  }

  vector<bool> is_red(n_x * n_y, false);
  vector<int> type_region = vector<int>(R, -1);
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
      cout << "ERROR: Impossible region " << r << "(" << x << "," << y << ")\n";
      exit(1);
    }
  }

  int nb_observations_total = nb_weeks;
  xt::xarray<int> sample =
      xt::zeros<int>({T, R, C, static_cast<ulong>(nb_observations_total)});

  vector<double> durations(T, 1);
  xt::xarray<int> nb_observations =
      nb_observations_total * xt::ones<int>({C, R, T});
  xt::xarray<int> nb_arrivals = xt::zeros<int>({C, R, T});
  xt::xarray<double> empirical_lambda = xt::zeros<double>({C, R, T});
  xt::xarray<double> theoretical_lambda = xt::zeros<double>({C, R, T});

  ifstream sample1_1(filename, ios::in);
  // printf("Total indexes = %ld\n", C * R * T * nb_weeks);
  for (ulong i = 0; i < C * R * T * nb_weeks; ++i) {
    int t, c, r, j, val;
    sample1_1 >> t >> c >> r >> j >> val;
    sample(t - 1, r - 1, c - 1, j - 1) = val;
    if (sample(t - 1, r - 1, c - 1, j - 1) > 0) {
      printf("sample c%d r%d t%d j%d = %d\n", c, r, t, j,
             sample(t - 1, r - 1, c - 1, j - 1));
    }
    nb_arrivals(c - 1, r - 1, t - 1) += val;
  }
  sample1_1.close();
  // cin.get();
  for (int c = 0; c < C; ++c) {
    for (int r = 0; r < R; ++r) {
      for (int t = 0; t < T; ++t) {
        empirical_lambda(c, r, t) =
            nb_arrivals(c, r, t) /
            static_cast<double>(nb_observations(c, r, t) * durations[t]);
      }
    }
  }

  for (int t = 0; t < T; ++t) {
    for (int r = 0; r < R; ++r) {
      int x = r % 10;
      int y = r / 10;
      double cx = x + 0.5;
      double cy = y + 0.5;
      for (int c = 0; c < C; ++c) {
        if ((!is_red[r] && t % 2 == 0) || (is_red[r] && t % 2 != 0)) {
          theoretical_lambda(c, r, t) = (constant_lambdas) ? 0.1 : (cx + cy);
        } else if ((!is_red[r] && t % 2 != 0) || (is_red[r] && t % 2 == 0)) {
          theoretical_lambda(c, r, t) =
              (constant_lambdas) ? 0.5 : 5 * (cx + cy);
        } else {
          std::cout << "ERROR: Unpredictable lambda case at t = ";
          std::cout << t << " and r = " << r << "(" << cx << "," << cy << ")\n";
          exit(1);
        }
      }
    }
  }

  vector<vector<int>> neighbors = vector<vector<int>>(R, vector<int>());
  xt::xarray<double> distance = GRB_INFINITY * xt::ones<double>({R, R});
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

    // if ((Yi + 1 <= n_y) && (Xi - 1 >= 1)) {
    //   neighbors[r].push_back(Yi * n_x + Xi - 1 - 1);
    // }

    // if ((Yi - 1 >= 1) && (Xi - 1) >= 1) {
    //   neighbors[r].push_back((Yi - 2) * n_x + Xi - 1 - 1);
    // }

    // if ((Yi + 1 <= n_y) && (Xi + 1) <= n_x) {
    //   neighbors[r].push_back(Yi * n_x + Xi + 1 - 1);
    // }

    // if ((Yi - 1 >= 1) && (Xi + 1) <= n_x) {
    //   neighbors[r].push_back((Yi - 2) * n_x + Xi + 1 - 1);
    // }
  }

  for (int r = 0; r < R; ++r) {
    for (int s : neighbors[r]) {
      distance(r, s) = 1;
    }
  }

  using laspated::Param;
  using laspated::projected_gradient_armijo_feasible;
  using laspated::RegularizedModel;

  Param param;
  double min_err = 1e100;
  int min_w = -1;
  double mean_emp = 0;
  vector<double> err_by_weight;
  xt::xarray<double> min_lambda;
  param.EPS = 0.001;
  param.sigma = 0.5;
  param.beta_bar = 1;
  param.max_iter = 30;
  param.upper_lambda = 1e3;
  param.lower_lambda = param.EPS;
  xt::xarray<double> x0 = param.EPS * xt::ones<double>({C, R, T});

  for (size_t i = 0; i < test_weights.size(); ++i) {
    if (constant_lambdas) {
      param.max_iter = 30;
      x0 = param.EPS * xt::ones<double>({C, R, T});
    }
    double w = test_weights[i];
    xt::xarray<double> alphas = w * neighbor_factor * xt::ones<double>({R, R});
    vector<double> weights = vector<double>(groups.size(), w);
    RegularizedModel model(nb_observations, nb_arrivals, durations, groups,
                           weights, alphas, distance, type_region, neighbors,
                           param);

    xt::xarray<double> x =
        projected_gradient_armijo_feasible<RegularizedModel>(model, param, x0);
    double err = model.average_rate_difference(theoretical_lambda, x);
    cout << "\tw = " << w << ", err = " << err << "\n";
    // cin.get();
    err_by_weight.push_back(err);
    mean_emp =
        model.average_rate_difference(theoretical_lambda, empirical_lambda);
    if (err < min_err) {
      min_err = err;
      min_w = i;
      min_lambda = x;
    }
  }

  double mean_cv = -1;
  double w_cv = -1;
  if (do_cross_validation) {
    // alphas and weights will be set by cross_validation
    xt::xarray<double> alphas = xt::zeros<double>({R, R});
    vector<double> weights = vector<double>(groups.size(), 0);

    param.EPS = 0.001;
    param.sigma = 0.5;
    param.beta_bar = 1;
    param.max_iter = 30;
    param.upper_lambda = 1e3;
    param.lower_lambda = param.EPS;
    xt::xarray<double> x0 = param.EPS * xt::ones<double>({C, R, T});

    RegularizedModel model(nb_observations, nb_arrivals, durations, groups,
                           weights, alphas, distance, type_region, neighbors,
                           param);
    param.cv_proportion = 0.2;
    std::vector<double> test_alphas = test_weights;
    auto result = cross_validation(param, model, sample, test_weights);
    mean_cv = model.average_rate_difference(theoretical_lambda, result.lambda);
    w_cv = result.weight;
  }

  return Result1{
      mean_emp,         min_err,   min_w, mean_cv, w_cv,
      err_by_weight,    C,         R,     T,       theoretical_lambda,
      empirical_lambda, min_lambda};
}