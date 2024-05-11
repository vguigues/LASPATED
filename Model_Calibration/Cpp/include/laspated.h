#ifndef LASPATED_H
#define LASPATED_H

#define USE_GUROBI 1
/* Source code for C++ calibration functions of paper titled: LASPATED: a
   Library for the Analysis of SPAtio-TEmporal Discrete data.

   Dependencies: gurobi, xtl, xtensor
*/
#ifdef USE_GUROBI
#include "gurobi_c++.h"
#else
#define GRB_INFINITY 1e100
#endif

#include <algorithm>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xview.hpp>

#include <xtensor/xarray.hpp>

namespace laspated {
class Param {
public:
  bool debug = false;
  // Path for custom data
  std::string data_folder = "calibration/NewHex_km";

  // Params for toy test
  std::pair<double, double> square_size{10, 10};
  std::pair<int, int> sub_regions{10, 10};
  int nb_groups = 2;
  int nb_weeks = 1000;
  bool full_neighbors = false;
  int neighbor_factor = 1;
  bool read_covariates = false;
  bool constant_lambdas = true;

  // Projected gradient params
  double EPS = 0.00001;
  double sigma = 0.5;
  double accuracy = 0.03;
  int max_iter = 30;
  double lower_lambda = 1e-6;
  double upper_lambda = 1e6;
  double beta_bar = 2;

  // Cross validation proportion
  double cv_proportion = 0.2;

  // Constructor methods
  Param() {}
  Param(const Param &p)
      : debug(p.debug), data_folder(p.data_folder), square_size(p.square_size),
        sub_regions(p.sub_regions), nb_groups(p.nb_groups),
        nb_weeks(p.nb_weeks), full_neighbors(p.full_neighbors),
        neighbor_factor(p.neighbor_factor), read_covariates(p.read_covariates),
        constant_lambdas(p.constant_lambdas), EPS(p.EPS), sigma(p.sigma),
        accuracy(p.accuracy), max_iter(p.max_iter),
        lower_lambda(p.lower_lambda), upper_lambda(p.upper_lambda),
        beta_bar(p.beta_bar), cv_proportion(p.cv_proportion) {}

  Param &operator=(const Param &p) {
    if (this == &p)
      return *this;

    debug = p.debug;
    // Path for custom data
    data_folder = p.data_folder;

    // Params for toy test
    square_size;
    sub_regions;
    nb_groups = p.nb_groups;
    nb_weeks = p.nb_weeks;
    full_neighbors = p.full_neighbors;
    neighbor_factor = p.neighbor_factor;
    read_covariates = p.read_covariates;
    constant_lambdas = p.constant_lambdas;

    // Projected gradient params
    EPS = p.EPS;
    sigma = p.sigma;
    accuracy = p.accuracy;
    max_iter = p.max_iter;
    lower_lambda = p.lower_lambda;
    upper_lambda = p.upper_lambda;
    beta_bar = p.beta_bar;

    // Cross validation proportion
    cv_proportion = p.cv_proportion;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &out, const Param &p) {
    out << "Parameters:\n";
    if (p.debug) {
      out << "DEBUG MODE\n";
    }
    out << "Data directory: " << p.data_folder << "\n";
    out << "Max Iter = " << p.max_iter << "\n";
    out << "EPS = " << p.EPS << "\n";
    out << "sigma = " << p.sigma << "\n";
    out << "beta_bar = " << p.beta_bar << "\n";
    return out;
  }
  ~Param();
};

class RegularizedModel {
public:
  Param &param;
  xt::xarray<int> nb_observations;
  xt::xarray<int> nb_arrivals;
  std::vector<double> durations;
  std::vector<std::vector<int>> groups;
  std::vector<double> weights;
  xt::xarray<double> alpha;
  xt::xarray<double> distance;
  std::vector<int> type_region;
  std::vector<std::vector<int>> neighbors;
  std::vector<int> which_group;

  ulong C, R, T;

  RegularizedModel(xt::xarray<int> &N, xt::xarray<int> &M,
                   std::vector<double> &a_durations,
                   std::vector<std::vector<int>> &a_groups,
                   std::vector<double> &a_weights, xt::xarray<double> &a_alphas,
                   xt::xarray<double> &a_distance, std::vector<int> &a_type,
                   std::vector<std::vector<int>> &a_neighbors, Param &a_param)
      : param(a_param) {
    if (N.dimension() != 3) { // N should be C,R,T
      std::cout << "Error: N has " << N.dimension()
                << " dimensions but should be 3\n";
      // fmt::print("Error: N has {} dimensions but should be 3. "
      //            "Problem was not set.\n",
      //            N.dimension());
      exit(1);
    }
    if (M.dimension() != 3) { // M also should be C,R,T
      std::cout << "Error: M has " << N.dimension()
                << " dimensions but should be 3\n";
      // fmt::print(
      //     "Error: M has {} dimensions but should be 3. Problem was not
      //     set.\n", M.dimension());
      exit(1);
    }

    C = N.shape(0);
    R = N.shape(1);
    T = N.shape(2);

    nb_observations = N;
    nb_arrivals = M;
    durations = a_durations;
    alpha = a_alphas;
    weights = a_weights;
    distance = a_distance;
    neighbors = a_neighbors;
    type_region = a_type;
    groups = a_groups;
    int nb_groups = groups.size();

    which_group = std::vector<int>(T, 0);
    for (int g = 0; g < nb_groups; ++g) {
      for (auto i : groups[g]) {
        which_group[i] = g;
      }
    }
  }

  double f(xt::xarray<double> &x) {
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

  xt::xarray<double> gradient(xt::xarray<double> &x) {

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
            if (type_region[r] == type_region[s]) {
              grad_component += 2 * alpha(r, s) * (x(c, r, t) - x(c, s, t)) /
                                (distance(r, s));
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

  xt::xarray<double> projection(xt::xarray<double> &x) {
    xt::xarray<double> z = xt::zeros<double>(x.shape());
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          z(c, r, t) = std::max(x(c, r, t), param.lower_lambda);
        }
      }
    }

    return z;
  }

  bool is_feasible(xt::xarray<double> &x) {
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          if (x(c, r, t) < param.lower_lambda) {
            return false;
          }
        }
      }
    }
    return true;
  }

  double get_rhs(xt::xarray<double> &grad, xt::xarray<double> &dir) {
    double sum = 0;
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          sum += grad(c, r, t) * dir(c, r, t);
        }
      }
    }
    return sum;
  }

  double get_lower_bound(xt::xarray<double> &x, xt::xarray<double> &grad) {
    double change = 0.0;
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          if (grad(c, r, t) < 0.0) {
            change += grad(c, r, t) * (param.upper_lambda - x(c, r, t));
          } else if (grad(c, r, t) > 0.0) {
            change += grad(c, r, t) * (param.lower_lambda - x(c, r, t));
          }
        }
      }
    }
    return change;
  }
};

#ifdef USE_GUROBI
class CovariatesModel {
public:
  GRBEnv env;
  Param &param;
  xt::xarray<int> nb_observations;
  xt::xarray<int> nb_arrivals;
  xt::xarray<double> regressors;
  ulong C, D, T, R, nb_regressors;

  CovariatesModel(xt::xarray<int> &N, xt::xarray<int> &M,
                  xt::xarray<double> &reg, Param &param)
      : env(), param(param) {
    if (N.dimension() != 4) { // N should be C,D,T,R
      std::cout << "Error: N has " << N.dimension()
                << " dimensions but should be 4\n";
      // fmt::print(
      //     "Error: N has {} dimensions but should be 4. Problem was not
      //     set.\n", N.dimension());
      exit(1);
    }
    if (M.dimension() != 4) { // M also should be C,D,T,R
      std::cout << "Error: M has " << N.dimension()
                << " dimensions but should be 4\n";
      // fmt::print(
      //     "Error: M has {} dimensions but should be 4. Problem was not
      //     set.\n", M.dimension());
      exit(1);
    }

    if (reg.dimension() != 2) { // R should be nb_regressors,R
      std::cout << "Error: regressors has " << reg.dimension()
                << "dimensions but should be 2.\n";
      // fmt::print("Error: regressor array has {} dimensons but should be 2. "
      //            "Problem was not set.\n",
      //            reg.dimension());
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
  }

  double f(xt::xarray<double> &x) {
    double f = 0;

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            double rates = 0;
            for (int j = 0; j < nb_regressors; ++j) {
              rates += x(c, d, t, j) * regressors(j, r);
            }
            f += nb_observations(c, d, t, r) * rates -
                 nb_arrivals(c, d, t, r) * log(rates);
          }
        }
      }
    }
    return f;
  }

  xt::xarray<double> gradient(xt::xarray<double> &x) {
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
            }
          }
        }
      }
    }
    return gradient;
  }

  bool is_feasible(xt::xarray<double> &x) {
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            double rates = 0;
            for (int j = 0; j < nb_regressors; ++j) {
              rates += x(c, d, t, j) * regressors(j, r);
            }
            if (rates < param.EPS) {
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  double get_rhs(xt::xarray<double> &grad, xt::xarray<double> &dir) {
    double sum = 0;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            sum += grad(c, d, t, j) * dir(c, d, t, j);
          }
        }
      }
    }
    return sum;
  }

  double get_lower_bound(xt::xarray<double> &x, xt::xarray<double> &grad) {
    GRBModel model(env);
    std::stringstream name;
    xt::xarray<GRBVar> y({C, D, T, nb_regressors});

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            name << "y_" << c << "_" << d << "_" << t << "_" << j;
            double ub = (j == 0) ? 1 : param.upper_lambda;
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
            // fmt::print("x(c{},d{},t{},j{}) = {:.7f}, grad = {:.3f}\n",
            // c,d,t,j, 	x(c,d,t,j),grad(c,d,t,j));
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
            model.addConstr(con1, GRB_GREATER_EQUAL, param.EPS, name.str());
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
      std::cout << "Status = " << status << "\n";
      // fmt::print("Status = {}\n", status);
      model.write("lower_bound_regressors.lp");
      std::cin.get();
      return GRB_INFINITY;
    }
  }
};
#endif

template <typename Model>
xt::xarray<double> projected_gradient_armijo_feasible(Model &model,
                                                      Param &param,
                                                      xt::xarray<double> &x) {
  int k = 0;
  double beta_k = param.beta_bar;
  double accuracy = param.accuracy;
  double eps = param.EPS;
  double upper_lambda = param.upper_lambda;
  double upper_bound = GRB_INFINITY;
  double sigma = param.sigma;
  int max_iter = param.max_iter;

  int j = 0;

  if (!model.is_feasible(x)) {
    x = model.projection(x);
  }
  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());

  double fold = model.f(x);
  xt::xarray<double> gradient = model.gradient(x);
  std::vector<double> f_val;

  while (k < max_iter) {
    x_aux = x - beta_k * gradient;
    z = model.projection(x_aux);
    if (!model.is_feasible(z)) {
      std::cout << "Infeasible z\n";
      // fmt::print("Infeasible z\n");
      exit(1);
    }
    diff_aux = x - z;
    double rhs = model.get_rhs(gradient, diff_aux);
    double f = model.f(z);
    if (rhs > 0.0 && f > fold - sigma * rhs) {
      bool stop = false;
      z_aux = xt::zeros<double>(x.shape());
      double this_pow = 1.0;
      int count = 0;
      double best_pow = -1.0;
      double best_val = f;
      while (!stop) {
        z_aux = x + (1.0 / this_pow) * (z - x);
        if (!model.is_feasible(z_aux)) {
          std::cout << "Infeasible z_aux\n";
          // fmt::print("Infeasible z_aux\n");
          std::cin.get();
        }
        f = model.f(z_aux);
        if (f < best_val) {
          best_val = f;
          best_pow = this_pow;
        }
        // fmt::print("\tj = {}, f = {}, 1/2j = {}\n", count, f, (1.0 /
        // this_pow));
        if (f <= fold - (sigma / this_pow) * rhs) {
          stop = true;
        } else {
          this_pow *= 2;
        }
        ++count;
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
    gradient = model.f(x);
    upper_bound = std::min(f, upper_bound);
    double lb = model.get_lower_bound(x, gradient);
    double lower_bound = std::max(lower_bound, f + lb);
    double gap = abs((upper_bound - lower_bound) / upper_bound);
    if (gap < accuracy) {
      break;
    }
    ++k;
  }

  return x;
}

typedef struct {
  double cpu_time;
  double weight;
  xt::xarray<double> lambda;
} CrossValidationResult;

CrossValidationResult cross_validation(Param &param, RegularizedModel &model,
                                       xt::xarray<int> &sample,
                                       std::vector<double> &alphas,
                                       std::vector<double> &group_weights) {
  auto t0 = std::chrono::high_resolution_clock::now();
  int nb_observations_total = sample.shape(3);
  int nb_groups = model.groups.size();
  double min_loss = GRB_INFINITY;
  double cpu_time = 0;
  int nb_in_block = floor(nb_observations_total * param.cv_proportion);
  xt::xarray<int> initial_nb_obs = model.nb_observations;
  xt::xarray<int> initial_nb_arrivals = model.nb_arrivals;
  // fmt::print("cross validation weights.size = {}\n", group_weights.size());
  double best_alpha = GRB_INFINITY;
  double best_weight = GRB_INFINITY;
  // fmt::print("Running cross validation with proportion = {} and weights =
  // {}\n", proportion, group_weights);
  for (int index_alpha = 0; index_alpha < alphas.size(); ++index_alpha) {
    double likelihood = 0;
    model.alpha = alphas[index_alpha] * xt::ones<double>({model.R, model.R});
    model.weights =
        std::vector<double>(model.groups.size(), group_weights[index_alpha]);
    // fmt::print("Testing weight = {}\n", group_weights[index_alpha]);
    for (int index_cross = 0; index_cross < floor(1 / param.cv_proportion);
         ++index_cross) {
      xt::xarray<int> nb_observations_current =
          xt::zeros<int>({model.C, model.R, model.T});
      xt::xarray<int> nb_calls_current =
          xt::zeros<int>({model.C, model.R, model.T});
      for (int index = index_cross * nb_in_block;
           index < (index_cross + 1) * nb_in_block; ++index) {
        for (int c = 0; c < model.C; ++c) {
          for (int r = 0; r < model.R; ++r) {
            for (int t = 0; t < model.T; ++t) {
              ++nb_observations_current(c, r, t);
              nb_calls_current(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      xt::xarray<double> x =
          param.EPS * xt::ones<double>({model.C, model.R, model.T});
      model.nb_observations = nb_observations_current;
      model.nb_arrivals = nb_calls_current;
      auto f_val =
          projected_gradient_armijo_feasible<RegularizedModel>(model, param, x);
      xt::xarray<int> nb_calls_remaining =
          xt::zeros<int>({model.C, model.R, model.T});
      for (int index = 0; index < index_cross * nb_in_block; ++index) {
        for (int c = 0; c < model.C; ++c) {
          for (int r = 0; r < model.R; ++r) {
            for (int t = 0; t < model.T; ++t) {
              nb_calls_remaining(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      for (int index = (index_cross + 1) * nb_in_block;
           index < nb_observations_total; ++index) {
        for (int c = 0; c < model.C; ++c) {
          for (int r = 0; r < model.R; ++r) {
            for (int t = 0; t < model.T; ++t) {
              nb_calls_remaining(c, r, t) += sample(t, r, c, index);
            }
          }
        }
      }
      double f = 0;
      for (int c = 0; c < model.C; ++c) {
        for (int r = 0; r < model.R; ++r) {
          for (int t = 0; t < model.T; ++t) {
            double current_lambda = x(c, r, t);
            f += (nb_observations_total - nb_in_block) * current_lambda *
                     model.durations[t] -
                 nb_calls_remaining(c, r, t) * log(current_lambda);
          }
        }
      }
      // fmt::print("\tf test set {} = {}, first f_val = {:.1f}, last f_val =
      // {:.1f}\n",index_cross, f, f_val[0], f_val.back());
      likelihood += f;
    }
    likelihood = likelihood / floor(1 / param.cv_proportion);
    // fmt::print("Likelihood current_alpha {} = {}\n", alphas[index_alpha],
    // likelihood);
    if (likelihood < min_loss) {
      min_loss = likelihood;
      best_alpha = alphas[index_alpha];
      best_weight = group_weights[index_alpha];
    }
  }
  // alpha = best_alpha;
  // weights = std::vector<double>(groups.size(), best_weight);
  model.alpha = best_alpha * xt::ones<double>({model.R, model.R});
  model.weights = std::vector<double>(model.groups.size(), best_weight);
  xt::xarray<int> nb_observations_current =
      xt::zeros<int>({model.C, model.R, model.T});
  xt::xarray<int> nb_calls_current =
      xt::zeros<int>({model.C, model.R, model.T});
  for (int c = 0; c < model.C; ++c) {
    for (int r = 0; r < model.R; ++r) {
      for (int t = 0; t < model.T; ++t) {
        for (int index = 0; index < nb_observations_total; ++index) {
          nb_calls_current(c, r, t) += sample(t, r, c, index);
          ++nb_observations_current(c, r, t);
        }
      }
    }
  }

  model.nb_observations = nb_observations_current;
  model.nb_arrivals = nb_calls_current;
  xt::xarray<double> x =
      param.EPS * xt::ones<double>({model.C, model.R, model.T});
  auto f_val =
      projected_gradient_armijo_feasible<RegularizedModel>(model, param, x);
  auto dt = std::chrono::high_resolution_clock::now();
  cpu_time = std::chrono::duration_cast<std::chrono::seconds>(dt - t0).count();
  model.nb_observations = initial_nb_obs;
  model.nb_arrivals = initial_nb_arrivals;
  // fmt::print("best_weight = {}\n", best_weight);
  return {cpu_time, best_weight, x};
}

} // namespace laspated

#endif
