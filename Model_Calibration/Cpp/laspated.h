#ifndef LASPATED_H
#define LASPATED_H

/* Source code for C++ calibration functions of paper titled: LASPATED: a
   Library for the Analysis of SPAtio-TEmporal Discrete data.

   Dependencies: gurobi, xtl and xtensor (we provide xtl and xtensor in the
   repository)

   Classes:

   Param: parameters for the projected gradient and cross validation functions.
      Attributes:
        EPS > 0: tolerance
        sigma in (0,1): parameter used in the line search of the projected
            gradient method
        gap in (0,1): gap used as stopping criterion for the projected
            gradient method
        max_iter > 0: maximum number of iterations for the projected
            gradient method
        lower_lambda: lower bound for decision variables in the projected
            gradient method
        upper_lambda: upper bound for decision variables in the projected
            gradient method
        beta_bar > 0: initial step size for the projected gradient
            method
        cv_proportion: proportion of samples used in each sub-iteration of the
            cross validation.

    AppParameters: Application-wide parameters. This object contains all
        attributes from Param and also the following:
      Attributes:
        model_type: The type of model (reg | no_reg)
        method: The calibration method (calibration | cross_validation)
        algorithm: The calibration algorithm (feasible | boundary)
        info_file: The general information file
        arrivals_file: The sample arrivals file
        neighbors_file: The zones neighborhood file
        alpha_regions_file: The spatial regularization matrix file
        time_groups_file:  The temporal regularization file
        duration:  The duration of each period
        cv_weights_file: The weights used in cross validation
        output_file: The intensities output file


    RegularizedModel: Model without covariates.
      Attributes:
        param: reference for parameter object.
        nb_observations: number of observations for each class-space-time index.
        nb_arrivals: number of events for each class-space-time index.
        durations: duration, in hours, for each time index.
        groups: description of time groups. groups[i][j] is the j-th time index
            in the i-th group.
        weights: time regularization weight for each time group.
        alpha: space regularization weight matrix of subregions.
        distance: distance matrix of subregions.
        type_region: space groups. type_region[r] is the type of subregion r.
        neighbors: neighborhood description. neighbors[i][j] is the j-th
            neighbor of subregion i.
        which_group: groups of each time index.
        C: number of events classes.
        R: number of subregions.
        T: number of time periods.

      Methods:
        f(x): objective function at x
        gradient(x): gradient at x
        projection(x): projection of x in the feasible set
        is_feasible(x): returns true if x is in the feasible set
        get_rhs(grad, dir): returns the directional derivative given by gradient
            grad and direction dir.
        lower_bound(x, grad): returns objective function lower bound given x and
            gradient grad.
        average_rate_difference(x1,x2): returns average difference in rates
            between x1 and x2.
    CovariatesModel: Model with covariates
      Attributes:
        env: Gurobi environment object
        param: reference for parameter object.
        nb_observations: number of observations for each class-space-time index.
        nb_arrivals: number of events for each class-space-time index.
        regressors: regressors array. regressor(j,r) is the j-th regressor for
              subregion r.
        C: number of events classes.
        D: number of days.
        T: number of time periods.
        R: number of subregions.
        nb_regressors: number of regressors.

      Methods:
        f(x): objective function at x
        gradient(x): gradient at x
        projection(x): projection of x in the feasible set
        is_feasible(x): returns true if x is in the feasible set
        get_rhs(grad, dir): returns the directional derivative given by gradient
            grad and direction dir.
        lower_bound(x, grad): returns objective function lower bound given x and
            gradient grad.
        average_rate_difference(x1,x2): returns average difference in rates
            between x1 and x2.
    CrossValidationResult: result struct for cross validation
      Attributes:
        wall_time: time spent, in seconds, to run cross validation
        weight: best weight found in cross validation
        lambda: solution given by best weight found in cross validation
    Free functions:

    projected_gradient_armijo_feasible<Model>(model, param, x0): runs projected
        gradient with initial solution x0 using model (<Model> can be
        <RegularizedModel> or <CovariatesModel>) and param objects.

    cross_validation(param, regularized_model, sample, alphas, group_weights):
        runs cross validation given parameters, sample array, weights alphas and
        corresponding time group weights. alphas and group_weights must be the
        same size.
*/
#if USE_GUROBI == 1
#include "gurobi_c++.h"
#else
#define GRB_INFINITY 1e100
#endif

#include <algorithm>
#include <boost/program_options.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <xtensor/xarray.hpp>

namespace po = boost::program_options;

namespace laspated {

class AppParameters {
 public:
  // Projected gradient params
  double EPS;           // Epsilon for feasibility/convergence
  double sigma;         // Sigma for armijo step
  double gap;           // Stopping criterion gap.
  int max_iter;         // Max number of iterations.
  double lower_lambda;  // lower bound on decision variables
  double upper_lambda;  // upper bound on decision variables
  double beta_bar;      // Initial step size in projected gradient

  // Cross validation proportion
  double cv_proportion = 0.2;  // Cross validation proportion

  // Model parameters
  std::string model_type;          // reg | no_reg
  std::string method;              // calibration | cross_validation
  std::string algorithm;           // feasible | boundary
  std::string info_file;           // path to info.dat
  std::string arrivals_file;       // path to arrivals.dat
  std::string neighbors_file;      // path to neighbors.dat
  std::string alpha_regions_file;  // path to alpha matrix
  std::string time_groups_file;    // path to time groups description
  std::string durations_file;
  std::string cv_weights_file;  // path to cross_validation weights
  std::string output_file;      // Path to results file
};

class Param {
 public:
  // Projected gradient params
  double EPS = 1e-5;
  double sigma = 0.5;
  double gap = 0.01;
  int max_iter = 30;
  double lower_lambda = 1e-6;
  double upper_lambda = 1.0;
  double beta_bar = 2.0;

  // Cross validation proportion
  double cv_proportion = 0.2;

  // Constructor methods
  Param() {}
  Param(AppParameters &app_params)
      : EPS(app_params.EPS),
        sigma(app_params.sigma),
        gap(app_params.gap),
        max_iter(app_params.max_iter),
        lower_lambda(app_params.lower_lambda),
        upper_lambda(app_params.upper_lambda),
        beta_bar(app_params.beta_bar),
        cv_proportion(app_params.cv_proportion) {}
  Param(const Param &p)
      : EPS(p.EPS),
        sigma(p.sigma),
        gap(p.gap),
        max_iter(p.max_iter),
        lower_lambda(p.lower_lambda),
        upper_lambda(p.upper_lambda),
        beta_bar(p.beta_bar),
        cv_proportion(p.cv_proportion) {}

  Param &operator=(const Param &p) {
    if (this == &p) return *this;

    // Projected gradient params
    EPS = p.EPS;
    sigma = p.sigma;
    gap = p.gap;
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

    out << "Max Iter = " << p.max_iter << "\n";
    out << "EPS = " << p.EPS << "\n";
    out << "sigma = " << p.sigma << "\n";
    out << "beta_bar = " << p.beta_bar << "\n";
    return out;
  }
  ~Param() = default;
};

class RegularizedModel {
 public:
  Param &param;
  std::string name = "Regularized";
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
    if (N.dimension() != 3) {  // N should be C,R,T
      std::cout << "Error: N has " << N.dimension()
                << " dimensions but should be 3\n";
      // fmt::print("Error: N has {} dimensions but should be 3. "
      //            "Problem was not set.\n",
      //            N.dimension());
      exit(1);
    }
    if (M.dimension() != 3) {  // M also should be C,R,T
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
    double obj = 0;

    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          double current_lambda = x(c, r, t);
          obj += nb_observations(c, r, t) * current_lambda * durations[t] -
                 nb_arrivals(c, r, t) * log(current_lambda * durations[t]);

          for (int s : neighbors[r]) {
            if (true /*type_region[r] == type_region[s] */) {
              obj += (0.5 * alpha(r, s)) * nb_observations(c, r, t) *
                     nb_observations(c, s, t) *
                     pow(x(c, r, t) - x(c, s, t), 2) / distance(r, s);
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
                obj += (0.5 * weights[grindex]) * nb_observations(c, r, t) *
                       nb_observations(c, r, tp) *
                       (pow(x(c, r, t) - x(c, r, tp), 2));
              }
            }
          }
        }
      }
    }

    return obj;
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
            if (true /*type_region[r] == type_region[s]*/) {
              grad_component += 2 * alpha(r, s) * nb_observations(c, r, t) *
                                nb_observations(c, s, t) *
                                (x(c, r, t) - x(c, s, t)) / (distance(r, s));
            }
          }
          gradient(c, r, t) = grad_component;
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
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
          if (x(c, r, t) < param.lower_lambda - param.EPS) {
            std::cout << "x(" << c << "," << r << "," << t << ") = ";
            std::cout << x(c, r, t) << " is less than lower_bound "
                      << param.lower_lambda << "\n";
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

  double average_rate_difference(xt::xarray<double> &x1,
                                 xt::xarray<double> &x2) {
    using namespace std;
    vector<double> difference_l2;
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
          difference_l2.push_back(abs(x1(c, r, t) - x2(c, r, t)) / x1(c, r, t));
        }
      }
    }
    // fmt::print("Error: {}\n", difference_l2);
    return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) /
           difference_l2.size();
  }
};

#if USE_GUROBI == 1
class CovariatesModel {
 public:
  GRBEnv env;
  Param &param;
  std::string name = "Covariates";
  xt::xarray<int> nb_observations;
  xt::xarray<int> nb_arrivals;
  xt::xarray<double> regressors;
  ulong C, D, T, R, nb_regressors;

  CovariatesModel(xt::xarray<int> &N, xt::xarray<int> &M,
                  xt::xarray<double> &reg, Param &param)
      : param(param) {
    if (N.dimension() != 4) {  // N should be C,D,T,R
      std::cout << "Error: N has " << N.dimension()
                << " dimensions but should be 4\n";
      // fmt::print(
      //     "Error: N has {} dimensions but should be 4. Problem was not
      //     set.\n", N.dimension());
      exit(1);
    }
    if (M.dimension() != 4) {  // M also should be C,D,T,R
      std::cout << "Error: M has " << N.dimension()
                << " dimensions but should be 4\n";
      // fmt::print(
      //     "Error: M has {} dimensions but should be 4. Problem was not
      //     set.\n", M.dimension());
      exit(1);
    }

    if (reg.dimension() != 2) {  // R should be nb_regressors,R
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
            // printf("rates %d %d %d %d\n", c, d, t, r);
            for (int j = 0; j < nb_regressors; ++j) {
              rates += x(c, d, t, j) * regressors(j, r);
              // printf("\tj = %d, x = %f, reg = %f\n", j, x(c, d, t, j),
              //        regressors(j, r));
            }
            // printf("rates = %f", rates);
            // std::cin.get();
            f += nb_observations(c, d, t, r) * rates -
                 nb_arrivals(c, d, t, r) * log(rates);
          }
        }
      }
    }

    return f;
  }

  xt::xarray<double> gradient(xt::xarray<double> &x) {
    using namespace std;
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

    xt::xarray<double> v1 = xt::zeros<double>(x.shape());
    xt::xarray<double> v2 = xt::zeros<double>(x.shape());
    xt::xarray<double> rates = xt::zeros<double>({C, D, T, R});
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            for (int r = 0; r < R; ++r) {
              rates(c, d, t, r) += x(c, d, t, j) * regressors(j, r);
            }
          }
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            vector<double> v1_aux;
            v1_aux.reserve(R);
            for (int r = 0; r < R; ++r) {
              v1_aux.push_back(nb_observations(c, d, t, r) * regressors(j, r));
            }
            sort(v1_aux.begin(), v1_aux.end());
            for (int r = 0; r < R; ++r) {
              v1(c, d, t, j) += v1_aux[r];
            }
            vector<double> v2_aux;
            v2_aux.reserve(R);
            for (int r = 0; r < R; ++r) {
              if (rates(c, d, t, r) == 0.0) {
                cout << "ERROR, rates = 0 in gradient.\n";
                cin.get();
              } else {
                v1_aux.push_back(nb_arrivals(c, d, t, r) * regressors(j, r) /
                                 rates(c, d, t, r));
              }
            }
            sort(v2_aux.begin(), v2_aux.end());
            for (int r = 0; r < R; ++r) {
              v2(c, d, t, j) += v2_aux[r];
            }
          }
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            gradient(c, d, t, j) = v1(c, d, t, j) - v2(c, d, t, j);
          }
        }
      }
    }

    return gradient;
  }

  xt::xarray<double> projection(xt::xarray<double> &x) {
    using namespace std;
    xt::xarray<GRBVar> y({C, D, T, nb_regressors});
    GRBModel model(env);
    stringstream name;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            name << "y_" << c << "_" << d << "_" << t << "_" << j;
            double ub = (j == 0) ? 1 : param.upper_lambda;
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
            if (x(c, d, t, j) != x(c, d, t, j)) {
              printf("x(c%d,d%d,t%d,j%d) = %f\n", c, d, t, j, x(c, d, t, j));
              cin.get();
            }
            obj += 0.5 * y(c, d, t, j) * y(c, d, t, j) -
                   x(c, d, t, j) * y(c, d, t, j);
          }
        }
      }
    }
    // cin.get();
    try {
      model.setObjective(obj, GRB_MINIMIZE);
    } catch (GRBException &ex) {
      cout << ex.getMessage() << "\n";
      cin.get();
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
            model.addConstr(con1, GRB_GREATER_EQUAL, param.EPS, name.str());
            name.str("");
            con1 = 0;
          }
        }
      }
    }

    model.update();
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
            }
          }
        }
      }
      // cin.get();
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }
    return y_val;
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
    using namespace std;
    vector<double> aux1;
    vector<double> aux2;
    double sum = 0;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < nb_regressors; ++j) {
            // if (grad(c, d, t, j) * dir(c, d, t, j) >= 0.0) {
            //   aux1.push_back(grad(c, d, t, j) * dir(c, d, t, j));
            // } else {
            //   aux2.push_back(-1 * grad(c, d, t, j) * dir(c, d, t, j));
            // }
            sum += grad(c, d, t, j) * dir(c, d, t, j);
          }
        }
      }
    }

    return sum;

    sort(aux1.begin(), aux1.end());
    sort(aux2.begin(), aux2.end());
    double rhs = 0.0;
    for (auto val : aux1) {
      rhs += val;
    }
    double rhs2 = 0.0;
    for (auto val : aux2) {
      rhs2 += val;
    }

    return rhs - rhs2;
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
      return model.get(GRB_DoubleAttr_ObjVal);
    } else {
      std::cout << "Error lower bound problem solved with status = " << status
                << "\n";
      model.write("lower_bound_regressors.lp");
      std::cout << "Wrote model at lower_bound.lp\n";
      exit(status);
      return GRB_INFINITY;
    }
  }

  double average_rate_difference(xt::xarray<double> &x1,
                                 xt::xarray<double> &x2) {
    using namespace std;
    vector<double> difference_l2;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            double rate1 = 0;
            double rate2 = 0;
            for (int j = 0; j < nb_regressors; ++j) {
              rate1 += x1(c, d, t, j) * regressors(j, r);
              rate2 += x2(c, d, t, j) * regressors(j, r);
            }
            difference_l2.push_back(abs(rate1 - rate2) / rate1);
          }
        }
      }
    }

    return accumulate(difference_l2.begin(), difference_l2.end(), 0.0) /
           difference_l2.size();
  }

  ~CovariatesModel() = default;
};
#endif

template <typename Model>
xt::xarray<double> projected_gradient_armijo_feasible(Model &model,
                                                      Param &param,
                                                      xt::xarray<double> &x) {
  if (model.name == "Regularized") {
    return regularized(model, param, x);
  } else if ((model.name == "Covariates")) {
    return covariates(model, param, x);
  }

  int k = 0;
  double beta_k = param.beta_bar;
  double gap = param.gap;
  double eps = param.EPS;
  double upper_lambda = param.upper_lambda;
  double upper_bound = GRB_INFINITY;
  double lower_bound = -GRB_INFINITY;
  double sigma = param.sigma;
  int max_iter = param.max_iter;

  int j = 0;
  bool is_feasible = model.is_feasible(x);
  if (!is_feasible) {
    x = model.projection(x);
  }
  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());

  double fold = model.f(x);
  xt::xarray<double> gradient = model.gradient(x);
  std::vector<double> f_val;
  double this_gap = 1e100;
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
    // std::cout << "\tk = " << k << ", fold = " << fold << ", beta_k = " <<
    // beta_k
    //           << ", rhs = " << rhs << ", f(z) = " << f << "\n";
    if (rhs > 0.0 && f > fold - sigma * rhs) {
      bool stop = false;
      // z_aux = xt::zeros<double>(x.shape());
      double this_pow = 1.0;
      int count = 0;
      double best_pow = -1.0;
      double best_val = 1e100;
      while (!stop) {
        z_aux = x + (1.0 / this_pow) * (z - x);
        if (!model.is_feasible(z_aux)) {
          std::cout << "Error: Infeasible z_aux\n";
          exit(100);
        }
        f = model.f(z_aux);
        if (f < best_val) {
          best_val = f;
          best_pow = this_pow;
        }
        // std::cout << "\t\tj = " << count << ", f = " << f
        //           << ", 1/(2j) = " << 1 / this_pow
        //           << "| best_pow = " << best_pow << "\n";
        if (f <= fold - (sigma / this_pow) * rhs) {
          stop = true;
        } else {
          this_pow *= 2;
        }
        ++count;
      }
      // printf("best_pow = %f\n", best_pow);
      // std::cin.get();
      if (best_pow < 0.0) {
        f = fold;
        beta_k *= 2.0;
      } else {
        f = best_val;
        beta_k = 2.0 / best_pow;
        x = x + (1.0 / best_pow) * (z - x);
      }
      f_val.push_back(f);
    } else {
      if (rhs > 0.0) {
        x = z;
      } else {
        f = fold;
      }
      beta_k *= 2.0;
    }

    fold = f;
    gradient = model.gradient(x);
    upper_bound = std::min(f, upper_bound);
    double lb = model.get_lower_bound(x, gradient);
    lower_bound = std::max(lower_bound, f + lb);
    // std::cout << "\tfold = " << fold << ", lower_bound = " << lower_bound
    //           << ", upper = " << upper_bound
    //           << " | upper - lower = " << abs(upper_bound - lower_bound)
    //           << " abs(ub) = " << abs(upper_bound)
    //           << " gap = " << abs(upper_bound - lower_bound) /
    //           abs(upper_bound)
    //           << "\n";
    this_gap = fabs((upper_bound / 10000) - (lower_bound / 10000)) /
               fabs(upper_bound / 10000);
    // printf("\tgap = %.3f\n", gap);
    if (gap < 0.0) {
      printf("Negative gap = %f", gap);
      std::cout << "\tfold = " << fold << ", lower_bound = " << lower_bound
                << ", upper = " << upper_bound
                << " | upper - lower = " << abs(upper_bound - lower_bound)
                << " abs(ub) = " << abs(upper_bound) << " gap = " << gap
                << "\n";
      std::cin.get();
    }
    if (this_gap < gap) {
      // printf("Gap closed!\n");
      break;
    }
    ++k;
  }
  // printf("k = %d, gap = %f\n", k, gap);
  return x;
}

template <typename Model>
xt::xarray<double> regularized(Model &model, Param &param,
                               xt::xarray<double> &x) {
  using namespace std;
  int k = 0;
  vector<double> f_val;
  double b_param = 2;
  double beta_k = b_param;
  double eps = param.EPS;
  double upper_lambda = param.upper_lambda;
  double upper_bound = GRB_INFINITY;
  double lower_bound = -GRB_INFINITY;
  double sigma = param.sigma;
  int max_iter = param.max_iter;

  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> gradient = xt::zeros<double>(x.shape());

  double fold = model.f(x);
  gradient = model.gradient(x);
  while (k < max_iter) {
    x_aux = x - beta_k * gradient;
    z = model.projection(x_aux);

    int j = 1;
    diff_aux = x - z;
    double f = model.f(z);
    double rhs = model.get_rhs(gradient, diff_aux);
    // double diff = f - fold + sigma * rhs;
    // printf("\tk = %d, fold = %f, f = %f, rhs = %f diff = %f\n", k, fold, f,
    // rhs,
    //        diff);
    if (f > fold - sigma * rhs) {
      bool stop = false;
      while (!stop) {
        z_aux = x + (1 / pow(2.0, j)) * (z - x);
        f = model.f(z_aux);
        // printf("\t\tj = %d, f = %f z_aux = (%f,%f)\n", j, f, z_aux(0),
        //        z_aux(1));
        if (f <= fold - (sigma / pow(2.0, j)) * rhs) {
          stop = true;
        } else {
          ++j;
        }
      }
      f_val.push_back(f);
      x = z_aux;
      beta_k = b_param / pow(2.0, j - 1);
    } else {
      x = z;
      beta_k *= 2;
    }
    fold = f;
    gradient = model.gradient(x);
    ++k;
    // std::cin.get();
  }
  return x;
}

template <typename Model>
xt::xarray<double> covariates(Model &model, Param &param,
                              xt::xarray<double> &x) {
  using namespace std;
  int k = 0;
  double b_param = param.beta_bar;
  double beta_k = b_param;
  double eps = param.EPS;
  double upper_lambda = param.upper_lambda;
  double sigma = param.sigma;
  vector<double> f_val;
  int max_iter = param.max_iter;
  x = model.projection(x);

  if (!model.is_feasible(x)) {
    printf("Projected x is not feasible!\n");
    std::cin.get();
  }

  xt::xarray<double> z = xt::zeros<double>(x.shape());
  xt::xarray<double> z_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
  xt::xarray<double> x_aux = xt::zeros<double>(x.shape());

  while (k < max_iter) {
    double fold = model.f(x);
    xt::xarray<double> gradient = model.gradient(x);
    xt::xarray<double> x_aux = x - beta_k * gradient;
    z = model.projection(x_aux);
    bool stop = false;
    int j = 0;
    diff_aux = x - z;
    double f = model.f(z);
    double rhs = model.get_rhs(gradient, diff_aux);
    // printf("k = %d, fold = %.8f, f = %.8f, rhs = %f\n", k, fold, f, rhs);
    if (f > fold - sigma * rhs + param.EPS) {
      bool stop = false;
      while (!stop) {
        z_aux = x + (1 / pow(2.0, j)) * (z - x);
        f = model.f(z_aux);
        // printf("\tj = %d, f = %f, fold = %f, sigma = %f, rhs = %f\n", j, f,
        //        fold, sigma, rhs);
        if (f <= fold - (sigma / pow(2.0, j)) * rhs) {
          stop = true;
        } else {
          ++j;
        }
      }
      f_val.push_back(f);
      x = z_aux;
      beta_k = b_param / pow(2.0, j - 1);
    } else {
      x = z;
      beta_k *= 2;
    }
    ++k;
  }
  // cin.get();
  return x;
}

template <typename Model>
xt::xarray<double> projected_gradient_armijo_boundary(Model &model,
                                                      Param &param,
                                                      xt::xarray<double> &x) {
  int k = 0;
  double beta_bar = param.beta_bar;
  // x = model.projection(x);
  while (k < param.max_iter) {
    double fold = model.f(x);
    xt::xarray<double> gradient = model.gradient(x);
    bool stop = false;
    int j = 0;
    xt::xarray<double> z = xt::zeros<double>(x.shape());
    xt::xarray<double> x_aux = xt::zeros<double>(x.shape());
    xt::xarray<double> diff_aux = xt::zeros<double>(x.shape());
    printf("k = %d, fold = %f\n", k, fold);
    while (!stop) {
      x_aux = x - (beta_bar / pow(2.0, j)) * gradient;
      z = model.projection(x_aux);
      double f = model.f(z);
      diff_aux = x - z;
      double rhs = model.get_rhs(gradient, diff_aux);
      printf("\tj = %d, f = %f, rhs = %f, sigma = %f\n", j, f, rhs,
             param.sigma);
      if (f <= fold - param.sigma * rhs) {
        stop = true;
      } else {
        ++j;
      }
    }
    x = z;
    ++k;
  }

  return x;
}

typedef struct {
  double wall_time;
  double weight;
  xt::xarray<double> lambda;
} CrossValidationResult;

CrossValidationResult cross_validation(Param &param, RegularizedModel &model,
                                       xt::xarray<int> &sample,
                                       std::vector<double> &group_weights) {
  if (sample.dimension() != 4 || sample.shape(0) != model.T ||
      sample.shape(1) != model.R || sample.shape(2) != model.C) {
    std::cout << "Incorrect sample shape. Sample must be shaped {T,R,C,K}.\n";
    exit(1);
  }
  std::vector<double> alphas = group_weights;
  auto t0 = std::chrono::high_resolution_clock::now();
  int nb_observations_total = sample.shape(3);
  int nb_groups = model.groups.size();
  double min_loss = GRB_INFINITY;
  double wall_time = 0;
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
    // std::cout << "Testing weight = " << group_weights[index_alpha] << "\n";
    for (int index_cross = 0; index_cross < floor(1 / param.cv_proportion);
         ++index_cross) {
      // std::cout << "\tcross = " << index_cross << "\n";
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
      xt::xarray<double> x0 =
          param.EPS * xt::ones<double>({model.C, model.R, model.T});
      model.nb_observations = nb_observations_current;
      model.nb_arrivals = nb_calls_current;
      xt::xarray<double> x =
          projected_gradient_armijo_feasible<RegularizedModel>(model, param,
                                                               x0);
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
      // std::cout << "\t\t f = " << f << "\n";
      likelihood += f;
    }
    likelihood = likelihood / floor(1 / param.cv_proportion);
    // printf("indexAlpha = %d, w = %f, likelihood = %f\n", index_alpha,
    //        alphas[index_alpha], likelihood);
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
  param.max_iter = 30;
  xt::xarray<double> x =
      param.EPS * xt::ones<double>({model.C, model.R, model.T});
  auto f_val =
      projected_gradient_armijo_feasible<RegularizedModel>(model, param, x);
  auto dt = std::chrono::high_resolution_clock::now();
  wall_time = std::chrono::duration_cast<std::chrono::seconds>(dt - t0).count();
  model.nb_observations = initial_nb_obs;
  model.nb_arrivals = initial_nb_arrivals;
  // printf("best_weight = %f\n", best_weight);
  return {wall_time, best_weight, x};
}

// CrossValidationResult cross_validation2(Param &param, RegularizedModel
// &model,
//                                         xt::xarray<int> &sample,
//                                         std::vector<double> &group_weights) {
//   if (sample.dimension() != 4 || sample.shape(0) != model.T ||
//       sample.shape(1) != model.R || sample.shape(2) != model.C) {
//     std::cout
//         << "Incorrect sample shape. Sample must be shaped {T, R, C, K}.\n ";
//     exit(1);
//   }
//   std::vector<double> alphas = group_weights;
//   auto t0 = std::chrono::high_resolution_clock::now();
//   int nb_observations_total = sample.shape(3);
//   int nb_groups = model.groups.size();
//   double min_loss = GRB_INFINITY;
//   double wall_time = 0;
//   int nb_in_block = floor(nb_observations_total * param.cv_proportion);
//   xt::xarray<int> initial_nb_obs = model.nb_observations;
//   xt::xarray<int> initial_nb_arrivals = model.nb_arrivals;
//   // fmt::print("cross validation weights.size = {}\n",
//   group_weights.size()); double best_alpha = GRB_INFINITY; double best_weight
//   = GRB_INFINITY;
//   // fmt::print("Running cross validation with proportion = {} and weights =
//   // {}\n", proportion, group_weights);

//   for (int index_alpha = 0; index_alpha < alphas.size(); ++index_alpha) {
//     double likelihood = 0;
//     model.alpha = alphas[index_alpha] * xt::ones<double>({model.R, model.R});
//     model.weights =
//         std::vector<double>(model.groups.size(), group_weights[index_alpha]);

//     for (size_t index_obs = 0; index_obs < nb_observations_total;
//     ++index_obs) {
//       for (size_t index_cross = 0; index_cross < nb_observations_total;
//            ++index_cross) {
//         for (int c = 0; c < model.C; ++c) {
//           for (int r = 0; r < model.R; ++r) {
//             for (int t = 0; t < model.T; ++t) {
//               if (index_cross != nb_obs)
//                 ++nb_observations_current(c, r, t);
//               nb_calls_current(c, r, t) += sample(t, r, c, index_cross);
//             }
//           }
//         }
//       }
//       xt::xarray<double> x =
//           param.EPS * xt::ones<double>({model.C, model.R, model.T});
//       model.nb_observations = nb_observations_current;
//       model.nb_arrivals = nb_calls_current;
//       xt::xarray<double> result_x =
//           projected_gradient_armijo_feasible<RegularizedModel>(model, param,
//           x);
//       xt::xarray<int> nb_calls_remaining =
//           xt::zeros<int>({model.C, model.R, model.T});
//       for (int c = 0; c < model.C; ++c) {
//         for (int r = 0; r < model.R; ++r) {
//           for (int t = 0; t < model.T; ++t) {
//             nb_calls_remaining(c, r, t) += sample(t, r, c, index_obs);
//           }
//         }
//       }

//       double f = 0;
//       for (int c = 0; c < model.C; ++c) {
//         for (int r = 0; r < model.R; ++r) {
//           for (int t = 0; t < model.T; ++t) {
//             double current_lambda = result_x(c, r, t);
//             f += current_lambda * model.durations[t] -
//                  nb_calls_remaining(c, r, t) * log(current_lambda);
//           }
//         }
//       }
//     }
//   }
// }

AppParameters load_options(int argc, char *argv[], po::variables_map &vm) {
  std::string config_file;
  // Declare a group of options that will be
  // allowed only on command line
  po::options_description generic("Generic Options");
  generic.add_options()("help,h", "Display this help message.")(
      "file,f", po::value<std::string>()->default_value(""),
      "Path to configuration file.");

  // Declare a group of options that will be
  // allowed both on command line and in
  // config file
  po::options_description config("Configuration");
  config.add_options()(
      "EPS,E", po::value<double>()->default_value(1e-5),
      "Epsilon for feasibility and convergence checks. Default = 1e-5")(
      "sigma,s", po::value<double>()->default_value(0.5),
      "Sigma parameter of armijo step. Default = 0.5")(
      "gap,g", po::value<double>()->default_value(0.01),
      "Gap used in stopping criterion. Default = 0.01")(
      "max_iter,I", po::value<int>()->default_value(30),
      "Max number of iterations used in stopping criterion. Default = 30")(
      "lower_lambda, L", po::value<double>()->default_value(1e-6),
      "lower bound on decision variables for both models. Default = 1e-6")(
      "upper_lambda, U", po::value<double>()->default_value(1.0),
      "upper bound on decision variables for both models. Default = 1.0")(
      "beta_bar, B", po::value<double>()->default_value(2.0),
      "Initial step size for projected gradient. Default = 2.0")(
      "cv_proportion,C", po::value<double>()->default_value(0.2),
      "Proportion of samples used as training in cross validation. "
      "Default = 0.2             ")(
      "model_type,M", po::value<std::string>()->default_value("no_reg"),
      "Chooses between RegularizedModel (option no_reg) or "
      "CovariatesModel(option reg).")(
      "method,m", po::value<std::string>()->default_value("calibration"),
      "Specifies the method being used. Default = calibration. Options are "
      "calibration, and cross_validation. If calibration and no_reg is set, "
      "parameter "
      "alpha_regions_file and time_groups_file must also be set. If "
      "cross_validation is "
      "set, cv_proportion and "
      "cv_weights_file must also be set.")(
      "algorithm,A", po::value<std::string>()->default_value("feasible"),
      "Specifies the projected gradient algorithm. Default = feasible. Options "
      "are feasible (projected gradient along the feasible region) and "
      "boundary (projected gradient along the boundary).")(
      "info_file,i", po::value<std::string>()->default_value(""),
      "Path to file with general information about the model. Default = "
      "''                         ")(
      "arrivals_file,a", po::value<std::string>()->default_value(""),
      "Path to file with arrivals data. Default = ''")(
      "neighbors_file,n", po::value<std::string>()->default_value(""),
      "Path to file with neighbors data. Default = ''")(
      "alpha_regions_file", po::value<std::string>()->default_value(""),
      "Path to file containing weight matrix for space regularization. Default "
      "= ''                ")(
      "time_groups_file", po::value<std::string>()->default_value("groups.txt"),
      "Path to file containing time groups information.")(
      "durations_file",
      po::value<std::string>()->default_value("durations.txt"),
      "Path to durations file. default = durations.txt")(
      "duration,d", po::value<double>()->default_value(1.0),
      "Duration of each period")(
      "cv_weights_file,W", po::value<std::string>()->default_value(""),
      "Path to file with weights that must be tested. Weights must be provided "
      "in one line, separated by spaces. Default = ''")(
      "output_file,O", po::value<std::string>()->default_value("output.txt"),
      "Path where to write output. Default = output.txt");

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
  AppParameters app_params;
  app_params.EPS = vm["EPS"].as<double>();
  app_params.sigma = vm["sigma"].as<double>();
  app_params.gap = vm["gap"].as<double>();
  app_params.max_iter = vm["max_iter"].as<int>();
  app_params.lower_lambda = vm["lower_lambda"].as<double>();
  app_params.upper_lambda = vm["upper_lambda"].as<double>();
  app_params.beta_bar = vm["beta_bar"].as<double>();
  app_params.cv_proportion = vm["cv_proportion"].as<double>();
  app_params.model_type = vm["model_type"].as<std::string>();
  app_params.method = vm["method"].as<std::string>();
  app_params.algorithm = vm["algorithm"].as<std::string>();
  app_params.info_file = vm["info_file"].as<std::string>();
  app_params.arrivals_file = vm["arrivals_file"].as<std::string>();
  app_params.neighbors_file = vm["neighbors_file"].as<std::string>();
  app_params.alpha_regions_file = vm["alpha_regions_file"].as<std::string>();
  app_params.time_groups_file = vm["time_groups_file"].as<std::string>();
  app_params.durations_file = vm["durations_file"].as<std::string>();
  app_params.cv_weights_file = vm["cv_weights_file"].as<std::string>();
  app_params.output_file = vm["output_file"].as<std::string>();
  return app_params;
}

}  // namespace laspated

#endif
