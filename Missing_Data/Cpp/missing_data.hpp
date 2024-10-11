#include "laspated.h"

class MissingLambdaModel {
 public:
  std::string name = "Regularized";
  ulong C;
  ulong D;
  ulong T;
  ulong R;

  xt::xarray<int> nb_observations;
  xt::xarray<int> nb_arrivals;
  xt::xarray<int> nb_missing_arrivals;
  std::vector<std::vector<std::pair<int, int>>> groups;
  xt::xarray<int> which_group;
  xt::xarray<double> durations;
  xt::xarray<double> alpha;
  std::vector<std::vector<int>> neighbors;
  std::vector<double> weight;
  laspated::Param& param;

  MissingLambdaModel(xt::xarray<int>& a_nb_observations,
                     xt::xarray<int>& a_nb_arrivals,
                     xt::xarray<int>& a_nb_missing_arrivals,
                     xt::xarray<double>& alphas, std::vector<double>& weights,
                     std::vector<std::vector<std::pair<int, int>>>& a_groups,
                     std::vector<std::vector<int>>& a_neighbors,
                     xt::xarray<double>& a_durations, laspated::Param& a_param)
      : param(a_param) {
    nb_observations = a_nb_observations;
    C = nb_observations.shape(0);
    D = nb_observations.shape(1);
    T = nb_observations.shape(2);
    R = nb_observations.shape(3);
    nb_arrivals = a_nb_arrivals;
    nb_missing_arrivals = a_nb_missing_arrivals;
    groups = a_groups;
    xt::xarray<int> which_group = -1 * xt::ones<int>({D, T});
    for (int g = 0; g < groups.size(); ++g) {
      for (auto elem : groups[g]) {
        int d = elem.first;
        int t = elem.second;
        which_group(d, t) = g;
      }
    }
    durations = a_durations;
    alpha = alphas;
    neighbors = a_neighbors;
    weight = weights;
  }

  double f(xt::xarray<double>& x) {
    double obj = 0.0;

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            obj +=
                nb_observations(c, d, t, r) * x(c, d, t, r) * durations(d, t) -
                nb_arrivals(c, d, t, r) * log(x(c, d, t, r));
            for (auto s : neighbors[r]) {
              obj += 0.5 * alpha(r, s) * nb_observations(c, d, t, r) *
                     nb_observations(c, d, t, s) *
                     pow(x(c, d, t, r) - x(c, d, t, s), 2);
            }
          }
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          double sum = 0.0;
          for (int r = 0; r < R; ++r) {
            sum += x(c, d, t, r);
          }
          obj -= (nb_missing_arrivals(c, d, t) * log(sum));
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < R; ++r) {
        for (int grindex = 0; grindex < groups.size(); ++grindex) {
          auto& group = groups[grindex];
          for (auto& elem : group) {
            int d = elem.first;
            int t = elem.second;
            for (auto& elem1 : group) {
              int d1 = elem1.first;
              int t1 = elem1.second;
              obj += 0.5 * weight[grindex] * nb_observations(c, d, t, r) *
                     nb_observations(c, d1, t1, r) *
                     pow(x(c, d, t, r) - x(c, d1, t1, r), 2);
            }
          }
        }
      }
    }

    return obj;
  }

  xt::xarray<double> gradient(xt::xarray<double>& x) {
    xt::xarray<double> grad = xt::zeros<double>(x.shape());
    xt::xarray<double> sum = xt::zeros<double>({C, D, T});
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          double aux = 0.0;
          for (int r = 0; r < R; ++r) {
            aux += x(c, d, t, r);
          }
          sum(c, d, t) = aux;
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          double grad_component =
              nb_observations(c, d, t, 0) * durations(d, t) -
              nb_missing_arrivals(c, d, t) / sum(c, d, t);
          for (int r = 0; r < R; ++r) {
            double grad_component1 =
                grad_component - (nb_arrivals(c, d, t, r) / x(c, d, t, r));
            for (int s : neighbors[r]) {
              grad_component1 += 2 * alpha(r, s) * nb_observations(c, d, t, r) *
                                 nb_observations(c, d, t, s) *
                                 (x(c, d, t, r) - x(c, d, t, s));
            }
            grad(c, d, t, r) = grad_component1;
          }
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            for (auto& elem : groups[which_group(d, t)]) {
              int d1 = elem.first;
              int t1 = elem.second;
              grad(c, d, t, r) += 2 * weight[which_group(d, t)] *
                                  nb_observations(c, d, t, r) *
                                  nb_observations(c, d1, t1, r) *
                                  (x(c, d, t, r) - x(c, d1, t1, r));
            }
          }
        }
      }
    }

    return grad;
  }

  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    double rhs = 0.0;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            rhs += grad(c, d, t, r) * dir(c, d, t, r);
          }
        }
      }
    }

    return rhs;
  }
  xt::xarray<double> projection(xt::xarray<double>& x) {
    xt::xarray<double> z = xt::zeros<double>(x.shape());
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            z(c, d, t, r) = std::max(param.EPS, x(c, d, t, r));
          }
        }
      }
    }
    return z;
  }
  bool is_feasible(xt::xarray<double>& x) { return true; }
  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    return 0.0;
  }
};

int factorial(int n) {
  int prod = 1;
  for (int i = 1; i <= n; ++i) {
    prod *= i;
  }

  return prod;
}
class MissingModel2 {
 public:
  std::string name = "Regularized";
  ulong C;
  ulong D;
  ulong T;
  ulong R;
  ulong S;

  xt::xarray<int> nb_observations;
  xt::xarray<int> nb_arrivals;
  xt::xarray<int> nb_missing_arrivals;
  xt::xarray<int> sample_arrivals;
  xt::xarray<int> sample_missing_arrivals;
  xt::xarray<int> mn_samples;
  xt::xarray<double> durations;
  laspated::Param& param;

  MissingModel2(xt::xarray<int> a_nb_observations,
                xt::xarray<int> a_nb_arrivals,
                xt::xarray<int> a_nb_missing_arrivals,
                xt::xarray<int> a_sample_arrivals,
                xt::xarray<int> a_sample_missing_arrivals,
                xt::xarray<int> a_mn_samples, xt::xarray<double> a_durations,
                laspated::Param& a_param)
      : param(a_param) {
    nb_observations = a_nb_observations;
    nb_arrivals = a_nb_arrivals;
    nb_missing_arrivals = a_nb_missing_arrivals;
    sample_arrivals = a_sample_arrivals;
    sample_missing_arrivals = a_sample_missing_arrivals;
    mn_samples = a_mn_samples;
    durations = a_durations;

    C = nb_observations.shape(0);
    D = nb_observations.shape(1);
    T = nb_observations.shape(2);
    R = nb_observations.shape(3);
    S = mn_samples.shape(4);
  }

  double f(xt::xarray<double>& x) {
    double obj = 0.0;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            obj +=
                nb_observations(c, d, t, r) * x(c, d, t, r) * durations(d, t);
          }
        }
      }
    }

    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int n = 0; n < nb_observations(c, d, t, 0); ++n) {
            double u = 0;
            for (int s = 0; s < S; ++s) {
              double aux = 1;
              for (int r = 0; r < R; ++r) {
                int power = mn_samples(c, d, t, n, s, r) +
                            sample_arrivals(c, d, t, r, n);
                double iaux = exp(-durations(d, t) * x(c, d, t, r));
                iaux = iaux * pow(durations(d, t) * x(c, d, t, r), power);
                double denom = factorial(power);
                iaux = iaux / denom;
                aux = aux * iaux;
              }
              u += aux;
            }
            obj -= log(u / S);
          }
        }
      }
    }
    return obj;
  }

  xt::xarray<double> gradient(xt::xarray<double>& x) {
    xt::xarray<double> grad = xt::zeros<double>(x.shape());
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            grad(c, d, t, r) = durations(d, t) * nb_observations(c, d, t, r);
            for (int n = 0; n < nb_observations(c, d, t, 0); ++n) {
              double u = 0;
              for (int s = 0; s < S; ++s) {
                double aux = 1;
                for (int k = 0; k < R; ++k) {
                  int power = mn_samples(c, d, t, n, s, k) +
                              sample_arrivals(c, d, t, k, n);
                  double iaux = exp(-durations(d, t) * x(c, d, t, k));
                  iaux = iaux * pow(durations(d, t) * x(c, d, t, k), power);
                  double denom = factorial(power);
                  iaux = iaux / denom;
                  aux = aux * iaux;
                }
                u += aux;
              }
              u = u / S;
              double up = 0;
              for (int s = 0; s < S; ++s) {
                int power_r = mn_samples(c, d, t, n, s, r) +
                              sample_arrivals(c, d, t, r, n);
                double faux = exp(-x(c, d, t, r) * durations(d, t)) *
                              pow(x(c, d, t, r), power_r) *
                              (-durations(d, t) * x(c, d, t, r) + power_r);
                double iaux1 = 1;
                for (int k = 0; k < R; ++r) {
                  int power_k = mn_samples(c, d, t, n, s, k) +
                                sample_arrivals(c, d, t, k, n);
                  iaux1 = iaux1 * pow(durations(d, t), power_k);
                  double denom = factorial(power_k);
                  iaux1 = iaux1 / denom;
                }
                double iaux2 = 1;
                for (int k = 0; k < R; ++r) {
                  if (k != r) {
                    int power_k = mn_samples(c, d, t, n, s, k) +
                                  sample_arrivals(c, d, t, k, n);
                    iaux2 = iaux2 * exp(-durations(d, t) * x(c, d, t, k)) *
                            pow(x(c, d, t, k), power_k);
                  }
                }
                up = up + faux * iaux1 * iaux2;
              }
              up = up / S;
              grad(c, d, t, r) -= up / u;
            }
          }
        }
      }
    }
    return grad;
  }

  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    double rhs = 0.0;
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            rhs += grad(c, d, t, r) * dir(c, d, t, r);
          }
        }
      }
    }

    return rhs;
  }
  xt::xarray<double> projection(xt::xarray<double>& x) {
    xt::xarray<double> z = xt::zeros<double>(x.shape());
    for (int c = 0; c < C; ++c) {
      for (int d = 0; d < D; ++d) {
        for (int t = 0; t < T; ++t) {
          for (int r = 0; r < R; ++r) {
            z(c, d, t, r) = std::max(param.EPS, x(c, d, t, r));
          }
        }
      }
    }
    return z;
  }
  bool is_feasible(xt::xarray<double>& x) { return true; }
  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    return 0.0;
  }
};
