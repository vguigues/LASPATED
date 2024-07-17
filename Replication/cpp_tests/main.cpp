#include <sstream>
#include <string>
#include <vector>

#include "tests.h"

void ex_no_reg() {
  using namespace std;
  vector<vector<double>> initial_weights{
      {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 50, 100},
      {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 50, 100}};

  vector<int> test_nb_groups{2, 4};
  vector<int> test_nb_weeks{1, 10, 50, 500};
  vector<int> test_neighbor_factors{0, 1};
  vector<int> test_constant_lambdas{0, 1};

  // vector<int> test_nb_groups{2};
  // vector<int> test_nb_weeks{10};
  // vector<int> test_neighbor_factors{1};
  // vector<int> test_constant_lambdas{1};

  for (auto constant_lambdas : test_constant_lambdas) {
    const size_t cl = constant_lambdas;
    for (auto nb_weeks : test_nb_weeks) {
      // nb_weeks *= 52; // Consistent with test 2
      vector<double> test_weights;
      for (size_t i = 0; i < initial_weights[cl].size(); ++i) {
        test_weights.push_back(initial_weights[cl][i] / nb_weeks);
      }
      for (auto nb_groups : test_nb_groups) {
        for (auto neighbor_factor : test_neighbor_factors) {
          cout << "Running test1: ";
          cout << "G = " << nb_groups << ", ";
          cout << "weeks = " << nb_weeks << ", ";
          cout << "neighbors = " << neighbor_factor << ", ";
          cout << "model_type = " << constant_lambdas << ", ";
          bool do_cross_validation =
              nb_groups == 2 && neighbor_factor == 1 && nb_weeks > 1;
          cout << "Cross validation = " << do_cross_validation << "\n";
          auto result =
              test1(nb_weeks, nb_groups, neighbor_factor, constant_lambdas,
                    test_weights, do_cross_validation);

          // Register empirical error
          // Plot relative error by penalties
          // Register min_error and corresponding weight
          stringstream err_arq_name;
          int ex_number = (constant_lambdas == 1) ? 1 : 2;
          err_arq_name << "cpp_tests/results/ex" << ex_number
                       << "/err_by_weight_w" << nb_weeks << "_g" << nb_groups
                       << "_a" << neighbor_factor << "_m" << constant_lambdas
                       << ".txt";
          ofstream err_arq(err_arq_name.str(), ios::out);
          err_arq << result.min_error << " " << result.min_w << " "
                  << test_weights.size() << " " << result.mean_emp << "\n";
          for (size_t i = 0; i < test_weights.size(); ++i) {
            err_arq << result.error_by_weight[i] << "\n";
          }
          if (do_cross_validation) {
            err_arq << result.mean_cv << " " << result.w_cv << "\n";
          }
          err_arq.close();
          if (neighbor_factor == 1 && nb_groups == 2 &&
              (nb_weeks == 1 || nb_weeks == 10) && constant_lambdas == 1) {
            stringstream lambda_by_t;
            lambda_by_t << "cpp_tests/results/ex" << ex_number
                        << "/lambdas_r1_w" << nb_weeks << ".txt";
            ofstream lambda_arq1(lambda_by_t.str(), ios::out);
            lambda_by_t.str("");
            lambda_by_t << "cpp_tests/results/ex" << ex_number
                        << "/lambdas_r6_w" << nb_weeks << ".txt";
            ofstream lambda_arq6(lambda_by_t.str(), ios::out);
            for (ulong t = 0; t < result.T; ++t) {
              double sum_theo1 = 0.0;
              double sum_emp1 = 0.0;
              double sum_est1 = 0.0;
              double sum_theo6 = 0.0;
              double sum_emp6 = 0.0;
              double sum_est6 = 0.0;
              for (ulong c = 0; c < result.C; ++c) {
                sum_theo1 += result.theoretical_lambda(c, 0, t);
                sum_emp1 += result.empirical_lambda(c, 0, t);
                sum_est1 += result.estimated_lambda(c, 0, t);
                sum_theo6 += result.theoretical_lambda(c, 5, t);
                sum_emp6 += result.empirical_lambda(c, 5, t);
                sum_est6 += result.estimated_lambda(c, 5, t);
              }
              lambda_arq1 << t << " " << sum_theo1 << " " << sum_emp1 << " "
                          << sum_est1 << "\n";
              lambda_arq6 << t << " " << sum_theo6 << " " << sum_emp6 << " "
                          << sum_est6 << "\n";
            }
            lambda_arq1.close();
            lambda_arq6.close();
          }
        }
      }
    }
  }
}

#if USE_GUROBI == 1
void ex2_reg() {
  using namespace std;
  // vector<int> test_nb_weeks{10, 50, 500};
  // vector<int> test_nb_weeks{50, 500};
  vector<int> test_use_holidays{0, 1};
  vector<int> test_nb_years{1, 10, 15};
  // vector<int> test_use_holidays{1};

  ofstream result_file("cpp_tests/results/ex3/table_ex3.txt", ios::out);
  result_file
      << "nb_weeks\tholidays\tCovariates\tReg 1\tReg 2\tReg 3\tEmpirical\n";
  cout << "nb_weeks\tholidays\tCovariates\tReg 1\tReg 2\tReg 3\tEmpirical\n";
  for (auto use_holidays : test_use_holidays) {
    for (auto nb_years : test_nb_years) {
      cout << "Running test2: ";
      cout << "weeks = " << nb_years << ", ";
      cout << "use_holidays = " << use_holidays << "\n";
      Result2 result = test2(nb_years, use_holidays);
      result_file << nb_years * 52 << " & " << use_holidays << " & "
                  << result.err_cov << " & " << result.min_err1 << " & "
                  << result.min_err2 << " & " << result.min_err3 << " & "
                  << result.err_emp << "\n";
      cout << nb_years << " " << use_holidays << " " << result.err_cov << " "
           << result.min_err1 << " " << result.min_err2 << " "
           << result.min_err3 << " " << result.err_emp << "\n";
      // cin.get();
    }
  }
  result_file.close();
}

void real_data(const std::string& data_dir) {
  using namespace std;
  const std::string base_dir = data_dir + "/discretizations";
  vector<string> test_base_paths{base_dir + "/rect", base_dir + "/hex",
                                 base_dir + "/district"};

  for (auto& base_path : test_base_paths) {
    auto result = test3(base_path);
    ulong C = result.C;
    ulong R = result.R;
    ulong T = result.T;
    stringstream t_rates_filename;
    t_rates_filename << "cpp_tests/results/real_data/rates_by_t_r" << R
                     << ".txt";
    ofstream rates_by_t(t_rates_filename.str(), ios::out);
    for (int t = 0; t < T; ++t) {
      vector<double> sum_c_emp(C, 0);
      vector<double> sum_c_reg(C, 0);
      vector<double> sum_c_cov(C, 0);
      rates_by_t << t << " ";
      for (int c = 0; c < C; ++c) {
        for (int r = 0; r < R; ++r) {
          sum_c_emp[c] += result.empirical_rates(c, r, t);
          sum_c_reg[c] += result.regularized_rates(c, r, t);
          sum_c_cov[c] += result.covariates_rates(c, r, t);
        }
        rates_by_t << sum_c_emp[c] << " " << sum_c_reg[c] << " " << sum_c_cov[c]
                   << " ";
      }
      rates_by_t << "\n";
    }
    rates_by_t.close();
    stringstream r_rates_filename;
    r_rates_filename << "cpp_tests/results/real_data/rates_by_region_r" << R
                     << ".txt";
    ofstream rates_by_r(r_rates_filename.str(), ios::out);
    for (int r = 0; r < R; ++r) {
      vector<double> sum_c_emp(C, 0);
      vector<double> sum_c_reg(C, 0);
      vector<double> sum_c_cov(C, 0);
      rates_by_r << r << " ";
      for (int c = 0; c < C; ++c) {
        for (int t = 0; t < T; ++t) {
          sum_c_emp[c] += result.empirical_rates(c, r, t);
          sum_c_reg[c] += result.regularized_rates(c, r, t);
          sum_c_cov[c] += result.covariates_rates(c, r, t);
        }
        rates_by_r << sum_c_emp[c] << " " << sum_c_reg[c] << " "
                   << sum_c_cov[c] / 0.5 << " ";
      }
      rates_by_r << "\n";
    }
    rates_by_r.close();
  }
}
#endif
void test_det() {
  using namespace std;
  vector<vector<double>> initial_weights{
      {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 50, 100},
      {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 50, 100}};

  vector<int> test_nb_groups{2};
  // vector<int> test_nb_weeks{1, 10, 50, 500};
  vector<int> test_nb_weeks{10};
  vector<string> test_samples{"sample1_10.txt"};
  vector<int> test_neighbor_factors{1};
  // vector<int> test_constant_lambdas{0, 1};
  vector<int> test_constant_lambdas{1};

  for (auto constant_lambdas : test_constant_lambdas) {
    const size_t cl = constant_lambdas;
    int index_week = 0;
    for (auto nb_weeks : test_nb_weeks) {
      // nb_weeks *= 52; // Consistency with test 2
      vector<double> test_weights;
      for (size_t i = 0; i < initial_weights[cl].size(); ++i) {
        test_weights.push_back(initial_weights[cl][i] / nb_weeks);
      }
      for (auto nb_groups : test_nb_groups) {
        for (auto neighbor_factor : test_neighbor_factors) {
          cout << "Running test1: ";
          cout << "G = " << nb_groups << ", ";
          cout << "weeks = " << nb_weeks << ", ";
          cout << "neighbors = " << neighbor_factor << ", ";
          cout << "model_type = " << constant_lambdas << ", ";
          bool do_cross_validation =
              nb_groups == 2 && neighbor_factor == 1 && nb_weeks > 1;
          cout << "Cross validation = " << do_cross_validation << "\n";
          auto result = test1_deterministic(
              test_samples[index_week], nb_weeks, nb_groups, neighbor_factor,
              constant_lambdas, test_weights, do_cross_validation);

          // Register empirical error
          // Plot relative error by penalties
          // Register min_error and corresponding weight
          stringstream err_arq_name;
          int ex_number = (constant_lambdas == 1) ? 1 : 2;
          err_arq_name << "cpp_tests/results/ex" << ex_number
                       << "/err_by_weight_w" << nb_weeks << "_g" << nb_groups
                       << "_a" << neighbor_factor << "_m" << constant_lambdas
                       << ".txt";
          ofstream err_arq(err_arq_name.str(), ios::out);
          err_arq << result.min_error << " " << result.min_w << " "
                  << test_weights.size() << " " << result.mean_emp << "\n";
          for (size_t i = 0; i < test_weights.size(); ++i) {
            err_arq << result.error_by_weight[i] << "\n";
          }
          if (do_cross_validation) {
            err_arq << result.mean_cv << " " << result.w_cv << "\n";
          }
          err_arq.close();
          if (neighbor_factor == 1 && nb_groups == 2 &&
              (nb_weeks == 1 || nb_weeks == 10) && constant_lambdas == 1) {
            stringstream lambda_by_t;
            lambda_by_t << "cpp_tests/results/ex" << ex_number
                        << "/lambdas_r1_w" << nb_weeks << ".txt";
            ofstream lambda_arq1(lambda_by_t.str(), ios::out);
            lambda_by_t.str("");
            lambda_by_t << "cpp_tests/results/ex" << ex_number
                        << "/lambdas_r6_w" << nb_weeks << ".txt";
            ofstream lambda_arq6(lambda_by_t.str(), ios::out);
            for (ulong t = 0; t < result.T; ++t) {
              double sum_theo1 = 0.0;
              double sum_emp1 = 0.0;
              double sum_est1 = 0.0;
              double sum_theo6 = 0.0;
              double sum_emp6 = 0.0;
              double sum_est6 = 0.0;
              for (ulong c = 0; c < result.C; ++c) {
                sum_theo1 += result.theoretical_lambda(c, 0, t);
                sum_emp1 += result.empirical_lambda(c, 0, t);
                sum_est1 += result.estimated_lambda(c, 0, t);
                sum_theo6 += result.theoretical_lambda(c, 5, t);
                sum_emp6 += result.empirical_lambda(c, 5, t);
                sum_est6 += result.estimated_lambda(c, 5, t);
              }
              lambda_arq1 << t << " " << sum_theo1 << " " << sum_emp1 << " "
                          << sum_est1 << "\n";
              lambda_arq6 << t << " " << sum_theo6 << " " << sum_emp6 << " "
                          << sum_est6 << "\n";
            }
            lambda_arq1.close();
            lambda_arq6.close();
          }
        }
      }
      ++index_week;
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc < 3 || argc > 5) {
    printf(
        "Error: wrong usage\nUsage: laspated -e "
        "[ex_no_reg|ex2_reg|real_data|all] (--data_dir [data_path])\nThe "
        "data_dir must be provided in case of testing real data\n");
    exit(1);
  }
  std::string example = argv[2];
  if ((example == "real_data" || example == "all") && argc != 5) {
    printf(
        "Error: wrong usage with real data.\nUsage: laspated -e "
        "[ex_no_reg|ex2_reg|real_data|all] --data_dir [data_path]\n");
    exit(1);
  }
  std::string data_path = argv[4];
  if (example == "ex_no_reg") {
    ex_no_reg();
  } else if (example == "ex2_reg") {
#if USE_GUROBI == 1
    ex2_reg();
#else
    std::cout << "Cannot run ex2_reg because gurobi is not being used\n";
#endif
  } else if (example == "real_data") {
#if USE_GUROBI == 1
    real_data(data_path);
#else
    std::cout << "Cannot run real_data because gurobi is not being used\n";
#endif
  } else if (example == "all") {
    ex_no_reg();
#if USE_GUROBI == 1
    ex2_reg();
    real_data(data_path);
#else
    printf(
        "Warning: Skipped ex2_reg and real_data because gurobi is not being "
        "used.\n");
#endif
  } else {
    std::cout << "Unknown example " << example << "\n";
  }
  return 0;
}
