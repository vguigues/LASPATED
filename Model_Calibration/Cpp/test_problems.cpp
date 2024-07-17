#include "test_problems.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#define RED 31
#define GREEN 32

using laspated::projected_gradient_armijo_feasible;

bool verify(const xt::xarray<double>& x, const xt::xarray<double>& expected,
            ulong N, double EPS) {
  bool optimal = true;
  for (ulong i = 0; i < N; ++i) {
    if (fabs(x(i) - expected(i)) > EPS) {
      optimal = false;
      break;
    }
  }
  return optimal;
}

void diff(const xt::xarray<double>& x, const xt::xarray<double>& expected,
          ulong N) {
  using std::cout;
  for (ulong i = 0; i < N; ++i) {
    printf("\t%ld: exp %f, got %f\n", i, expected(i), x(i));
  }
}

int test_problem_unc3(Param& param) {
  ModelUncProblem3 model(param);
  ulong N = 2;
  xt::xarray<double> expected({1, 16});
  xt::xarray<double> x0 = xt::zeros<double>({N});
  xt::xarray<double> x =
      projected_gradient_armijo_feasible<ModelUncProblem3>(model, param, x0);
  bool optimal = fabs(model.f(x) - model.f(expected)) < param.EPS ||
                 verify(x, expected, N, param.EPS);
  if (!optimal) {
    printf("\033[%dmTEST Unc3 FAILED.\033[m Expected/Result values:\n", RED);
    printf("f(expected) = %f, f(x) = %f\n", model.f(expected), model.f(x));
    diff(x, expected, N);
    return 1;
  } else {
    printf("\033[%dmTest problem Unc3 passed\033[m\n", GREEN);
  }
  return 0;
}
#if USE_GUROBI == 1
GRBEnv env;
int test_problem3(Param& param) {
  ModelProblem3 model(env, param);
  ulong N = 5;
  xt::xarray<double> expected({-0.76744, 0.25581, 0.62791, -0.11628, 0.25581});
  xt::xarray<double> x0 = 2 * xt::zeros<double>({N});
  xt::xarray<double> x =
      projected_gradient_armijo_feasible<ModelProblem3>(model, param, x0);
  bool optimal = fabs(model.f(x) - model.f(expected)) < param.EPS ||
                 verify(x, expected, N, param.EPS);

  if (!optimal) {
    printf("\033[%dmTEST 3 FAILED.\033[m Expected/Result values:\n", RED);
    printf("f(expected) = %f, f(x) = %f\n", model.f(expected), model.f(x));
    diff(x, expected, N);
    return 1;
  } else {
    printf("\033[%dmTest problem 3 passed\033[m\n", GREEN);
  }
  return 0;
}

int test_problem14(Param& param) {
  ModelProblem14 model(env, param);
  ulong N = 2;
  xt::xarray<double> expected({0.5 * (sqrt(7) - 1), 0.25 * (sqrt(7) + 1)});
  xt::xarray<double> x0({2, 2});
  xt::xarray<double> x =
      projected_gradient_armijo_feasible<ModelProblem14>(model, param, x0);
  bool optimal = fabs(model.f(x) - model.f(expected)) < param.EPS ||
                 verify(x, expected, N, param.EPS);
  if (!optimal) {
    printf("\033[%dmTEST 14 FAILED.\033[m Expected/Result values:\n", RED);
    printf("f(expected) = %f, f(x) = %f\n", model.f(expected), model.f(x));
    diff(x, expected, N);
    return 1;
  } else {
    printf("\033[%dmTest problem 14 passed\033[m\n", GREEN);
  }
  return 0;
}

int test_problem28(Param& param) {
  ModelProblem28 model(env, param);
  ulong N = 3;
  xt::xarray<double> expected({0.5, -0.5, 0.5});
  xt::xarray<double> x0({-4, 1, 1});
  xt::xarray<double> x =
      projected_gradient_armijo_feasible<ModelProblem28>(model, param, x0);
  bool optimal = fabs(model.f(x) - model.f(expected)) < param.EPS ||
                 verify(x, expected, N, param.EPS);

  if (!optimal) {
    printf("\033[%dmTEST 28 FAILED.\033[m Expected/Result values:\n", RED);
    printf("f(expected) = %f, f(x) = %f\n", model.f(expected), model.f(x));
    diff(x, expected, N);
    return 1;
  } else {
    printf("\033[%dmTest problem 28 passed\033[m\n", GREEN);
  }
  return 0;
}

int test_problem30(Param& param) {
  ModelProblem30 model(env, param);
  ulong N = 3;
  xt::xarray<double> expected({1, 0, 0});
  xt::xarray<double> x0({1, 1, 1});
  xt::xarray<double> x =
      projected_gradient_armijo_feasible<ModelProblem30>(model, param, x0);
  bool optimal = fabs(model.f(x) - model.f(expected)) < param.EPS ||
                 verify(x, expected, N, param.EPS);

  if (!optimal) {
    printf("\033[%dmTEST 30 FAILED.\033[m Expected/Result values:\n", RED);
    printf("f(expected) = %f, f(x) = %f\n", model.f(expected), model.f(x));
    diff(x, expected, N);
    return 1;
  } else {
    printf("\033[%dmTest problem 30 passed\033[m\n", GREEN);
  }

  return 0;
}

#endif

int main(int argc, char const* argv[]) {
  Param param;

  param.lower_lambda = -1000.0;
  param.upper_lambda = 1000.0;
  param.max_iter = 40;
  param.EPS = 0.001;

  test_problem_unc3(param);
#if USE_GUROBI == 1
  test_problem3(param);
  test_problem14(param);
  test_problem28(param);
  test_problem30(param);
#endif
}
