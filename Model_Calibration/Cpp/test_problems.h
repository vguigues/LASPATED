#ifndef M_PROBLEM3_H
#define M_PROBLEM3_H

#include <iostream>
#include <sstream>

#include "laspated.h"

using laspated::Param;
using std::cout;

double rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
  double rhs = 0.0;
  ulong n = grad.shape(0);
  for (ulong i = 0; i < n; ++i) {
    rhs += grad(i) * dir(i);
  }
  return rhs;
}

#if USE_GUROBI == 1
#include "gurobi_c++.h"
class ModelUncProblem3 {
 public:
  Param& params;
  std::string name = "Regularized";
  ModelUncProblem3(Param& a_param) : params(a_param) {}
  double f(xt::xarray<double>& x) {
    return -1 *
           ((70 - 2 * (x(0) + x(1))) * (x(0) + x(1)) - pow(x(0), 2) - 2 * x(1));
  }

  xt::xarray<double> gradient(xt::xarray<double>& x) {
    return xt::xarray<double>({-1 * (70 - 4 * (x(0) + x(1)) - 2 * x(0)),
                               -1 * (70 - 4 * (x(0) + x(1)) - 2)});
  }
  xt::xarray<double> projection(xt::xarray<double>& x) { return x; }

  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }

  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    double change = 0.0;
    ulong n = x.shape(0);
    for (ulong i = 0; i < n; ++i) {
      if (grad(i) < 0.0) {
        change += grad(i) * (params.upper_lambda - x(i));
      } else {
        change += grad(i) * (params.lower_lambda - x(i));
      }
    }

    return change;
  }

  bool is_feasible(xt::xarray<double>& x) { return x(0) > 0.0 && x(1) > 0.0; }
};

class ModelProblem3 {
 public:
  GRBEnv& env;
  Param& params;
  std::string name = "Covariates";
  ModelProblem3(GRBEnv& a_env, Param& a_param) : env(a_env), params(a_param) {}

  double f(xt::xarray<double>& x) {
    return pow(x(0) - x(1), 2) + pow(x(1) + x(2) - 2, 2) + pow(x(3) - 1, 2) +
           pow(x(4) - 1, 2);
  }

  xt::xarray<double> gradient(xt::xarray<double>& x) {
    xt::xarray<double> gradient = xt::zeros<double>(x.shape());
    gradient(0) = 2 * (x(0) - x(1));
    gradient(1) = -2 * (x(0) - x(1)) + 2 * (x(1) + x(2) - 2);
    gradient(2) = 2 * (x(1) + x(2) - 2);
    gradient(3) = 2 * (x(3) - 1);
    gradient(4) = 2 * (x(4) - 1);

    return gradient;
  }

  xt::xarray<double> projection(xt::xarray<double>& x) {
    GRBModel model(env);
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});
    std::stringstream name;
    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      y(i) = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                          name.str());
      name.str("");
    }

    GRBQuadExpr obj;

    for (ulong i = 0; i < n; ++i) {
      obj += 0.5 * y(i) * y(i) - x(i) * y(i);
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBLinExpr g1 = y(0) + 3 * y(1);
    model.addConstr(g1, GRB_EQUAL, 0.0, "g1");
    GRBLinExpr g2 = y(2) + y(3) - 2 * y(4);
    model.addConstr(g2, GRB_EQUAL, 0.0, "g2");
    GRBLinExpr g3 = y(1) - y(4);
    model.addConstr(g3, GRB_EQUAL, 0.0, "g3");

    model.update();
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_DualReductions, 0);
    model.optimize();

    auto status = model.get(GRB_IntAttr_Status);
    xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());

    if (status == GRB_OPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }

    return y_val;
  }

  bool is_feasible(xt::xarray<double>& x) {
    if (fabs(x(0) + 3 * x(1)) > params.EPS) {
      printf("Unfeasible at g1: %f", fabs(x(0) + 3 * x(1)));
      return false;
    }

    if (fabs(x(2) + x(3) - 2 * x(4)) > params.EPS) {
      printf("Unfeasible at g2: %f", fabs(x(2) + x(3) - 2 * x(4)));
      return false;
    }

    if (fabs(x(1) - x(4)) > params.EPS) {
      printf("Unfeasible at g3: %f", fabs(x(1) - x(4)));
      return false;
    }

    return true;
  }

  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }

  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    GRBModel model(env);
    std::stringstream name;
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});

    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      y(i) = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                          name.str());
      name.str("");
    }

    GRBLinExpr obj = 0;
    for (ulong i = 0; i < n; ++i) {
      obj += grad(i) * (y(i) - x(i));
    }

    model.setObjective(obj, GRB_MINIMIZE);
    GRBLinExpr g1 = y(0) + 3 * y(1);
    model.addConstr(g1, GRB_EQUAL, 0.0, "g1");
    GRBLinExpr g2 = y(2) + y(3) - 2 * y(4);
    model.addConstr(g2, GRB_EQUAL, 0.0, "g2");
    GRBLinExpr g3 = y(1) - y(4);
    model.addConstr(g3, GRB_EQUAL, 0.0, "g3");

    model.update();
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
};

class ModelProblem14 {
 public:
  GRBEnv& env;
  Param& params;
  std::string name = "Covariates";
  ModelProblem14(GRBEnv& a_env, Param& a_param) : env(a_env), params(a_param) {}

  double f(xt::xarray<double>& x) {
    return pow(x(0) - 2, 2) + pow(x(1) - 1, 2);
  }

  xt::xarray<double> gradient(xt::xarray<double>& x) {
    return xt::xarray<double>({2 * (x(0) - 2), 2 * (x(1) - 1)});
  }

  xt::xarray<double> projection(xt::xarray<double>& x) {
    GRBModel model(env);
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});
    std::stringstream name;
    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      y(i) = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                          name.str());
      name.str("");
    }

    GRBQuadExpr obj;

    for (ulong i = 0; i < n; ++i) {
      obj += 0.5 * y(i) * y(i) - x(i) * y(i);
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBQuadExpr g1 = -0.25 * y(0) * y(0) - y(1) * y(1) + 1;
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");
    GRBLinExpr g2 = y(0) - 2 * y(1) + 1;
    model.addConstr(g2, GRB_EQUAL, 0.0, "g2");

    model.update();
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_DualReductions, 0);
    model.optimize();

    auto status = model.get(GRB_IntAttr_Status);
    xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());

    if (status == GRB_OPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
    } else if (status == GRB_SUBOPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
      if (!is_feasible(y_val)) {
        cout << "Error. Suboptimal projection solution is not feasible!\n";
        model.write("projection_model.lp");
        exit(status);
      }
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }

    return y_val;
  }

  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }

  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    GRBModel model(env);
    std::stringstream name;
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});

    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      if (i == 0) {
        y(i) = model.addVar(1, 10, 0, GRB_CONTINUOUS, name.str());
      } else {
        y(i) = model.addVar(-10, 10, 0, GRB_CONTINUOUS, name.str());
      }
      name.str("");
    }

    GRBLinExpr obj = 0;
    for (ulong i = 0; i < n; ++i) {
      obj += grad(i) * (y(i) - x(i));
    }

    model.setObjective(obj, GRB_MINIMIZE);
    GRBQuadExpr g1 = -0.25 * y(0) * y(0) - y(1) * y(1) + 1;
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");
    GRBLinExpr g2 = y(0) - 2 * y(1) + 1;
    model.addConstr(g2, GRB_EQUAL, 0.0, "g2");

    model.update();
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

  bool is_feasible(xt::xarray<double>& x) {
    if (0.25 * pow(x(0), 2) - pow(x(1), 2) + 1 < -params.EPS) {
      return false;
    }

    if (fabs(x(0) - 2 * x(1) + 1.0) > params.EPS) {
      return false;
    }
    return true;
  }
};

class ModelProblem28 {
 public:
  GRBEnv& env;
  Param& params;
  std::string name = "Covariates";
  ModelProblem28(GRBEnv& a_env, Param& a_param) : env(a_env), params(a_param) {}

  double f(xt::xarray<double>& x) {
    return pow(x(0) + x(1), 2) + pow(x(1) + x(2), 2);
  }
  xt::xarray<double> gradient(xt::xarray<double>& x) {
    xt::xarray<double> gradient = xt::zeros<double>(x.shape());
    gradient(0) = 2 * (x(0) + x(1));
    gradient(1) = 2 * (x(0) + x(1)) + 2 * (x(1) + x(2));
    gradient(2) = 2 * (x(1) + x(2));

    return gradient;
  }
  xt::xarray<double> projection(xt::xarray<double>& x) {
    GRBModel model(env);
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});
    std::stringstream name;
    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      y(i) = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                          name.str());
      name.str("");
    }

    GRBQuadExpr obj;

    for (ulong i = 0; i < n; ++i) {
      obj += 0.5 * y(i) * y(i) - x(i) * y(i);
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBLinExpr g1 = y(0) + 2 * y(1) + 3 * y(2) - 1;
    model.addConstr(g1, GRB_EQUAL, 0.0, "g1");

    model.update();
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_DualReductions, 0);
    model.optimize();

    auto status = model.get(GRB_IntAttr_Status);
    xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());

    if (status == GRB_OPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }

    return y_val;
  }
  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }
  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    GRBModel model(env);
    std::stringstream name;
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});

    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      y(i) = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                          name.str());
      name.str("");
    }

    GRBLinExpr obj = 0;
    for (ulong i = 0; i < n; ++i) {
      obj += grad(i) * (y(i) - x(i));
    }

    model.setObjective(obj, GRB_MINIMIZE);
    GRBLinExpr g1 = y(0) + 2 * y(1) + 3 * y(2) - 1;
    model.addConstr(g1, GRB_EQUAL, 0.0, "g1");

    model.update();
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
  bool is_feasible(xt::xarray<double>& x) {
    return abs(x(0) + 2 * x(1) + 3 * x(2) - 1.0) <= params.EPS;
  }
};

class ModelProblem30 {
 public:
  GRBEnv& env;
  Param& params;
  std::string name = "Covariates";
  ModelProblem30(GRBEnv& a_env, Param& a_param) : env(a_env), params(a_param) {}

  double f(xt::xarray<double>& x) {
    return pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2);
  }
  xt::xarray<double> gradient(xt::xarray<double>& x) { return 2 * x; }
  xt::xarray<double> projection(xt::xarray<double>& x) {
    GRBModel model(env);
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});
    std::stringstream name;
    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      if (i == 0) {
        y(i) = model.addVar(1, 10, 0, GRB_CONTINUOUS, name.str());
      } else {
        y(i) = model.addVar(-10, 10, 0, GRB_CONTINUOUS, name.str());
      }
      name.str("");
    }

    GRBQuadExpr obj;

    for (ulong i = 0; i < n; ++i) {
      obj += 0.5 * y(i) * y(i) - x(i) * y(i);
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBQuadExpr g1 = y(0) * y(0) + y(1) * y(1) - 1;
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");

    model.update();
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_DualReductions, 0);
    model.optimize();

    auto status = model.get(GRB_IntAttr_Status);
    xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());

    if (status == GRB_OPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }

    return y_val;
  }
  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }
  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    GRBModel model(env);
    std::stringstream name;
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});

    for (ulong i = 0; i < n; ++i) {
      name << "y_" << i;
      if (i == 0) {
        y(i) = model.addVar(1, 10, 0, GRB_CONTINUOUS, name.str());
      } else {
        y(i) = model.addVar(-10, 10, 0, GRB_CONTINUOUS, name.str());
      }
      name.str("");
    }

    GRBLinExpr obj = 0;
    for (ulong i = 0; i < n; ++i) {
      obj += grad(i) * (y(i) - x(i));
    }

    model.setObjective(obj, GRB_MINIMIZE);
    GRBQuadExpr g1 = y(0) * y(0) + y(1) * y(1) - 1;
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");

    model.update();
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
  bool is_feasible(xt::xarray<double>& x) {
    if (pow(x(0), 2) + pow(x(1), 2) - 1 < -params.EPS) {
      return false;
    }
    if (x(0) < 1 - params.EPS || x(0) > 10 + params.EPS) {
      return false;
    }
    if (fabs(x(1)) > 10 + params.EPS || fabs(x(2)) > 10 + params.EPS) {
      return false;
    }
    return true;
  }
};

class ModelProblem31 {
 public:
  GRBEnv& env;
  Param& params;
  std::string name = "Covariates";
  ModelProblem31(GRBEnv& a_env, Param& a_param) : env(a_env), params(a_param) {}

  double f(xt::xarray<double>& x) {
    return 9 * pow(x(0), 2) + pow(x(1), 2) + 9 * pow(x(2), 2);
  }
  xt::xarray<double> gradient(xt::xarray<double>& x) {
    return xt::xarray<double>({18 * x(0), 2 * x(1), 18 * x(2)});
  }
  xt::xarray<double> projection(xt::xarray<double>& x) {
    GRBModel model(env);
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});
    y(0) = model.addVar(-10, 10, 0, GRB_CONTINUOUS, "y_0");
    y(1) = model.addVar(1, 10, 0, GRB_CONTINUOUS, "y_1");
    y(2) = model.addVar(-10, 1, 0, GRB_CONTINUOUS, "y_2");

    GRBQuadExpr obj;

    for (ulong i = 0; i < n; ++i) {
      obj += 0.5 * y(i) * y(i) - x(i) * y(i);
    }
    model.setObjective(obj, GRB_MINIMIZE);

    GRBQuadExpr g1 = (1.0 / 4) * ((x(0) + x(1)) * (x(0) + x(1)) -
                                  (x(1) - x(0)) * (x(1) - x(0)));
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");

    model.update();
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_DualReductions, 0);
    model.optimize();

    auto status = model.get(GRB_IntAttr_Status);
    xt::xarray<double> y_val = GRB_INFINITY * xt::ones<double>(y.shape());

    if (status == GRB_OPTIMAL) {
      for (ulong i = 0; i < n; ++i) {
        y_val(i) = y(i).get(GRB_DoubleAttr_X);
      }
    } else {
      cout << "Error. Projection problem solved with status = " << status
           << "\n";
      model.write("projection_model.lp");
      cout << "Wrote model at projection_model.lp\n";
      exit(status);
    }

    return y_val;
  }
  double get_rhs(xt::xarray<double>& grad, xt::xarray<double>& dir) {
    return rhs(grad, dir);
  }
  double get_lower_bound(xt::xarray<double>& x, xt::xarray<double>& grad) {
    GRBModel model(env);
    std::stringstream name;
    ulong n = x.shape(0);
    xt::xarray<GRBVar> y({n});

    y(0) = model.addVar(-10, 10, 0, GRB_CONTINUOUS, "y_0");
    y(1) = model.addVar(1, 10, 0, GRB_CONTINUOUS, "y_1");
    y(2) = model.addVar(-10, 1, 0, GRB_CONTINUOUS, "y_2");
    name.str("");

    GRBLinExpr obj = 0;
    for (ulong i = 0; i < n; ++i) {
      obj += grad(i) * (y(i) - x(i));
    }

    model.setObjective(obj, GRB_MINIMIZE);
    GRBQuadExpr g1 = (1.0 / 4) * ((x(0) + x(1)) * (x(0) + x(1)) -
                                  (x(1) - x(0)) * (x(1) - x(0)));
    model.addQConstr(g1, GRB_GREATER_EQUAL, 0.0, "g1");

    model.update();
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
  bool is_feasible(xt::xarray<double>& x) {
    if (x(0) * x(1) - 1 < -params.EPS) {
      return false;
    }

    if (fabs(x(0)) > 10 + params.EPS) {
      return false;
    }

    if (x(1) < 1 - params.EPS || x(1) > 10 + params.EPS) {
      return false;
    }

    if (x(2) < -10 - params.EPS || x(1) > 1 + params.EPS) {
      return false;
    }

    return true;
  }
};
#endif

#endif
