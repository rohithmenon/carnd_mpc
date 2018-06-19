#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Reference velocity
double ref_v = 80;

// State variables and actuator variables start index
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double prev_throttle;
  double prev_steer;

  FG_eval(Eigen::VectorXd coeffs, double prev_steer, double prev_throttle) {
      this->coeffs = coeffs;
      this->prev_steer = prev_steer;
      this->prev_throttle = prev_throttle;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
      // The cost is stored is the first element of `fg`.
      // Any additions to the cost should be added to `fg[0]`.
      fg[0] = 0;

      // Cost for CTE and ECPI.
      for (int i = 0; i < N; i++) {
          fg[0] += CppAD::pow(vars[cte_start + i], 2);
          fg[0] += CppAD::pow(vars[epsi_start + i], 2);
          fg[0] += CppAD::pow(vars[v_start + i] - ref_v, 2);
      }

      // Cost for large values of actuations.
      for (int i = 0; i < N - 1; i++) {
          auto delta_val = CppAD::pow(vars[delta_start + i], 2);
          fg[0] += delta_val;
          fg[0] += CppAD::pow(vars[a_start + i], 2);
          // When steering value is high, apply a penalty for high speed.
          auto steer_factor = 1.0 / (1.0 + CppAD::exp(-25 * (delta_val - 0.05)));
          fg[0] += CppAD::pow(steer_factor * vars[v_start+i], 2);
      }

      // Cost for high difference in consecutive values
      for (int i = 0; i < N - 2; i++) {
          fg[0] += 2000 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
          fg[0] += 10 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
      }

      // Cost for high difference in previous steering actuation value and the first value of predicted
      // steering actuation.
      fg[0] += 2000 * CppAD::pow(vars[delta_start] - prev_steer, 2);
      //
      // Setup Constraints100
      //
      // NOTE: In this section you'll setup the model constraints.

      // Initial constraints
      //
      // We add 1 to each of the starting indices due to cost being located at
      // index 0 of `fg`.
      // This bumps up the position of all the other values.
      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];

      // The rest of the constraints
      for (int t = 1; t < N; t++) {
          // The state at time t+1 .
          AD<double> x1 = vars[x_start + t];
          AD<double> y1 = vars[y_start + t];
          AD<double> psi1 = vars[psi_start + t];
          AD<double> v1 = vars[v_start + t];
          AD<double> cte1 = vars[cte_start + t];
          AD<double> epsi1 = vars[epsi_start + t];

          // The state at time t.
          AD<double> x0 = vars[x_start + t - 1];
          AD<double> y0 = vars[y_start + t - 1];
          AD<double> psi0 = vars[psi_start + t - 1];
          AD<double> v0 = vars[v_start + t - 1];
          AD<double> cte0 = vars[cte_start + t - 1];
          AD<double> epsi0 = vars[epsi_start + t - 1];

          // Only consider the actuation at time t.
          AD<double> delta0 = vars[delta_start + t - 1];
          AD<double> a0 = vars[a_start + t - 1];

          AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
          AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * x0 * coeffs[2] + 3 * x0 * x0 * coeffs[3]);

          // Here's `x` to get you started.
          // The idea here is to constraint this value to be 0.
          //
          // Recall the equations for the model:
          // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
          // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
          // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
          // v_[t+1] = v[t] + a[t] * dt
          // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
          // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
          fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
          fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
          fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
          fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
          fg[1 + cte_start + t] =
                  cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
          fg[1 + epsi_start + t] =
                  epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
      }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

// Calculate the next state after dt seconds.
Eigen::VectorXd StateAt(const Eigen::VectorXd& state, double steer, double throttle, double dt) {
    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    Eigen::VectorXd new_state(6);
    new_state << x + v * cos(psi) * dt,
                 y + v * sin(psi) * dt,
                 psi + v * steer / Lf * dt,
                 v + throttle * dt,
                 cte + v * sin(psi) * dt,
                 epsi + v * steer / Lf * dt;
    return new_state;
}

vector<double> MPC::Solve(Eigen::VectorXd cur_state, Eigen::VectorXd coeffs) {
    typedef CPPAD_TESTVECTOR(double) Dvector;

    // Compute the state after latency seconds. 0.1 seconds = 100millis
    Eigen::VectorXd state = StateAt(cur_state, prev_steer, prev_throttle, 0.1);

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    // number of independent variables
    // N timesteps == N - 1 actuations
    size_t n_vars = N * 6 + (N - 1) * 2;
    // Number of constraints
    size_t n_constraints = N * 6;

    // Initial value of the independent variables.
    // Should be 0 except for the initial values.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++) {
        vars[i] = 0.0;
    }
    // Set the initial variable values
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

    // Lower and upper limits for x
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (int i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    // NOTE: Feel free to change this to something else.
    for (int i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits.
    // NOTE: Feel free to change this to something else.
    for (int i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Lower and upper limits for constraints
    // All of these should be 0 except the initial
    // state indices.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (int i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // Object that computes objective and constraints
    FG_eval fg_eval(coeffs, prev_steer, prev_throttle);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    //
    // Check some of the solution values
    //
    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    vector<double> result;
    double new_steer = solution.x[delta_start];
    double new_throttle = solution.x[a_start];
    prev_steer = new_steer;
    prev_throttle = new_throttle;
    result.push_back(new_steer);
    result.push_back(new_throttle);

    for (int i = 0; i < N-1; ++i) {
        result.push_back(solution.x[x_start + i + 1]);
        result.push_back(solution.x[y_start + i + 1]);
    }
    return result;
}
