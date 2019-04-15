# MPC Controller

This project consist in model and programming a MPC controller to control the steering angle and speed of a car.

## MPC Controller

The MPC Controller is controller that predicts the possible outcomes and find
the best solution for the next step of the process.

#### Transform from world coordinates to car coordinates
First step is to transform the points from the world coordinate to car coordinate system.

```Cpp
for (size_t i = 0; i < ptsx.size(); i++) {
    double ptx_diff = ptsx[i] - px;
    double pty_diff = ptsy[i] - py;
    ptsx_car.push_back(ptx_diff * cos(psi) + pty_diff * sin(psi));
    ptsy_car.push_back(pty_diff * cos(psi) - ptx_diff * sin(psi));
}

Eigen::VectorXd ptx = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        ptsx_car.data(), static_cast<long>(ptsx_car.size()));
Eigen::VectorXd pty = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        ptsy_car.data(), static_cast<long>(ptsy_car.size()));
```

 
#### Find Coefficients of the polynomial

We find the coefficients of a 3rd order polynomial with the polyfit functions.

```Cpp
Eigen::VectorXd coeffs = polyfit(ptx, pty, 3);
```

#### Defining initial state

```Cpp
// Initial state.
double x0 = 0;
double y0 = 0;
double psi0 = 0;
double cte0 = coeffs[0];
double epsi0 = -atan(coeffs[1]);
```
#### Calculating the position after 100ms (dt), to anticipate the 100ms delay of the system. 

```CPP
// State after delay.
double x_delay = x0 + ( v * cos(psi0) * dt );
double y_delay = y0 + ( v * sin(psi0) * dt );
double psi_delay = psi0 - ( v * delta * dt / Lf );
double v_delay = v + a * dt;
double cte_delay = cte0 + ( v * sin(epsi0) * dt );
double epsi_delay = epsi0 - ( v * atan( (3*coeffs[3]*pow(x_delay, 2) + (2*coeffs[2]*x_delay) + coeffs[1] )) * dt / Lf );
```

#### Calculating the cross track error and the heading error

```Cpp
double cte = polyeval(coeffs, x0) - x0;
double epsi =  atan(coeffs[1] + (2*coeffs[2]*x0) + (3*coeffs[3]*x0*x0));
```

#### Timestep Length and Elapsed Duration (N & dt)
1 second prediction seem to work well, tried different periods, but they all are very worse than this.

Working with N = 8 is also ok, doesn't make much difference, as normally the last points doesn't act much in practice.

```
size_t N = 10;
double dt = 0.1;
```

### The Model

#### Initial position in array

Defining the initial position in the array for each variable
```
static size_t x_start = 0;
static size_t y_start = x_start + N;
static size_t psi_start = y_start + N;
static size_t v_start = psi_start + N;
static size_t cte_start = v_start + N;
static size_t epsi_start = cte_start + N;
static size_t delta_start = epsi_start + N;
static size_t a_start = delta_start + N - 1;
```

#### Cost Functions

The cost function for the model is defined as taught in class, with the difference that
the cross track error has to be much more prominent on the cost value, so we scale with the
cte_weight variable. The psi cost weight is also defined, so the speed is not as important as
this other variable.

```
int cte_weight = 10000;
int psi_weight = 2000;

for (size_t t = 0; t < N; t++) {
  fg[0] += 10000*CppAD::pow(vars[cte_start + t], 2); // Make CTE more important on the cost
  fg[0] += 2000*CppAD::pow(vars[epsi_start + t], 2);
  fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2); // Maintain velocity
}

// Minimize the use of actuators.
for (size_t t = 0; t < N - 1; t++) {
  fg[0] += CppAD::pow(vars[delta_start + t], 2);
  fg[0] += CppAD::pow(vars[a_start + t], 2);
}

// Minimize the value gap between sequential actuations.
for (size_t t = 0; t < N - 2; t++) {
  fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
  fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
}
    
```
#### Model and constraints to update the predicted values

Here we model the constrains to update the predicted variables, again, as taught in class,
we set the equations equal to zero.

```
fg[1 + x_start] = vars[x_start];
fg[1 + y_start] = vars[y_start];
fg[1 + psi_start] = vars[psi_start];
fg[1 + v_start] = vars[v_start];
fg[1 + cte_start] = vars[cte_start];
fg[1 + epsi_start] = vars[epsi_start];

// The rest of the constraints
for (size_t t = 1; t < N; ++t) {
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
  AD<double> f0 = coeffs[3] * CppAD::pow(x0, 3) + coeffs[2] * CppAD::pow(x0, 2)+ coeffs[1] * x0 + coeffs[0];
  AD<double> psides0 = CppAD::atan( (3*coeffs[3]*CppAD::pow(x0,2)) + (2*coeffs[2]*x0) + coeffs[1]);

  // Here's `x` to get you started.
  // The idea here is to constraint this value to be 0.
  //
  // Recall the equations for the model:
  // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
  // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
  // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
  // v_[t] = v[t-1] + a[t-1] * dt
  // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
  // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
  fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
  fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
  fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
  fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
  fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
}
```

#### More constraints settings

We define the variables based on number of samples we will predict, and set
some initial values and upper and lower limits.

```
size_t n_vars =  N * 6 + (N - 1) * 2;
size_t n_constraints =  N * 6;

double x = state[0];
double y = state[1];
double psi = state[2];
double v = state[3];
double cte = state[4];
double epsi = state[5];

// Initial value of the independent variables.
// SHOULD BE 0 besides initial state.
Dvector vars(n_vars);
for (size_t i = 0; i < n_vars; ++i) {
vars[i] = 0;
}
// Set the initial variable values
vars[x_start] = x;
vars[y_start] = y;
vars[psi_start] = psi;
vars[v_start] = v;
vars[cte_start] = cte;
vars[epsi_start] = epsi;

Dvector vars_lowerbound(n_vars);
Dvector vars_upperbound(n_vars);

// Set all non-actuators upper and lowerlimits
// to the max negative and positive values.
for (size_t i = 0; i < delta_start; ++i) {
  vars_lowerbound[i] = -1.0e19;
  vars_upperbound[i] = 1.0e19;
}

// The upper and lower limits of delta are set to -25 and 25
// degrees (values in radians).
// NOTE: Feel free to change this to something else.
for (size_t i = delta_start; i < a_start; ++i) {
  vars_lowerbound[i] = -0.436332;
  vars_upperbound[i] = 0.436332;
}

// Acceleration/decceleration upper and lower limits.
// NOTE: Feel free to change this to something else.
for (size_t i = a_start; i < n_vars; ++i) {
  vars_lowerbound[i] = -1.0;
  vars_upperbound[i] = 1.0;
}


// Lower and upper limits for the constraints
// Should be 0 besides initial state.
Dvector constraints_lowerbound(n_constraints);
Dvector constraints_upperbound(n_constraints);
for (size_t i = 0; i < n_constraints; ++i) {
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
```

We use ipopt to solve the predictions and return the values on a standard vector.
 
```
// place to return solution
CppAD::ipopt::solve_result<Dvector> solution;

// solve the problem
CppAD::ipopt::solve<Dvector, FG_eval>(
    options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
    constraints_upperbound, fg_eval, solution);

// Check some of the solution values
ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

// Cost
auto cost = solution.obj_value;
std::cout << "Cost " << cost << std::endl;
std::cout << "Result " << solution.x << " : " << solution.x[a_start] << std::endl;
/**
* Returning the whole solution as a vector.
*/
std::vector<double> res;
for (size_t i = 0; i < solution.x.size(); ++i) {
res.push_back(solution.x[i]);
}

return res;
```

## Authors

* **Ian Esper**