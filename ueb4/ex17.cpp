#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"


class Network_ODE_Function : public ODE_Function
{
private:
  double _omega;
  double _cap1;
  double _cap2;
  double _res1;
  double _res2;

public:
  Network_ODE_Function (double omega, double cap1, double cap2,
                        double res1, double res2)
    : _omega(omega), _cap1(cap1), _cap2(cap2), _res1(res1), _res2(res2)
  { }

  virtual void Eval (double time, const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    y_new(0) = (sin(_omega * time) - y_old(0) - _res1 * _cap2 * y_old(1)) /
               (_res1 * _cap1);
    y_new(1) = (y_old(0) - y_old(1)) / (_res2 * _cap2);
  }

  virtual void EvalDfDy (double t,
                         const ngbla::Vector<> & y,
                         ngbla::Matrix<> & dfdy) const
  {
    dfdy = 0.0;
    dfdy(0,0) = - 1 / (_res1 * _cap1);
    dfdy(0,1) = - _cap2 / _cap1;
    dfdy(1,0) =   1 / (_res2 * _cap2);
    dfdy(1,1) = - 1 / (_res2 * _cap2);
  }

};


int main ()
{
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;

  double t_0 = 0;
  double t_end = 10;
  double alpha_min = 0.2;
  double alpha_max = 2.0;
  double beta = 0.9;
  double h_start = 1e-3;
  double h_min = 1e-7;
  double h_max = 1;
  double save_every = 1e-2;

  Network_ODE_Function func(1, 0.001, 0.001, 1000, 1);
  ngbla::Vector<> y_0(2);
  y_0(0) = -1;
  y_0(1) = 1;

  std::ofstream out("data/ex17.out");
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, 100, y_0, 1e-3,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out);

  std::ofstream out1("data/ex17_tol_e2.out");
  std::ofstream out2("data/ex17_tol_e3.out");
  std::ofstream out3("data/ex17_tol_e4.out");
  std::ofstream out4("data/ex17_tol_e5.out");
  std::ofstream out5("data/ex17_tol_e6.out");
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, 1e-2,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out1);
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, 1e-3,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out2);
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, 1e-4,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out3);
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, 1e-5,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out4);
  ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, 1e-6,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out5);

  return 0;
}

