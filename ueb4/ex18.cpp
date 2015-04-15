#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ode/impl_rk.hpp"


class Chemical_Reaction5 : public ODE_Function
{
private:
  double _k1;
  double _k2;
  double _k3;
  double _k4;
  double _k5;

public:
  Chemical_Reaction5 ( double k1, double k2, double k3, double k4, double k5)
    : _k1(k1), _k2(k2), _k3(k3), _k4(k4), _k5(k5)
  { }

  virtual void Eval (double time, const ngbla::Vector<> & y,
                     ngbla::Vector<> & f) const
  {
    f(0) =  - _k1 * y(0) * y(1) - _k3 * y(0) * y(2)
            + _k4 * y(2) * y(2);

    f(1) =  - _k1 * y(0) * y(1) - _k2 * y(1) * y(2)
            + _k5 * y(4);

    f(2) =    _k1 * y(0) * y(1) - _k2 * y(1) * y(2)
            + _k3 * y(0) * y(2) - 2 * _k4 * y(2) * y(2);

    f(3) =    _k1 * y(0) * y(1) + _k2 * y(1) * y(2)
            + _k4 * y(2) * y(2);

    f(4) =    _k3 * y(0) * y(2) - _k5 * y(4);

  }

  virtual void EvalDfDy (double time,
                         const ngbla::Vector<> & y,
                         ngbla::Matrix<> & dfdy) const
  {
    dfdy = 0.0;

    dfdy(0,0) = - _k1 * y(1) - _k3 * y(2);
    dfdy(0,1) = - _k1 * y(0);
    dfdy(0,2) = - _k3 * y(0) + 2 * _k4 * y(2);

    dfdy(1,0) = - _k1 * y(1);
    dfdy(1,1) = - _k1 * y(0) - _k2 * y(2);
    dfdy(1,2) = - _k2 * y(1);
    dfdy(1,4) =   _k5;

    dfdy(2,0) =   _k1 * y(1) + _k3 * y(2);
    dfdy(2,1) =   _k1 * y(0) - _k2 * y(2);
    dfdy(2,2) = - _k2 * y(1) + _k3 * y(0) - 4 * _k4 * y(2);

    dfdy(3,0) =   _k1 * y(1);
    dfdy(3,1) =   _k1 * y(0) + _k2 * y(2);
    dfdy(3,2) =   _k2 * y(1) + 2 * _k4 * y(2);

    dfdy(4,0) =   _k3 * y(2);
    dfdy(4,2) =   _k3 * y(0);
    dfdy(4,4) = - _k5;
  }

};


int main ()
{
  ImplicitMP impl_mp;
  TwoStepGauss two_step;
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;

  double t_0 = 0;
  double t_end = 200;
  double tol = 1e-3;
  double alpha_min = 0.2;
  double alpha_max = 2.0;
  double beta = 0.9;
  double h_start = 1e-3;
  double h_min = 1e-7;
  double h_max = 1;
  double save_every = 1e-3;

  Chemical_Reaction5 func(1.34, 1.6*1e9, 8.0*1e3, 4.0*1e7, 1.0);
  ngbla::Vector<> y_0(5);
  y_0(0) = 0.05;
  y_0(1) = 1e-4;
  y_0(2) = 1e-10;
  y_0(3) = 0.1;
  y_0(4) = 1e-4;

  std::ofstream out("data/ex18.out");
  ODESolverAdaptive (func, impl_mp, two_step, t_0, t_end, y_0, tol,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out);

  //ODESolverAdaptive (func, expl_euler, impr_euler, t_0, t_end, y_0, tol,
                     //alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     //save_every, out);

  return 0;
}

