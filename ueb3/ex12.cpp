#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ode/impl_rk.hpp"


class Pendulum_ODE_Function : public ODE_Function
{
private:
  double _length;
  double _mass;
  double _gravity;

public:
  Pendulum_ODE_Function (double length, double mass, double gravity)
    : _length(length), _mass(mass), _gravity(gravity)
  { }

  virtual void Eval (double time, const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    y_new(0) = y_old(1);
    y_new(1) = -sin(y_old(0));
  }

  virtual void EvalDfDy (double t,
                         const ngbla::Vector<> & y,
                         ngbla::Matrix<> & dfdy) const
  {
    dfdy = 0.0;
    dfdy(0,1) = 1.0;
    dfdy(1,0) = - cos(y(0));
  }

};


int main ()
{
  ImplicitMP impl_mp;
  TwoStepGauss two_step_gauss;
  double step = 1e-3;
  double t0 = 0;
  double t_end = 10000;
  int save_every = 2e2;

  std::ofstream out("data/ex12_impl_mp.out");
  Pendulum_ODE_Function func(1, 1, 1);
  ngbla::Vector<> y0(2);
  y0(0) = M_PI/4;
  y0(1) = 0;
  ODESolver (func, impl_mp, t0, y0, t_end, step, save_every, out);

  std::ofstream out2("data/ex12_two_step_gauss.out");
  ODESolver (func, two_step_gauss, t0, y0, t_end, step, save_every, out2);

  return 0;
}

