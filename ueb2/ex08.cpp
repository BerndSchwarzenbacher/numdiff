#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "implicit_euler.hpp"


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
  ImplicitEuler impl_euler;

  double step = 1e-3;
  double t0 = 0;
  double t_end = 1000;
  int save_every = 5e2;

  Network_ODE_Function func(1, 0.001, 0.001, 1000, 1);
  ngbla::Vector<> y0(2);
  y0(0) = M_PI/4;
  y0(1) = 0;

  std::ofstream out("ex08_impl_omeg1_1e-3.out");
  ODESolver (func, impl_euler, t0, y0, t_end, step, save_every, out);

  std::ofstream out2("ex08_expl_omeg1_1e-3.out");
  ODESolver (func, expl_euler, t0, y0, t_end, step, save_every, out2);

  Network_ODE_Function func2(1000, 0.001, 0.001, 1000, 1);
  std::ofstream out3("ex08_impl_omeg1000_1e-3.out");
  ODESolver (func2, impl_euler, t0, y0, t_end, step, save_every , out3);

  std::ofstream out4("ex08_expl_omeg1000_1e-4.out");
  ODESolver (func2, expl_euler, t0, y0, t_end, 1e-4, 5e3, out4);

  return 0;
}

