#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ueb2/implicit_euler.hpp"
#include "impl_rk.hpp"


class Network_ODE_Function : public ODE_Function
{
private:
  double _omega;
  double _cap;
  double _res;
  double _cur_s;
  double _volt_t;

public:
  Network_ODE_Function (double omega, double cap, double res, double cur_s,
                        double volt_t)
    : _omega(omega), _cap(cap), _res(res), _volt_t(volt_t)
  { }

  virtual void Eval (double time, const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    //double pot = (sin(_omega * time) - y_old(0)) / _volt_t;
    //y_new(0) = (_cur_s * (exp(pot) - 1) - y_old(0) / _res) / _cap;
    y_new(0) = 0.001 * exp((sin(time) - y_old(0))/0.2) - 0.001 - y_old(0);
  }

  virtual void EvalDfDy (double time,
                         const ngbla::Vector<> & y,
                         ngbla::Matrix<> & dfdy) const
  {
    //double pot = (sin(_omega * time) - y(0)) / _volt_t;
    //dfdy(0, 0) = -(_cur_s * exp(pot) / (_volt_t * _cap)) - 1.0 / (_res * _cap);

    dfdy(0, 0) = - 0.001 * exp((sin(time) - y(0))/0.2) * y(0) - 1;
  }

};


int main ()
{
  ImplicitEuler impl_euler;
  ImplicitEulerRK impl_euler_rk;
  ImplicitMP impl_mp;
  TwoStepGauss two_step_gauss;

  double t0 = 0;
  double t_end = 30;

  Network_ODE_Function func(1, 1, 1, 0.001, 0.2);
  ngbla::Vector<> y0(1);
  y0(0) = 0.0;

  std::ofstream out1("data/ex13_impl_mp_e-2.out");
  ODESolver (func, impl_mp, t0, y0, t_end, 1e-2, 1, out1);

  return 0;
}

