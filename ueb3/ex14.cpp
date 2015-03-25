#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ueb2/implicit_euler.hpp"
#include "impl_rk.hpp"


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
  ImplicitEuler impl_euler;
  ImplicitEuler impl_euler_rk;
  ImplicitMP impl_mp;

  double t0 = 0;
  double t_end = 1000;

  Network_ODE_Function func(1, 0.001, 0.001, 1000, 1);
  ngbla::Vector<> y0(2);
  y0(0) = -1;
  y0(1) = 1;

  std::ofstream out1("data/ex14_impl_eul_1.out");
  std::ofstream out2("data/ex14_impl_eul_e-1.out");
  std::ofstream out3("data/ex14_impl_eul_e-2.out");
  std::ofstream out4("data/ex14_impl_eul_e-3.out");
  ODESolver (func, impl_euler_rk, t0, y0, t_end, 1,    1, out1);
  ODESolver (func, impl_euler_rk, t0, y0, t_end, 1e-1, 1, out2);
  ODESolver (func, impl_euler_rk, t0, y0, t_end, 1e-2, 1, out3);
  ODESolver (func, impl_euler_rk, t0, y0, t_end, 1e-3, 1, out4);

  std::ofstream out5("data/ex14_impl_mp_1.out");
  std::ofstream out6("data/ex14_impl_mp_e-1.out");
  std::ofstream out7("data/ex14_impl_mp_e-2.out");
  std::ofstream out8("data/ex14_impl_mp_e-3.out");
  ODESolver (func, impl_mp, t0, y0, t_end, 1,    1, out5);
  ODESolver (func, impl_mp, t0, y0, t_end, 1e-1, 1, out6);
  ODESolver (func, impl_mp, t0, y0, t_end, 1e-2, 1, out7);
  ODESolver (func, impl_mp, t0, y0, t_end, 1e-3, 1, out8);

  return 0;
}

