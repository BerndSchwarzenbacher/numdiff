#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "implicit_euler.hpp"


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
  ExplicitEuler expl_euler;
  //ImprovedEuler impr_euler;
  ImplicitEuler impl_euler;

  double step = 1e-5;
  double t0 = 0;
  double t_end = 10000;
  int save_every = 5e4;

  std::ofstream out("ex06_impl_1e-5.out");
  Pendulum_ODE_Function func(1, 1, 1);
  ngbla::Vector<> y0(2);
  y0(0) = M_PI/4;
  y0(1) = 0;
  ODESolver (func, impl_euler, t0, y0, t_end, step, save_every, out);

  //std::ofstream out2("ex06_expl_euler.out");
  //ODESolver (func, expl_euler, t0, y0, t_end, step, save_every, out2);

  return 0;
}

