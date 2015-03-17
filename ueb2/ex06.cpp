#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "implicit_euler.hpp"


class Impl_Pendulum_ODE_Function : public ODE_Function
{
private:
  double _length;
  double _mass;
  double _gravity;
  double _step;

public:
  Impl_Pendulum_ODE_Function (double length, double mass, double gravity,
                              double step)
    : _length(length), _mass(mass), _gravity(gravity), _step(step)
  { }

  virtual void Eval (double time, const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    y_new(0) = y_new(0) - y_old(0) - _step * y_old(1);
    y_new(1) = y_new(1) - y_old(1) - _step * sin(y_old(0));
  }

};


int main ()
{
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;
  ImplicitEuler impl_euler;

  double step = 1e-3;
  double t0 = 0;
  double t_end = 1;

  std::ofstream out("ex06.out");
  Impl_Pendulum_ODE_Function func(1, 1, 1, step);
  ngbla::Vector<> y0(2);
  y0(0) = M_PI/4;
  y0(1) = 0;
  ODESolver (func, impl_euler, t0, y0, t_end, step, out);

  //ofstream out2("mass_spring.out");
  //MassSpring_ODE_Function ms(10, 1);
  //Vector<> y0ms(2);
  //y0ms(0) = 1.0;
  //y0ms(1) = 0;
  //ODESolver (ms, impr_euler, 0, y0ms, 1000, 0.1, out2);

  return 0;
}

