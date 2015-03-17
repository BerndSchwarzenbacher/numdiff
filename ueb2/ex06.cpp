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

  virtual void Eval (double time, const ngbla::Vector<> & y,
                     ngbla::Vector<> & f) const
  {
    f(0) = y(1);
    f(1) = - _gravity / _length * sin(y(0));
  }

};


int main ()
{
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;
  ImplicitEuler impl_euler;

  std::ofstream out("ex06.out");
  Pendulum_ODE_Function func(1, 1, 1);
  ngbla::Vector<> y0(2);
  y0(0) = M_PI/4;
  y0(1) = 0;
  ODESolver (func, impl_euler, 0, y0, 1, 1e-3, out);

  //ofstream out2("mass_spring.out");
  //MassSpring_ODE_Function ms(10, 1);
  //Vector<> y0ms(2);
  //y0ms(0) = 1.0;
  //y0ms(1) = 0;
  //ODESolver (ms, impr_euler, 0, y0ms, 1000, 0.1, out2);

  return 0;
}

