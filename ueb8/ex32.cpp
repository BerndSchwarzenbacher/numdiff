#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ode/RK_orig.hpp"


class Pendulum : public ODE_Function
{
private:
  double _mass;
  double _gravity;
  double _length;
public:
  Pendulum (double mass, double gravity, double length)
    : _mass(mass), _gravity(gravity), _length(length)
  {}

  virtual void Eval (double time, const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    double p_norm = y_old(2) * y_old(2) + y_old(3) * y_old(3);
    double q_norm = y_old(0) * y_old(0) + y_old(1) * y_old(1);
    double lambda = (_mass * _gravity * y_old(1) - 1.0/_mass * p_norm) / q_norm;
    y_new(0) = 1.0/_mass * y_old(2);
    y_new(1) = 1.0/_mass * y_old(3);
    y_new(2) = lambda * y_old(0);
    y_new(3) = - _mass * _gravity + lambda * y_old(1);
  }
};


int main ()
{
  ClassicalRK classical_rk;

  double t0 = 0;
  double t_end = 100;
  double step = 1e-4;
  int save_every = 1e3;

  Pendulum func(1, 9.81, 1);
  ngbla::Vector<> y0(4);
  y0(0) = 0;
  y0(1) = 1;
  y0(2) = 1;
  y0(3) = 0;

  std::ofstream out1("data/ex32.out");
  ODESolver (func, classical_rk, t0, y0, t_end, step, save_every, out1);

  return 0;
}

