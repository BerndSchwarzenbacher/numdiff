#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"
#include "../ode/RK_orig.hpp"


class Orbit_ODE_Function : public ODE_Function
{
private:
  double _mu1;
  double _mu2;

public:
  Orbit_ODE_Function (double mu)
    : _mu1(mu), _mu2(1-mu)
  { }

  virtual void Eval (double t, const ngbla::Vector<> & y,
                     ngbla::Vector<> & f) const
  {
    double x1_mu1 = y(0) + _mu1;
    double x1_mu2 = y(0) - _mu2;

    double n1 = pow((x1_mu1)*(x1_mu1) + y(2)*y(2), 3.0/2);
    double n2 = pow((x1_mu2)*(x1_mu2) + y(2)*y(2), 3.0/2);

    f(0) = y(1);
    f(1) = y(0) + 2*y(3) - _mu2 * x1_mu1 / n1 - _mu1 * x1_mu2 / n2;
    f(2) = y(3);
    f(3) = y(2) - 2*y(1) - _mu2 * y(2) / n1 - _mu1 * y(2) / n2;
  }
};


int main ()
{
  Bogacki_Shampine3 bog_sh3;
  ClassicalRK class_rk;

  double t_0 = 0;
  double t_end = 17.0652166;
  double tol = 1/t_end;
  double alpha_min = 0.2;
  double alpha_max = 2.0;
  double beta = 0.9;
  double h_start = 1e-6;
  double h_min = 1e-10;
  double h_max = 1e-2;
  double save_every = 1e-4;

  Orbit_ODE_Function func(0.012277471);
  ngbla::Vector<> y0(4);
  y0(0) = 0.994;
  y0(1) = 0;
  y0(2) = 0;
  y0(3) = -2.0015851063790825;

  std::ofstream out("data/ex22.out");
  ODESolverAdaptive (func, bog_sh3, class_rk, t_0, t_end, y0, tol,
                     alpha_min, alpha_max, beta, h_start, h_min, h_max,
                     save_every, out);

  return 0;
}

