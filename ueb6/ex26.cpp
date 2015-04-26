#include "../ode/bla/bla.hpp"
#include "../ode/hamiltonian.hpp"


class Hamiltonian_ODE : public Hamiltonian
{
public:
  virtual double Eval (const ngbla::Vector<> & p, const ngbla::Vector<> & q)
    const
  {
    return 1.0/2.0 * p(0)*p(0) - cos(q(0));
  }

  virtual void EvalDp (const ngbla::Vector<> & p, const ngbla::Vector<> & q,
                       ngbla::Vector<> & dHdp) const
  {
    dHdp(0) = p(0);
  }

  virtual void EvalDq (const ngbla::Vector<> & p, const ngbla::Vector<> & q,
                       ngbla::Vector<> & dHdq) const
  {
    dHdq(0) = sin(q(0));
  }
};


int main ()
{
  ExplicitSimplecticEuler simpl_euler;
  ExplicitStoermerVerlet expl_sv;

  double t0 = 0;
  double tend = 20;
  double h = 1e-4;
  double save_every = 1e2;

  Hamiltonian_ODE func;
  ngbla::Vector<> p0(1);
  p0(0) = 1;
  ngbla::Vector<> q0(1);
  q0(0) = 0;

  std::ofstream out1("data/ex26_euler.out");
  ODESolver_Hamilton (func, simpl_euler, t0, p0, q0, tend, h, save_every, out1);

  std::ofstream out2("data/ex26_expl_sv.out");
  ODESolver_Hamilton (func, expl_sv, t0, p0, q0, tend, h, save_every, out2);

  return 0;
}

