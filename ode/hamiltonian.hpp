#ifndef __HAMILTONIAN__
#define __HAMILTONIAN__
// a simple ODE - solver library
// Joachim Schoeberl

#include "bla/bla.hpp"


class Hamiltonian
{
public:
  virtual double Eval (const ngbla::Vector<> & p, const ngbla::Vector<> & q)
    const = 0;

  virtual void EvalDp (const ngbla::Vector<> & p, const ngbla::Vector<> & q,
                       ngbla::Vector<> & dHdp) const = 0;

  virtual void EvalDq (const ngbla::Vector<> & p, const ngbla::Vector<> & q,
                       ngbla::Vector<> & dHdq) const = 0;
};



// base class for the single-step method
class SSM_Hamilton
{
public:
  // do the step
  virtual bool Step (double t, double h, const Hamiltonian & func,
                     const ngbla::Vector<> & pold, ngbla::Vector<> & pnew,
                     const ngbla::Vector<> & qold, ngbla::Vector<> & qnew)
    const = 0;
};





// the time integration loop
void ODESolver_Hamilton (const Hamiltonian & func, const SSM_Hamilton & ssm,
                         double t0, ngbla::Vector<> & p0, ngbla::Vector<> & q0,
                         double tend, double h, int save_every,
                         std::ostream & out)
{
  double t = t0;
  int n = p0.Size();
  int j = 0;

  ngbla::Vector<> pold(n);
  ngbla::Vector<> pnew(n);
  ngbla::Vector<> qold(n);
  ngbla::Vector<> qnew(n);

  pold = p0;
  qold = q0;
  while (t < tend)
  {
    if ((j % save_every) == 0)
    {
      out << t;
      out << " " << func.Eval(pold, qold);
      for (int i = 0; i < n; ++i)
        out << " " << pold(i);
      for (int i = 0; i < n; ++i)
        out << " " << qold(i);
      out << std::endl;
    }

    ssm.Step (t, h, func, pold, pnew, qold, qnew);
    pold = pnew;
    qold = qnew;
    t += h;
    ++j;
  }
}

/* *************** Here are the specific single-step methods *************** */
class ExplicitSimplecticEuler : public SSM_Hamilton
{
public:
  virtual bool Step (double t, double h, const Hamiltonian & func,
                     const ngbla::Vector<> & pold, ngbla::Vector<> & pnew,
                     const ngbla::Vector<> & qold, ngbla::Vector<> & qnew)
    const
  {
    ngbla::Vector<> dHdp(qold.Size());
    ngbla::Vector<> dHdq(qold.Size());

    func.EvalDq (pold, qold, dHdq); // only for seperable Hamiltonian
    pnew = pold - h * dHdq;
    func.EvalDp (pnew, qold, dHdp);
    qnew = qold + h * dHdp;

    return true;
  }
};

class ExplicitStoermerVerlet : public SSM_Hamilton
{
public:
  virtual bool Step (double t, double h, const Hamiltonian & func,
                     const ngbla::Vector<> & pold, ngbla::Vector<> & pnew,
                     const ngbla::Vector<> & qold, ngbla::Vector<> & qnew)
    const
  {
    ngbla::Vector<> phalf(pold.Size());
    ngbla::Vector<> qhalf(qold.Size());

    ngbla::Vector<> dHdp(qold.Size());
    ngbla::Vector<> dHdq(qold.Size());

    func.EvalDq (pold, qold, dHdq); // only for seperable Hamiltonian
    phalf = pold - h/2.0 * dHdq;
    func.EvalDp (phalf, qold, dHdp);
    qhalf = qold + h/2.0 * dHdp;

    func.EvalDp (phalf, qold, dHdp); // only for seperable Hamiltonian
    qnew = qhalf + h/2.0 * dHdp;
    func.EvalDq (phalf, qnew, dHdq);
    pnew = phalf - h/2.0 * dHdq;

    return true;
  }
};

#endif

