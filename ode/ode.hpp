#ifndef __ODE__
#define __ODE__
// a simple ODE - solver library
// Joachim Schoeberl

#include "bla/bla.hpp"


// the base class for the right-hand-side f(t,y)
class ODE_Function
{
public:
  // must be overloaded by derived class
  virtual void Eval (double t, const ngbla::Vector<> & y, ngbla::Vector<> & f)
    const = 0;

  virtual void EvalDfDy (double t,
                         const ngbla::Vector<> & y,
                         ngbla::Matrix<> & dfdy) const
  {
    // numerical differentiation
    int n = y.Size();
    ngbla::Vector<> yr(n), yl(n), fr(n), fl(n);
    double eps = 1e-6;
    for (int i = 0; i < n; i++)
    {
      yl = y;  yl(i) -= eps;
      yr = y;  yr(i) += eps;
      Eval (t, yl, fl);
      Eval (t, yr, fr);

      dfdy.Col(i) = 1.0/(2*eps) * (fr - fl);
    }
  }
};



// base class for the single-step method
class SSM
{
public:
  // do the step
  virtual void Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew)
    const = 0;
};





// the time integration loop
void ODESolver (const ODE_Function & func, const SSM & ssm,
                double t0, ngbla::Vector<> & y0, double tend, double h,
                std::ostream & out)
{
  double t = t0;
  int n = y0.Size();

  ngbla::Vector<> yold(n), ynew(n);

  yold = y0;
  while (t < tend)
  {
    out << t;
    for (int i = 0; i < n; ++i)
      out << " " << yold(i);
    out << std::endl;

    ssm.Step (t, h, func, yold, ynew);
    yold = ynew;
    t += h;
  }
}








/* *************** Here are the specific single-step methods *************** */



class ExplicitEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew) const
  {
    ngbla::Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h * f;
  }
};


class ImprovedEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew) const
  {
    ngbla::Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h/2.0 * f;

    func.Eval (t+h/2.0, ynew, f);
    ynew = yold + h * f;
  }
};

#endif

