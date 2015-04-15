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
  virtual bool Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew)
    const = 0;
  virtual int Order () const { return -1; }
};





// the time integration loop
void ODESolver (const ODE_Function & func, const SSM & ssm,
                double t0, ngbla::Vector<> & y0, double tend, double h,
                int save_every, std::ostream & out)
{
  double t = t0;
  int n = y0.Size();
  int j = 0;

  ngbla::Vector<> yold(n), ynew(n);

  yold = y0;
  while (t < tend)
  {
    if ((j % save_every) == 0)
    {
      out << t;
      for (int i = 0; i < n; ++i)
        out << " " << yold(i);
      out << std::endl;
    }

    ssm.Step (t, h, func, yold, ynew);
    yold = ynew;
    t += h;
    ++j;
  }
}

// the time integration loop for adaptive
void ODESolverAdaptive (const ODE_Function & func, const SSM & ssm_low_order,
                        const SSM & ssm_high_order, double t_0, double t_end,
                        ngbla::Vector<> & y_0, double tol,
                        double alpha_min, double alpha_max, double beta,
                        double h_start, double h_min, double h_max,
                        double save_every, std::ostream & out )
{
  double t = t_0;
  double t_saved = t;
  int dim = y_0.Size();
  double h_new = h_start;
  double h = h_new;

  ngbla::Vector<> y_old(dim);
  ngbla::Vector<> y_new_low(dim);
  ngbla::Vector<> y_new_high(dim);

  int high_order = ssm_high_order.Order();

  double err_quot = 0;

  y_old = y_0;
  while (t < t_end)
  {
    if ((t - t_saved) >= save_every || t == t_0)
    {
      out << t;
      out << " " << h;
      for (int i = 0; i < dim; ++i)
        out << " " << y_old(i);
      out << std::endl;

      t_saved = t;
    }

    do
    {
      h = fmin (h_new, t_end - t);

      if (!ssm_low_order.Step (t, h, func, y_old, y_new_low)
          || !ssm_high_order.Step (t, h, func, y_old, y_new_high))
      {
        if ( h == h_min )
          throw ngstd::Exception();

        h_new = fmax (h / 2.0, h_min);
        continue;
      }

      double loc_err = L2Norm(y_new_low - y_new_high);
      err_quot = loc_err / (h * tol);

      double alpha = fmax (alpha_min, pow(err_quot, -1.0/high_order));
      alpha = fmin (alpha, alpha_max);

      h_new = fmax (h_min, beta * alpha * h);
      h_new = fmin (h_max, h_new);
    }
    while (err_quot > 1 && h != h_min);

    y_old = y_new_high;
    t += h;
  }

  out << t;
  out << " " << h;
  for (int i = 0; i < dim; ++i)
    out << " " << y_old(i);
  out << std::endl;

}







/* *************** Here are the specific single-step methods *************** */



class ExplicitEuler : public SSM
{
public:
  virtual bool Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew) const
  {
    ngbla::Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h * f;

    return true;
  }

  virtual int Order () const
  {
    return 1;
  }
};


class ImprovedEuler : public SSM
{
public:
  virtual bool Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew) const
  {
    ngbla::Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h/2.0 * f;

    func.Eval (t+h/2.0, ynew, f);
    ynew = yold + h * f;

    return true;
  }

  virtual int Order () const
  {
    return 2;
  }
};

#endif

