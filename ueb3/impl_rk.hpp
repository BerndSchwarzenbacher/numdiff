#ifndef __IMPL_RK__
#define __IMPL_RK__

#include "../ode/ode.hpp"
#include "../ode/bla/bla.hpp"
#include "../ode/RK_orig.hpp"

class ImplicitRKMethod : public RungeKuttaMethod
{
public:
  ImplicitRKMethod (int as) : RungeKuttaMethod (as) { ; }

  virtual void Step ( double time, double step, const ODE_Function & func,
                      const ngbla::Vector<> & yold, ngbla::Vector<> & ynew )
    const
  {
    int n = yold.Size();
    int N = n * stages;

    ngbla::Vector<> all_ki_old(N);
    ngbla::Vector<> all_ki_new(N);
    all_ki_old = 0.0;
    all_ki_new = 0.0;

    ngbla::Matrix<double> dfdy(n, n);
    ngbla::Matrix<double> dFdy(N, N);
    ngbla::Matrix<double> dFdy_inverse(N, N);
    ngbla::Matrix<double> id = ngbla::Identity(N);
    ngbla::Vector<> F_eval(N);
    ngbla::Vector<> f_eval(n);

    double qt = 0.0;
    double tol = 1e-4;

    ngbla::Vector<> k_delta(N);
    ngbla::Vector<> k_delta_old(N);

    for (int i = 0; i < stages; ++i)
    {
      double ti = time + step * c(i);

      ngbla::Vector<> yi(n);
      yi = yold;
      for (int l = 0; l < stages; ++l)
        yi += step * A(i, l) * all_ki_old.Range(l*n, (l+1)*n);

      func.Eval(ti, yi, f_eval);
      F_eval.Range(i*n, (i+1)*n) = f_eval;
      func.EvalDfDy(ti, yi, dfdy);

      for (int ki_entry = 0; ki_entry < n; ++ki_entry)
      {
        for (int l = 0; l < stages; ++l)
        {
          for (int kl_entry = 0; kl_entry < n; ++kl_entry)
          {
            dFdy(i*n + ki_entry, l*n + kl_entry)
              = A(i, l) * all_ki_old(l*n+kl_entry);
          }
        }
      }

    }
    F_eval -= all_ki_old;
    dFdy = step * dFdy - id;

    ngbla::CalcInverse(dFdy, dFdy_inverse);

    k_delta = dFdy_inverse * F_eval;
    all_ki_new = all_ki_old - k_delta;

    int cnt=0;
    do
    {
      all_ki_old = all_ki_new;

      for (int i = 0; i < stages; ++i)
      {
        double ti = time + step * c(i);

        ngbla::Vector<> yi(n);
        yi = yold;
        for (int l = 0; l < stages; ++l)
          yi += step * A(i, l) * all_ki_old.Range(l*n, (l+1)*n);

        func.Eval(ti, yi, f_eval);
        F_eval.Range(i*n, (i+1)*n) = f_eval;
        func.EvalDfDy(ti, yi, dfdy);

        for (int ki_entry = 0; ki_entry < n; ++ki_entry)
        {
          for (int l = 0; l < stages; ++l)
          {
            for (int kl_entry = 0; kl_entry < n; ++kl_entry)
            {
              dFdy(i*n + ki_entry, l*n + kl_entry)
                = A(i, l) * all_ki_old(l*n+kl_entry);
            }
          }
        }

      }
      F_eval -= all_ki_old;
      dFdy = step * dFdy - id;

      ngbla::CalcInverse(dFdy, dFdy_inverse);

      k_delta_old = k_delta;
      k_delta = dFdy_inverse * F_eval;
      all_ki_new -= k_delta;

      qt = L2Norm(k_delta) / L2Norm(k_delta_old);
      if (qt >= 1 && cnt >= 10)
      {
        std::cout << "Newton-Verfahren konvergiert nicht!" << std::endl;
        break;
      }
      ++cnt;
    }
    while ((qt / (1 - qt) * L2Norm(k_delta)) > tol);

    ynew = yold;
    for (int i = 0; i < stages; i++)
      ynew += step * b(i) * all_ki_new.Range (i*n, (i+1)*n);
  }
};

class ImplicitMP : public ImplicitRKMethod
{
public:
  ImplicitMP() : ImplicitRKMethod (1)
  {
    ngbla::Matrix<> A(1,1);
    ngbla::Vector<> b(1), c(1);

    c = { 0.5 };
    b = { 1.0 };

    A = 0.5;

    SetAbc (A, b, c);
  }
};

class TwoStepGauss : public ImplicitRKMethod
{
public:
  TwoStepGauss() : ImplicitRKMethod (2)
  {
    ngbla::Matrix<> A(2,2);
    ngbla::Vector<> b(2), c(2);

    c = { 0.5 - sqrt(3) / 6.0, 0.5 - sqrt(3) / 6.0 };
    b = { 0.5, 0.5 };

    A = 0.25;
    A(0, 1) = 0.25 - sqrt(3) / 6.0;
    A(1, 0) = 0.25 + sqrt(3) / 6.0;

    SetAbc (A, b, c);
  }
};

class ImplicitEulerRK : public ImplicitRKMethod
{
public:
  ImplicitEulerRK() : ImplicitRKMethod (1)
  {
    ngbla::Matrix<> A(1,1);
    ngbla::Vector<> b(1), c(1);

    c = 1.0;
    b = 1.0;

    A = 1.0;

    SetAbc (A, b, c);
  }
};

#endif

