#ifndef __IMPLICIT_EULER__
#define __IMPLICIT_EULER__

#include "../ode/bla/bla.hpp"
#include "../ode/ode.hpp"

class ImplicitEuler : public SSM
{
public:
  virtual void Step (double time, double step,
                     const ODE_Function & ode_func,
                     const ngbla::Vector<> & y_old,
                     ngbla::Vector<> & y_new) const
  {
    const double size = y_old.Size();
    ngbla::Matrix<double> f_jacobi(size, size);
    ngbla::Matrix<double> f_inverse(size, size);
    ngbla::Matrix<double> f_step(size, size);
    ngbla::Matrix<double> id = ngbla::Identity(size);
    ngbla::Vector<> f_eval(size);

    double qt = 0.0;
    double tol = 1e-2;

    ngbla::Vector<> y_delta(size);
    ngbla::Vector<> y_delta_old(size);

    // y_old should be a good start value for y
    ode_func.EvalDfDy(time, y_old, f_jacobi);
    f_step = step * f_jacobi - id;
    //std::cout << "F_step: " << std::endl << f_step << std::endl;
    ngbla::CalcInverse(f_step, f_inverse);
    //std::cout << "F_inverse: " << std::endl << f_inverse << std::endl;

    ode_func.Eval(time, y_old, f_eval);
    f_eval = step * f_eval;
    y_delta = - f_inverse * f_eval;
    y_new = y_old + y_delta;

    do
    {
      ode_func.EvalDfDy(time, y_new, f_jacobi);
      f_step = step * f_jacobi - id;
      ngbla::CalcInverse(f_step, f_inverse);

      ode_func.Eval(time, y_new, f_eval);
      f_eval = y_old - y_new + step * f_eval;

      y_delta_old = y_delta;
      y_delta = -f_inverse * f_eval;
      y_new += y_delta;
      qt = L2Norm(y_delta) / L2Norm(y_delta_old);

      if (qt >= 1)
      {
        std::cout << "Newton-Verfahren konvergiert nicht!" << std::endl;
        break;
      }
    }
    while ((qt / (1 - qt) * L2Norm(y_delta)) > tol);
  }

};

#endif

