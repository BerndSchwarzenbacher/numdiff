// compile with
// g++ ode_demo.cpp bla/calcinverse.cpp bla/exception.cpp bla/localheap.cpp -std=c++11 -o ode_demo -fmax-errors=2 -fdiagnostics-color=auto
//


// a simple linear algebra library
#include "bla/bla.hpp"
using namespace ngbla;

#include "ode.hpp"



/* ******************* The ODEs I want to solve  ******************** */


//  f(t,y) = a y
class My_First_ODE_Function : public ODE_Function
{
  double a;

public:
  My_First_ODE_Function (double aa) { a = aa; }

  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const
  {
    f = a * y;
  }
};





class MassSpring_ODE_Function : public ODE_Function
{
  double m;
  double k;
public:
  MassSpring_ODE_Function (double am, double ak) { m = am; k = ak; }

  // position, momentum
  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const
  {
    f(0) = 1.0/m * y(1);
    f(1) = -k * y(0);
  }
};







int main ()
{
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;

  ofstream out("data.out");
  My_First_ODE_Function func(1);
  Vector<> y0(1);  // initial conditions
  y0(0) = 1.0;
  ODESolver (func, expl_euler, 0, y0, 10, 0.1, out);


  ofstream out2("mass_spring.out");
  MassSpring_ODE_Function ms(10, 1);
  Vector<> y0ms(2);
  y0ms(0) = 1.0;
  y0ms(1) = 0;
  ODESolver (ms, impr_euler, 0, y0ms, 1000, 0.1, out2);

  return 0;
}
