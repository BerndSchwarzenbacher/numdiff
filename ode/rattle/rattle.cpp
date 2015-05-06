#include "../bla/bla.hpp"
using namespace ngbla;

#include "rattle.hpp"


/*
// the time integration loop
void ODESolver (const ConstrainedHamiltonian & func, 
const SSM & ssm, double t0, Vector<> & q0, Vector<> & p0, 
double tend, double h)
{
double t = t0;
int n = q0.Size();

Vector<> qold(n), qnew(n);
Vector<> pold(n), pnew(n);

qold = q0;
pold = p0;
int cnt = 0;
while (t < tend)
{
  if ( (cnt++ % 1) == 0)
    {
      cout << t;
      cout << " " << func.Eval (t, qold, pold);

      for (int i = 0; i < n; i++)
        cout << " " << qold(i);
      for (int i = 0; i < n; i++)
        cout << " " << pold(i);
      cout << endl;
    }
  ssm.Step (t, h, func, qold, qnew, pold, pnew);
  qold = qnew;
  pold = pnew;
  t += h;
 }
}
*/
                           
                           
                           
                           
                           
 /* ******************* Here are the specific single-step methods ******************** */

 /*
   We assume a seperated Hamiltonian H(p,q) = T(p) + V(q)  with quadratic kinetic energy T(q)
 */

void Rattle::Step (double t, double h, const ConstrainedHamiltonian & func,
                   const Vector<> & qold, Vector<> & qnew,
                   const Vector<> & pold, Vector<> & pnew) const
{
  int n = func.Dim();
  int m = func.DimConstraints();
  
  Vector<> pmid(n), dhdq(n), dhdp(n), lam(m), res(m);
  Matrix<> dgdq(m,n), dgdq2(m,n), dhdpdp(n), newtonmat(m,m), invnewton(m,m);


  func.EvalDq (t, qold, pold, dhdq);
  func.EvalDpDp (t, pold, qold, dhdpdp);

  func.ConstraintDq(t, qold, dgdq);
  lam = 0.0;

  int it = 0;
  do
    {
      if (it++ > 5) 
        {
          cerr << "Newton does not converge, decrease step-size" << endl;
          exit(1);
        }

      pmid = pold - h/2 * (dhdq + Trans (dgdq) * lam);

      func.EvalDp (t, qold, pmid, dhdp);
      qnew = qold + h * dhdp;
      func.Constraint (t, qnew, res);

      func.ConstraintDq(t, qnew, dgdq2);
      newtonmat = dgdq2 * dhdpdp * Trans (dgdq);
      newtonmat *= h*h/2;

      CalcInverse (newtonmat, invnewton);

      lam += invnewton * res;
    }
  while (L2Norm(res) > 1e-10);

  /*
  // checking step 1:
  Vector<> help(n), help2(n);

  func.EvalDq(t, pmid, qold, help);
  func.ConstraintDq (t, qold, dgdq);
  cout << "pmid = " << pmid << " =?= " << pold - h/2 * (help + Trans (dgdq) * lam);
    
  func.Constraint (t, qnew, res);
  cout << "constraint = " << res << endl;
  */
    
  func.EvalDq (t, qnew, pold, dhdq);    
  pnew = pmid - h/2 * dhdq;
  func.ConstraintDq(t, qnew, dgdq2);
  res = dgdq2 * dhdpdp * pnew;

  newtonmat = dgdq2 * dhdpdp * Trans (dgdq2);
  newtonmat *= h/2;
  CalcInverse (newtonmat, invnewton);
  lam = invnewton * res;

  pnew -= h/2 * Trans(dgdq2) * lam;

  /*
  // checking step 2:
  func.EvalDp (t, qnew, pnew, dhdp);
  cout << "pnew = " << pnew << " =?= " << pmid - h/2 * (help + Trans (dgdq2) * lam);
  cout << "constraint: = " << dgdq2 * dhdp << endl;
  */
}





