/*******************************************************/
/*                                                     */
/* A solver for constrained Hamiltonian systems        */
/* by MATH-CCES at RWTH Aachen                         */
/*                                                     */
/* *****************************************************/
   


// the base class for the Hamiltonian function
class ConstrainedHamiltonian
{
public:
  virtual int Dim () const = 0;
  virtual int DimConstraints () const = 0;
  virtual double Eval (double t, const Vector<> & q, const Vector<> & p) const = 0;
  virtual void EvalDq (double t, const Vector<> & q, const Vector<> & p, Vector<> & Dq) const = 0;
  virtual void EvalDp (double t, const Vector<> & q, const Vector<> & p, Vector<> & Dp) const = 0;
  virtual void EvalDpDp (double t, const Vector<> & q, const Vector<> & p, Matrix<> & DpDp) const = 0;

  virtual void Constraint (double t, const Vector<> & q, const Vector<> & g) const = 0;
  virtual void ConstraintDq (double t, const Vector<> & q, Matrix<> & Dq) const = 0;

  virtual void Redraw (const Vector <> & q) const { cout << "redraw not implemented" << endl; }
};




// a single-step method
class SSM
{
public:
  // do the step
  virtual void Step (double t, double h, 
                     const ConstrainedHamiltonian & func, 
                     const Vector<> & qold, Vector<> & qnew,
                     const Vector<> & pold, Vector<> & pnew) const = 0;
};


// the rattle algorithm
// following 'Geometric numerical integration illustrated by the Störmer/Verlet method'
// by E. Haier, Ch. Lubich and G. Wanner
// in Acta Numerica (2003), p1-51

class Rattle : public SSM
{
public:
  virtual void Step (double t, double h, const ConstrainedHamiltonian & func,
                     const Vector<> & qold, Vector<> & qnew,
                     const Vector<> & pold, Vector<> & pnew) const;
};





// the solver
class ODESolver 
{
  const ConstrainedHamiltonian & func;
  const SSM & ssm;
  double t, tend, h;
  Vector<> vecq, vecp;
public:
  ODESolver (const ConstrainedHamiltonian & afunc, const SSM & assm, 
             double t0, Vector<> & q0, Vector<> & p0, double atend, double ah)
    : func(afunc), ssm(assm),
      t(t0), tend(atend), h(ah),
      vecq(q0.Size()), vecp(q0.Size())
  { 
    vecq = q0;
    vecp = p0; 
  }

  void Solve()
  {
    while (t < tend)
      {
        Step();
        t += h;
      }
  }

  void Step()
  {
    Vector<> qnew(vecq.Size());
    Vector<> pnew(vecp.Size());
   
    cout << "\rt = " << setw(8) << t << ", H = " << func.Eval(t, vecq, vecp) << "                 " << flush;
    ssm.Step (t, h, func, vecq, qnew, vecp, pnew);
    vecq = qnew;
    vecp = pnew;
    t += h;
  }

  void Redraw () 
  {
    func.Redraw(vecq);
  }
};

