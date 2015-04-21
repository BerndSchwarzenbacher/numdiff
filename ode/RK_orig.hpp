#ifndef __RK_ORIG__
#define __RK_ORIG__
// Runge Kutta methods
// Joachim Schoeberl


class RungeKuttaMethod : public SSM
{
protected:
  int stages;
  ngbla::Matrix<> A;
  ngbla::Vector<> b,c;
public:
  RungeKuttaMethod (int s)
    : stages(s), A(s,s), b(s), c(s)
  { ; }

  void SetAbc (const ngbla::Matrix<> & aa,
               const ngbla::Vector<> & ab,
               const ngbla::Vector<> & ac)
  {
    A = aa;
    b = ab;
    c = ac;
  }
};



class ExplicitRKMethod : public RungeKuttaMethod
{
public:
  ExplicitRKMethod (int as) : RungeKuttaMethod (as){ ; }


  virtual bool Step (double t, double h, const ODE_Function & func,
                     const ngbla::Vector<> & yold, ngbla::Vector<> & ynew) const
  {
    int n = yold.Size();

    ngbla::Vector<> yi(n), ki(n);
    ngbla::Vector<> all_ki(stages*n);
    all_ki = 0.0;

    for (int i = 0; i < stages; i++)
    {
      double ti = t + h * c(i);
      yi = yold;
      for (int j = 0; j < i; j++)
        yi += h * A(i,j) * all_ki.Range (j*n, (j+1)*n);
      func.Eval (ti, yi, ki);
      all_ki.Range(i*n, (i+1)*n) = ki;
    }

    ynew = yold;
    for (int i = 0; i < stages; i++)
      ynew += h * b(i) * all_ki.Range (i*n, (i+1)*n);

    return true;
  }
};



class ClassicalRK : public ExplicitRKMethod
{
public:
  ClassicalRK() : ExplicitRKMethod (4)
  {
    ngbla::Matrix<> A(4,4);
    ngbla::Vector<> b(4), c(4);

    c = { 0, 0.5, 0.5, 1 };

    b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

    A = 0.0;
    A(1,0) = 0.5;
    A(2,1) = 0.5;
    A(3,2) = 1;

    SetAbc (A, b, c);
  }

  virtual int Order () const { return 4; }
};

class Bogacki_Shampine3 : public ExplicitRKMethod
{
public:
  Bogacki_Shampine3() : ExplicitRKMethod(4)
  {
    ngbla::Matrix<> A(4,4);
    ngbla::Vector<> b(4), c(4);

    c = { 0, 0.5, 0.75, 1 };

    b = { 2.0/9, 1.0/3, 4.0/9, 0};

    A = 0.0;

    A(1,0) = 0.5;
    A(2,1) = 0.75;

    A(3,0) = 2.0/9;
    A(3,1) = 1.0/3;
    A(3,2) = 4.0/9;

    SetAbc (A, b, c);
  }

  virtual int Order () const { return 3; }
};


class ImprovedEulerRK : public ExplicitRKMethod
{
public:
  ImprovedEulerRK() : ExplicitRKMethod (2)
  {
    ngbla::Matrix<> A(2,2);
    ngbla::Vector<> b(2), c(2);

    c = { 0, 0.5 };
    b = { 0.0, 1.0 };

    A = 0.0;
    A(1,0) = 0.5;

    SetAbc (A, b, c);
  }
};

#endif

