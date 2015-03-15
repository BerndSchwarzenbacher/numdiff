/*
  Basic Linear Algebra demos
*/

//g++ demo_bla.cpp bla/calcinverse.cpp bla/exception.cpp bla/localheap.cpp -std=c++11 -o demo_bla -fmax-errors=2 -fdiagnostics-color=auto

// ng-soft header files
#include "bla/bla.hpp"


//using namespace std;
using namespace ngbla;

int main ()
{
  Vector<double> u(3), v(5);
  Matrix<double> m(5,3);

  u = 1.0;

  for (int i = 0; i < Height(m); i++)
    for (int j = 0; j < Width(m); j++)
      m(i,j) = i+j;


  // do some math:
  v = m * u;

  cout << "u = " << endl << u << endl;
  cout << "v = " << endl << v << endl;
  cout << "m = " << endl << m << endl;

  // use result directly
  cout << "Trans(m) * v = " << endl << Trans(m) * v << endl;

  // fix-size objects:
  Vec<3,double> u3;
  Mat<3,3,double> m3;

  u3 = 2.7;
  m3 = 1.0;
  cout << "m3 * u3 = " << m3 * u3 << endl;

  {
    // C++11 : initializer lists
    Vector<> x = { 1, 2, 3 };
    Matrix<> m = { { 2, 2, 2 }, { 5, 5, 5 } };
    Vector<> y = m * x;
    cout << "x = " << x << endl;
    cout << "m = " << m << endl;
    cout << "y = " << y << endl;
  }


  // own memory management:
  double data[1000];

  FlatVector<double> fu(100, data);
  FlatVector<double> fv(100, data+100);

  // overlay
  FlatVector<double> fw(200, data);

  fu = 1.0;
  fv = fu;

  cout << "(u,v) = " << InnerProduct (fu, fv) << endl;
  cout << "(w,w) = " << InnerProduct (fw, fw) << endl;



  // access row/col of matrix:
  cout << "Row 2 of m is " << endl << m.Row(2) << endl;
  cout << "type is " << typeid(m.Row(2)).name() << endl;

  cout << "Col 1 of m is " << endl << m.Col(1) << endl;
  cout << "type is " << typeid(m.Col(1)).name() << endl;


  // more complicated vectors

  Vector<Vec<2,Complex> > sysu(4);
  Vector<Vec<3,Complex> > sysv(3);
  Matrix<Mat<3,2,Complex> > sysm(3,4);

  for (int i = 0; i < sysu.Size(); i++)
    {
      sysu(i)(0) = Complex (1,0);
      sysu(i)(1) = Complex (0,1);
    }

  cout << "sysu = " << sysu << endl;

  for (int i = 0; i < sysm.Height(); i++)
    for (int j = 0; j < sysm.Width(); j++)
      {
	sysm(i,j) = 0.0;
	sysm(i,j)(1,0) = 2.5;
      }

  cout << "sysm = " << sysm << endl;

  sysv = sysm * sysu;

  cout << "sysv = " << sysv << endl;

  Vector<double> a(1000);
  Vector<double> b(1000);
  a = 1.0;
  b = 1.0;

  double sum = 0;
  for (int i = 0; i < 1000000; i++)
    sum += InnerProduct (a, b);

  cout << "sum = " << sum << endl;

  // Vector<Complex> ac(1000);
  // Vector<Complex> bc(1000);
  Vector<complex<double> > ac(1000);
  Vector<complex<double> > bc(1000);
  // Vector<Vec<2> > ac(1000);
  // Vector<Vec<2> > bc(1000);
  ac = 1.0;
  bc = 1.0;

  Complex sumc = 0;
  for (int i = 0; i < 100000; i++)
    sumc += InnerProduct (ac, bc);

  cout << "sum = " << sumc << endl;

  return 0;
}

