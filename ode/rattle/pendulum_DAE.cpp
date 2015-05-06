// compile with
// g++ rattle.cpp pendulum_DAE.cpp ../bla/calcinverse.cpp ../bla/exception.cpp ../bla/localheap.cpp -std=c++11 -o test -lglut -lGLU -lGL
//

#include "../bla/bla.hpp"
using namespace ngbla;

#include <unistd.h>
#include "rattle.hpp"

#include <GL/glut.h>



/* ******************* The ODEs I want to solve  ******************** */

class Pendulum : public ConstrainedHamiltonian
{
  double m;
public:
  Pendulum () { m = 5; }

  virtual int Dim () const { return 2; }
  virtual int DimConstraints () const { return 1; }


  virtual double Eval (double t, const Vector<> & q, const Vector<> & p) const 
  {
    return 0.5 * m * ( p(0) * p(0) + p(1) * p(1) ) + q(1) + 1;
  }
  virtual void EvalDq (double t, const Vector<> & q, const Vector<> & p, Vector<> & Dq) const
  {
    Dq(0) = 0;
    Dq(1) = 1;
  }
  virtual void EvalDp (double t, const Vector<> & q, const Vector<> & p, Vector<> & Dp) const 
  {
    Dp(0) = m * p(0);
    Dp(1) = m * p(1);
  }
  virtual void EvalDpDp (double t, const Vector<> & q, const Vector<> & p, Matrix<> & DpDp) const 
  {
    DpDp(0,0) = m;
    DpDp(0,1) = 0;
    DpDp(1,0) = 0;
    DpDp(1,1) = m;
  }

  virtual void Constraint (double t, const Vector<> & q, const Vector<> & g) const 
  {
    g(0) = q(0) * q(0) + q(1) * q(1) - 1;
  }
  virtual void ConstraintDq (double t, const Vector<> & q, Matrix<> & Dq) const
  {
    Dq(0,0) = 2 * q(0);
    Dq(0,1) = 2 * q(1);
  }



  virtual void Redraw (const Vector<> & q) const
  {   
    glClear(GL_COLOR_BUFFER_BIT);		/* clear the display */

    double cx = q(0);
    double cy = q(1);

    DrawBar (0, 0, cx, cy);
    DrawMass (cx, cy, 0.1);
  }

  void DrawBar (double ax, double ay, double bx, double by) const
  {
    glColor3f(1.0, 0, 0);	
    glLineWidth(3);

    glBegin(GL_LINES);	
    glVertex2d(ax, ay);	
    glVertex2d(bx, by);
    glEnd();
  }

  void DrawMass (double cx, double cy, double r) const
  {
    glColor3f(0.0, 1.0, 0.0);	
    glBegin(GL_POLYGON);
    for (int j = 0; j < 20; j++)
      {
        glVertex2d (cx+r*sin(2*M_PI*j/10), 
                    cy+r*cos(2*M_PI*j/10));
      }
    glEnd();
  }
};











ODESolver * psolver;


void displayCB(void)		/* function called whenever redisplay needed */
{
  psolver -> Redraw();
  glutSwapBuffers();
}

bool working = 0;


void keyCB(unsigned char key, int x, int y)	/* called on key press */
{
  if( key == 'q' ) exit(0);
  if( key == 's' ) 
    {
      psolver -> Step();
      displayCB();
    }
      
  if( key == 'b' ) working = 1;
  if( key == 'e' ) working = 0;
}

void idlefunction ()         /* function called when nothing else to do */
{
  if (working)
    {
      for (int i = 0; i < 10; i++)
	psolver -> Step();
      usleep(2000);
      
      glutPostRedisplay();
      // displayCB();
    }
}





int main (int argc, char ** argv)
{
  Pendulum pendulum;
  Rattle ssm;

  Vector<> q0(2), p0(2);
  q0(0) = 1;
  q0(1) = 0;
  p0(0) = 0;
  p0(1) = 0;

  ODESolver solver(pendulum, ssm, 0, q0, p0, 100000, 1e-3);
  // solver.Solve();


  psolver = &solver;


  cout << "Press  'q'uit   'b'egin  'e'nd  's'tep" << endl;

  int win;

  glutInit(&argc, argv);		/* initialize GLUT system */

  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowSize(500,500);		/* width=400pixels height=500pixels */
  win = glutCreateWindow("Hamiltonian");	/* create window */

  /* from this point on the current window is win */
  glClearColor(0.0,0.0,0.0,0.0);	/* set background to black */
  gluOrtho2D(-2, 2, -2, 2);		/* how object is mapped to window */
  glutDisplayFunc(displayCB);		/* set window's display callback */
  glutKeyboardFunc(keyCB);		/* set window's key callback */

  glutIdleFunc (idlefunction);

  glutMainLoop();			/* start processing events... */




  return 0;
}
