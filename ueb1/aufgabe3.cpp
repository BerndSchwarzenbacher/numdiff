#include <fstream>
#include <math.h>
#include <iostream>

int main ()
{
  std::ofstream values_file;
  values_file.open ("values");

  // Optionen
  int max_time = 1e4;
  int steps_per_second = 1e5;
  double step = 1.0/steps_per_second;

  // Startwerte
  double alpha_old = M_PI/4;
  double alpha_dot_old = 0;

  double alpha_new = 0;
  double alpha_dot_new = 0;

  for (int i = 0; i < max_time * steps_per_second; ++i)
  {
    alpha_new = alpha_old + step * alpha_dot_old;
    alpha_dot_new = alpha_dot_old - step * sin(alpha_old);

    // nur jede Sekunde speichern
    if (i % steps_per_second == 0)
      values_file << alpha_new << " " << alpha_dot_new << std::endl;

    alpha_old = alpha_new;
    alpha_dot_old = alpha_dot_new;
  }

  values_file.close();
  return 0;
}

