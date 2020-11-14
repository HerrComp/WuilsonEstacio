#include "md-simul.h"

int main(void)
{
  Particle ball;

  // start system
  initial_conditions(ball);
  compute_force(ball);
  start_integration(ball, DT);
  print(ball, 0.0);

  // evolve
  for(int istep = 0; istep < NSTEPS; ++istep) {
    std::string fname ="datos-" + std::to_string(istep)+ ".csv"
    std::ofstream fout("datos-"+istep+".csv");
    time_integration(ball, DT);
    compute_force(ball);
   // print(ball, istep*DT);
   fout <<ball.Rx<< ", "
   << ball.Ry << ", "
   << ball.Rz << ", "
   << ball.mass << ", "
   <<  ball.rad << ", "
   fout.close();
  }