#include <iostream>
#include <fstream>
#include <vector>

/* este programa simula un cuerpo que cae vajo la accion de la gravedad y en el futuro va a rebotar contra el suelo y contra otros cuerpos */


// cuerpo
struct body {
  double mass;
  double r[3], v[3], f[3];
};

const int N = 1;  
const double G = 9.81;
const double DT = 0.1; 

void initial_conditions(std::vector<body> & bodies); 
void timestep(std::vector<body> & bodies, double dt); 
void start_timeintegration(std::vector<body> & bodies, double dt); 
void compute_force(std::vector<body> & bodies); 
void print_system(const std::vector<body> & bodies, double time); 

int main(void){

  std::vector<body> bodies(N); 

  initial_conditions(bodies); 
  compute_force(bodies); 
  print_system(bodies, 0); 

  return 0;

}

// helper function
void initial_conditions(std::vector<body> & bodies)
{
  /*
    z
    |
    |____ y
   /
  x  
  */
  
  bodies[0].mass = 1.23;
  bodies[0].r[2] = 7.86; 
  bodies[0].v[2] = 1.32; 
}

void timestep(std::vector<body> & bodies, double dt);
void start_timeintegration(std::vector<body> & bodies, double dt);

void compute_force(std::vector<body> & bodies)
{
  
  for (auto & cuerpo : bodies) {
    cuerpo.f[0] = cuerpo.f[1] = cuerpo.f[2] = 0.0;
  }

  for (auto & cuerpo : bodies) {
    // add gravity
    cuerpo.f[2] -= cuerpo.mass*G;
  }
  
}

void print_system(const std::vector<body> & bodies, double time)
{
  std::ofstream fout("datos.txt", std::ofstream::out);
  fout.precision(15); fout.setf(std::ios::scientific);


  for (const auto & cuerpo : bodies) {
    fout << cuerpo.r[0] << "  " << cuerpo.r[1] << "  " << cuerpo.r[2] << "  "  
         << cuerpo.v[0] << "  " << cuerpo.v[1] << "  " << cuerpo.v[2] << "  "
         << cuerpo.f[0] << "  " << cuerpo.f[1] << "  " << cuerpo.f[2] << "  "  
         << cuerpo.mass << "\n";
  }
}

