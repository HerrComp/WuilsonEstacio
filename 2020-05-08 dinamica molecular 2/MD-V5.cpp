//con Algoritmo de leap-frog
#include <iostream>
#include <fstream>
#include <vector>

/* este programa simula un cuerpo que cae vajo la accion de la gravedad y en el futuro va a rebotar contra el suelo y contra otros cuerpos */

// cuerpo
struct body {
  double mass;
  double r[3], v[3], f[3];
};

// simulation conditions o condiciones externas
const int N = 1;  // 1 porque es una sola particula
const double G = 9.81;
const double DT = 0.1; // derivada de tiempo

// helper function o funcion auxiliar
void initial_conditions(std::vector<body> & bodies); // crea las condiciones iniciales y recibe los arreglos
void timestep(std::vector<body> & bodies, double dt); // para cambio de posicion old a new
void start_timeintegration(std::vector<body> & bodies, double dt); //recibe los tiempo
void compute_force(std::vector<body> & bodies); // para las fuerzas
void print_system(const std::vector<body> & bodies, double time); // recibe los cuerpos y el tiempo y los imprime


int main(void){
// declaracion de los cuerpos
  std::vector<body> bodies(N); 

// pre-processing
  initial_conditions(bodies); // inicia las condiciones en los bodies r(t=0), v(t=0)
  compute_force(bodies);  //para compilar las fuerzas de los cuerpos f(t=0)
 // print_system(bodies, 0); // esto para imprimir el sistema de cuerpos
  start_timeintegration(bodies, DT);

 // std::cout << 0 << bodies[0].r[2] << std::endl; // imprime la pocicion inicial del cuerpo.

// processing
  for (int step = 0; step < 1000; ++step){
    double tstep = step*DT;
    std::cout << tstep << "  " << bodies[0].r[2] << std::endl;
    timestep(bodies, 0); // se actualiza r y v  -> actualizar f
    compute_force(bodies); // actualizar la fuerza para las nuevas pociciones 
  }

  print_system(bodies, 0); 

  return 0;
}

// helper function
void initial_conditions(std::vector<body> & bodies)
{
 
  bodies[0].mass = 1.23;
  bodies[0].r[2] = 7.86; // Altura inicial
  bodies[0].v[2] = 1.32; // velocidad inicial hacia arriba
}

void timestep(std::vector<body> & bodies, double dt)
{
  for (auto & cuerpo : bodies){
    for(int ii =0; ii <= 2; ++ii){
      cuerpo.v[ii] += dt*cuerpo.f[ii]/cuerpo.mass;
      cuerpo.r[ii] += cuerpo.v[ii]*dt;
    }
  }
}

void start_timeintegration(std::vector<body> & bodies, double dt)
{
  for (auto & cuerpo : bodies){
    for(int ii = 0; ii <=2; ++ii){
      cuerpo.v[ii]=cuerpo.v[ii] - dt*cuerpo.f[ii]/(2*cuerpo.mass);
    }
  }
}

void compute_force(std::vector<body> & bodies)
{
  // reset forces para no utilizar fuerzas de pasos anteriores
  for (auto & cuerpo : bodies) {
    cuerpo.f[0] = cuerpo.f[1] = cuerpo.f[2] = 0.0;
  }

  for (auto & cuerpo : bodies) { // add gravity
    cuerpo.f[2] -= cuerpo.mass*G;
  }
  
}

void print_system(const std::vector<body> & bodies, double time)
{
  std::ofstream fout("datos.txt", std::ofstream::out);
  fout.precision(15); fout.setf(std::ios::scientific);

// se crea un for automatico para recorrer automaticamente los arreglos y el & es para queno cree nuevas referencias o variables
  for (const auto & cuerpo : bodies) {
    fout << cuerpo.r[0] << "  " << cuerpo.r[1] << "  " << cuerpo.r[2] << "  "  
         << cuerpo.v[0] << "  " << cuerpo.v[1] << "  " << cuerpo.v[2] << "  "
         << cuerpo.f[0] << "  " << cuerpo.f[1] << "  " << cuerpo.f[2] << "  "  
         << cuerpo.mass << "\n";
  }
}

