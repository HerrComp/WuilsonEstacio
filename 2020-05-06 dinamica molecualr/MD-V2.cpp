#include <iostream>
#include <fstream>
#include <vector>

/* este programa simula un cuerpo que cae vajo la accion de la gravedad y en el futuro va a rebotar contra el suelo y contra otros cuerpos */

/*  cuerpo:
-masas(densidad)
- (forma: esfera)
- r[3],v[3],F[3]

condiciones externas
-gravedad
*/

/*
 Funciones
  - initial_conditions

  - timestep (nueva posicion y nueva velocidad)
  - start_timeintegration()
  - compute_force()
    - (implementar la fuerza de Hertz)
    - (Fuerza del rebote)
    - *(Fuerza de gravedad)
    - (Amortiguamiento)
    - (Fuerza de fricción, si es que se considera)
  - print_system
  
  - Visualization
*/

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

  std::vector<body> bodies(N); 

  initial_conditions(bodies); // inicia las condiciones en losbopdies
  compute_force(bodies);  //para compilar las fuerzas de los cuerpos
  print_system(bodies, 0); // esto para imprimir el sistema de cuerpos

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
  bodies[0].r[2] = 7.86; // Altura inicial
  bodies[0].v[2] = 1.32; // velocidad inicial hacia arriba
}

void timestep(std::vector<body> & bodies, double dt);
void start_timeintegration(std::vector<body> & bodies, double dt);

void compute_force(std::vector<body> & bodies)
{
  // reset forces para no utilizar fuerzas de pasos anteriores
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

// se crea un for automatico para recorrer automaticamente los arreglos y el & es para queno cree nuevas referencias o variables
  for (const auto & cuerpo : bodies) {
    fout << cuerpo.r[0] << "  " << cuerpo.r[1] << "  " << cuerpo.r[2] << "  "  
         << cuerpo.v[0] << "  " << cuerpo.v[1] << "  " << cuerpo.v[2] << "  "
         << cuerpo.f[0] << "  " << cuerpo.f[1] << "  " << cuerpo.f[2] << "  "  
         << cuerpo.mass << "\n";
  }
}

