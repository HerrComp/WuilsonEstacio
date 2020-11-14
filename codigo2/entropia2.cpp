#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random64>
using namespace std;

const double K=1e4;
const double Lx=1000,Ly=1000;
const int Nx=7, Ny=7, N=Nx*Ny;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

class Cuerpo;
class Colisionador;

//-----------Clase Cuerpo-----------
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double GetV(void){return norma(V);};
  double GetVx(void){return V.x();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}
void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::Mueva_r(double dt,double Constante){
  r+=V*(Constante*dt);
}
void Cuerpo::Mueva_V(double dt,double Constante){
  V+=F*(Constante*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//-----------Clase Colisionador----------
class Colisionador{
private:

public:
  void CalculeTodasLasFuerzas(Cuerpo* Grano);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2);
};
void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano){
  int i,j;
  //Borrar todas las fuerzas
  for(i=0;i<N+4;i++) Grano[i].BorreFuerza();
  //Calcular todas las fuerzas entre parejas de granos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeLaFuerzaEntre(Grano[i],Grano[j]);
}
void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2){
  vector3D F2, dr,r_unitario; double distancia,s;
  dr=Grano2.r-Grano1.r; distancia=norma(dr); r_unitario=dr/distancia;
  s=(Grano1.R+Grano2.R)-distancia;
  if(s>0){ //Si se chocan,
    F2=r_unitario*(K*pow(s,1.5));  Grano2.AgregueFuerza(F2);  Grano1.AgregueFuerza(F2*(-1));
  }
}
//-----------Funciones Globales---------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Gas.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:1010]"<<endl;
  cout<<"set yrange [-10:1010]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<1000/7<<"*t,0";
    cout<<" , "<<1000/7<<"*t,1000";
    cout<<" , 0,"<<1000/7<<"*t";
    cout<<" , 1000,"<<1000/7<<"*t";
}
void TermineCuadro(void){
    cout<<endl;
}

//-------Programa Principal------------
int main(void){
  double t, dt=1e-2;
  double tdibujo; int Ndibujos;
  Cuerpo Grano[N+4]; int i,j;
  Colisionador Newton;
  Crandom ran64(1);
  double theta;
  double m0=1,R0=1,V=1;
  double Rpared=10000,Mpared=10000; 
  double T=Lx/V, teq = 10*T, tmax=T;
  double dx = Lx/(Nx+1), dy = Ly/(Ny+1);
  
  InicieAnimacion(); Ndibujos=2000;
  
 
  //PARTICULAS
  

  for(i=0;i<N;i++){
      //--------------------(x0  ,y0  ,z0,  Vx0,  Vy0,Vz0,m0,R0);
      theta=2*M_PI*ran64.r();
      Grano[i].Inicie(80*cos(i*Phi),80*sin(i*Phi),0, V*cos(theta),V*sin(theta),0, m0,R0);
    }
  
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    //cout<<Grano[1].Getx()<<" "<<Grano[1].Gety()<<endl;
    //Muevase con Omelyan PEFRL
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Zeta);
    Newton.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,(1-2*Lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,1-2*(Xi+Zeta));
    Newton.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,(1-2*Lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Zeta);
  }  
  return 0;
}