#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;


 
// functions to describe simple projectile motion
// here use use ri,rj,rk to define directions to prevent confusion with
// standard ODE notation, where x=independent variable, \vec y=dependent variable(s)

// first we define our vaiables
// x    = time
// y[0] = position along i axis   ; f_ri = dri/dt => velocity along i axis x 
// y[1] = velocity along i axis   ; f_vi = dvi/dt => acceleration along i axis x
// y[2] = position along j axis   ; f_rj = drj/dt => velocity along j axis  y
// y[3] = velocity along j axis   ; f_vj = dvj/dt => acceleration along j axis y
// y[4] z
// y[5] Vz

const double PI = 3.1415926535897;

const double g=9.81*3.2808;    // [m/s^2]
const double m=0.145;     // [kg]  n.b. simple projectile motion does not depent on the mass
const double d = 0.075; //Diameter of ball
const double b=1.6/10000*d; // constant for air resistance
const double c=0.25*d*d; // constant for air resistance
const double B = 4.1*pow(10, -4); //constant B from (3.43) on page 127.

const TString name = "Fastball";   //What type of throw
double fi = 0;
//const double w = 1800.0/60; //angular velocity per SECOND.
const double w = 10000.0/60; //angular velocity per SECOND.


double F_v(double v) { //Calculate F(v) from equation(3.41) shown on page 126. specify for baseball.
  //const double vd = 35.0; //m/s
  const double vd = 35.0*3.2808; //m/s
  //const double delta = 5.0; //m/s
  const double delta = 5.0*3.2808; //m/s
  return 0.0039+0.0058/(1+exp((v-vd)/delta));
}
double vijk (double vi, double vj, double vk) { //overall velocity |v|
	return sqrt(vi*vi+vj*vj+vk*vk);
}

double f_ri(double x, const vector<double> &y){  // change in position along x axis
  (void) x;   // prevent unused variable warning
  return y[1];
}
double f_vi(double x, const vector<double> &y){  // change in velocity along x axis
  (void) x;
  double vv = vijk(y[1], y[3], y[5]);
  return -1*(F_v(vv))*vv*y[1]+B*w*(y[5]*sin(fi) - y[3]*cos(fi));
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}
double f_rj(double x, const vector<double> &y){  // change in position along y axis
  (void) x;   // prevent unused variable warning
  return y[3];
}
double f_vj(double x, const vector<double> &y){  // change in velocity along y axis
  (void) x;
  double vv = vijk(y[1], y[3], y[5]);
  return -1*(F_v(vv))*vv*y[3]+B*w*y[1]*cos(fi);
  // return g;    // if no air constant acceleration along -j direction: F/m = -g
}
double f_rk(double x, const vector<double> &y){  // change in position along z axis
  (void) x;   // prevent unused variable warning
  return y[5];
}
double f_vk(double x, const vector<double> &y){  // change in velocity along z axis
  (void) x;
  double vv = vijk(y[1], y[3], y[5]);
  return -1*(F_v(vv))*vv*y[5]-B*w*y[1]*sin(fi)-g;

}

double f_stop(double x, const vector<double> &y){
  (void) x;
  
  if (y[0]>60){ 
    return 1;  // stop calulation if the current step takes height to negative value
  }
  return 0;  // continue calculation
}


int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  //UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  //UInt_t dw = 1.1*dh;
  // ******************************************************************************
  

  // *** test 2: Use RK4SolveN to calculate simple projectile motion
  vector<pfunc_t> v_fun(6);   // 6 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  v_fun[4]=f_rk;
  v_fun[5]=f_vk;
  
  
  vector<double> y0(6);

  if(name == "Slider") fi = 0*PI/180; //angular velocity angle in radian.
  if(name == "Curveball") fi = 45*PI/180; //angular velocity angle in radian.
  if(name == "Screwball") fi = 135*PI/180; //angular velocity angle in radian.
  if(name == "Fastball") fi = 225*PI/180; //angular velocity angle in radian.
  
  double v_0 = 85.0*1.46667; //initial velocity in m/s. 
  if(name =="Fastball") v_0 = 95.0*1.46667; //initial velocity in m/s. 
  const double theta =1.0*PI/180; //initial velocity angle in radian.
  const double h = 1*pow(10, -4); //step length. 1/h is the step needed.
  
  y0[0]=0;   // init position on i-axis
  y0[1]=v_0*cos(theta);  // init velocity along i axis
  y0[2]=0;   // repeat for j-axis
  y0[3]=0;
  y0[4]=0;     //repreat for k-axis
  y0[5]=v_0*sin(theta);
  
  auto tgN = RK4SolveN(v_fun, y0, (int)(1.0/h), 0, 30, f_stop);
  /*
  TCanvas *c2 = new TCanvas("c2","ODE solutions 2",dw,dh);
  tgN[2].Draw("a*");
  c2->Draw();
  */
  
  // save our graphs
  TFile *tf=new TFile(name+".root","recreate");
  for (unsigned i=0; i<v_fun.size(); i++){
    tgN[i].Write();
  }
  tf->Close();
  
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
