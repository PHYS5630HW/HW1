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
// y[1] = velocity along i axis   ; f_vi = dvi/dt => acceleration along i axis  x
// y[2] = position along j axis   ; f_rj = drj/dt => velocity along j axis  y
// y[3] = velocity along j axis   ; f_vj = dvj/dt => acceleration along j axis y

const double g=9.81;    // [m/s^2]
const double m=0.145;     // [kg]  n.b. simple projectile motion does not depent on the mass
const double d = 0.075; //Diameter of ball
const double b=1.6/10000*d; // constant for air resistance
const double c=0.25*d*d; // constant for air resistance


double f_ri(double x, const vector<double> &y){  // change in position along i axis
  (void) x;   // prevent unused variable warning
  return y[1];
}
double f_vi(double x, const vector<double> &y){  // change in velocity along i axis
  (void) x;
  return -(b + c*sqrt(y[1]*y[1] + y[3]*y[3]))*y[1] / m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}
double f_rj(double x, const vector<double> &y){  // change in position along j axis
  (void) x;   // prevent unused variable warning
  return y[3];
}
double f_vj(double x, const vector<double> &y){  // change in velocity along j axis
  (void) x;
  return -(b + c*sqrt(y[1]*y[1] + y[3]*y[3]))*y[3] / m - g;
  // return g;    // if no air constant acceleration along -j direction: F/m = -g
}

double f_stop(double x, const vector<double> &y){
  (void) x;
  
  if (y[0]>18.5){ 
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
  vector<pfunc_t> v_fun(4);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  
  
  double v_0 = 40;
  for(int i = 0; i<1000; i++){
  v_0 = v_0 + .1; 
  vector<double> y0(4);
  // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
  double theta = 1*3.14159/180;
  y0[0]=0;   // init position on i-axis
  y0[1]=v_0*cos(theta);  // init velocity along i axis
  y0[2]=1.4;   // repeat for j-axis
  y0[3]=v_0*sin(theta);
  
  auto tgN = RK4SolveN(v_fun, y0, 500, 0, 3, f_stop);
  if(y0[0]>18.45 && (y0[2] > 0.88 && y0[2] < 0.92)) {//check if final condition matches
    cout<<"Final Height = "<<y0[2]<<" m \n";
    cout<<"Time = "<<TMath::MaxElement(tgN[0].GetN(),tgN[0].GetX())<<" s \n";
    cout<<"Initial Velocity = "<<v_0<<" m/s or "<<v_0*3600/1609.34<<" mph"<<endl;
    break;
  }
  //TCanvas *c2 = new TCanvas("c2","ODE solutions 2",dw,dh);
  //tgN[0].Draw("a*");
  //c2->Draw();
  }
  
  /*
  // save our graphs
  TFile *tf=new TFile("RKnDemo.root","recreate");
  for (unsigned i=0; i<v_fun.size(); i++){
    tgN[i].Write();
  }
  tf->Close();
  */
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(1,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

