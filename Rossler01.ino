// Rossler01.ino
// Tangible Networks
// Espen Knoop 25th Feb 2015
// 
// Implements a Rossler system, that can be coupled to neighbouring nodes 
// to demonstrate synchronisation, consensus and control.
// Uses a Runge-Kutta ODE solver.


#include <TN05.h>

static double range = 15.0;


TN Tn = TN(-range,range);

// Parameters
double a = 0.2;
double b = 0.2;
double lambda;
double lambdaRange[] = {1.0, 5.0};
double kRange[] = {0.0, 1.0};
double k;

// For RK
double k1x,k2x,k3x,k4x;
double k1y,k2y,k3y,k4y;
double k1z,k2z,k3z,k4z;

boolean in_con[] = {0,0,0,0,0,0};
double in_val[] = {0.0,0.0,0.0,0.0,0.0,0.0};

double x, y, z;
double dt = 0.005; // for eqn
//int timestep = 2; // ms
boolean sw = 0, msw = 0;
int beepCtr = 0;

void setup() {
  Serial.begin(115200);
  reset(); // go to random initial condition
}

void loop() {
  if ( (!Tn.sw() && sw) || (!Tn.masterSw() && msw) ) reset();
  sw = Tn.sw();
  msw = Tn.masterSw();
  lambda = lambdaRange[0] + Tn.pot()*(lambdaRange[1] - lambdaRange[0]);
  k = kRange[0] + Tn.masterRead()*(kRange[1] - kRange[0]);

  // Do 10 timesteps per colour/switch update
  for (int i=0; i<10; i++) {
    for (int j=0; j<6; j++) {
      in_con[j] = Tn.isConnected(j);
      in_val[j] = Tn.analogRead(j);
    }
    // Runge-kutta
    k1x = dxdt(x,y,z);
    k1y = dydt(x,y,z);
    k1z = dzdt(x,y,z);
    k2x = dxdt(x + 0.5*dt*k1x,y + 0.5*dt*k1y,z + 0.5*dt*k1z);
    k2y = dydt(x + 0.5*dt*k1x,y + 0.5*dt*k1y,z + 0.5*dt*k1z);
    k2z = dzdt(x + 0.5*dt*k1x,y + 0.5*dt*k1y,z + 0.5*dt*k1z);
    k3x = dxdt(x + 0.5*dt*k2x,y + 0.5*dt*k2y,z + 0.5*dt*k2z);
    k3y = dydt(x + 0.5*dt*k2x,y + 0.5*dt*k2y,z + 0.5*dt*k2z);
    k3z = dzdt(x + 0.5*dt*k2x,y + 0.5*dt*k2y,z + 0.5*dt*k2z);
    k4x = dxdt(x + dt*k3x,y + dt*k3y,z + dt*k3z);
    k4y = dydt(x + dt*k3x,y + dt*k3y,z + dt*k3z);
    k4z = dzdt(x + dt*k3x,y + dt*k3y,z + dt*k3z);
    
    x += dt*(k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0;
    y += dt*(k1y + 2.0*k2y + 2.0*k3y + k4y)/6.0;
    z += dt*(k1z + 2.0*k2z + 2.0*k3z + k4z)/6.0;
    
    // periodic boundary conditions in x and y
    while (x > range) x -= range*2;
    while (x < -range) x += range*2;
    while (y > range) y -= range*2;
    while (y < -range) y += range*2;
    
    Tn.analogWrite(2*x);
    //delay(timestep);  // To slow down simulation
  }

  // Beep 
  // If we're below thresh, reset
  // If we're reset, and hit thresh: start beep.
  // If we're beeping and hit count, stop beep.
  // Like poincare section!
  if (x < 0) {
    beepCtr = 0;
  }
  if (x > 0 && beepCtr == 0) {
    beepCtr = 1;
    analogWrite(SPKR,128);
  }
  if (beepCtr > 0 && beepCtr < 1000) beepCtr++;
  if (beepCtr == 15) analogWrite(SPKR,0);
  

  // Set the output colour  
  // Tn.colour(0.1*(x+5.0), 0.1*(y+5.0), 0.2*z);
  Tn.colour(0.1*(x+5.0),0.0,0.1*(5.0-x));
  // Serial.println(x);
  // Tn.printState();
}



// Pick a random starting point close to origin
void reset() {
  x = double(random(100))*0.01*3.0-1.5; // -5 - 5
  y = double(random(100))*0.01*3.0-1.5; // -5 - 5
  z = double(random(100))*0.01*5; // -5 - 5
}

double dxdt(double xin, double yin,double zin) {
  double out;
  out = -yin-zin;
  for (int j=0; j<6; j++) { 
    if (in_con[j]) out += k*(in_val[j] - xin);
  }
  return out;

}

double dydt(double xin,double yin, double zin){
  return xin + a*yin;
}

double dzdt(double xin, double yin, double zin){
  return b + zin*(xin-lambda);
}
  



