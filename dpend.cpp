#include <iostream>
#include <cstdlib>
#include <cmath>
#include "cpgplot.h"
#include <algorithm> 
//#include <thread>
//#include <chrono>

// A Program for simulating a double pendulum
// Filename: dpend.cpp Src: AR 12/2020
// File based on lorenz.cpp example
// Modified to simulate a double pendulum

  const float g = 9.81; //Gravity in meters per second squared
  const float twoPI = 2.0 * M_PI;


  
//-----------------------------------------------------------------------
//  This function draws a straight line using cpgplot
//  Parameters are the x,y for start and end of the line
//  also specifies the colour of the line to draw.
//-----------------------------------------------------------------------
void DrawStraightLine(float x0, float y0, float x1, float y1, int iCol )
{
    float lx[2],ly[2];
    
    lx[0] = x0;
    lx[1] = x1;
    ly[0] = y0;
    ly[1] = y1;
    
    cpgsci(iCol);
    cpgline(2,lx,ly);
}

//
//-----------------------------------------------------------------------
//  This function calculates the derivative of the angular velocity 
//  of m1 with respect to time
//  
//-----------------------------------------------------------------------
float dw1dt(float theta1, float theta2, float w1, float w2, float r1, float r2, float m1, float m2)
{
    float dwdt;
    //Calculation broken into smaller parts to help debug equations
    float top; 
    float top1; 
    float top2; 
    float top3; 
    float bot; 
  
    top1 = (-g) * ((2.0 * m1) + m2) * sin(theta1); 
    // std::cout << "top1 " << top1 << "\n";    //used for debugging
    
    top2 = ((-m2) * g * sin(theta1 - 2.0 * theta2));
    // std::cout << "top2 " << top2 << "\n";    //used for debugging
    
    top3 = ((-2.0) * sin(theta1 - theta2) * m2 * ((pow(w2, 2.0) * r2) + (pow(w1, 2.0) * r1 * cos(theta1 - theta2))));
    // std::cout << "top3 " << top3 << "\n";    //used for debugging
  
    bot = r1 * ((2.0 * m1) + m2 - (m2 * cos((2.0 * theta1) - (2.0 * theta2))));
    
    top = top1 + top2 + top3;
    
    //  std::cout << "top " << top << "\n";    //used for debugging
    //  std::cout << "bot " << bot << "\n";  //used for debugging
    
    dwdt = top / bot;
  
    return dwdt;
     
}
    

//
//-----------------------------------------------------------------------
//  This function calculates the derivative of the angular velocity 
//  of m2 with respect to time
//  
//-----------------------------------------------------------------------

float dw2dt(float theta1, float theta2, float w1, float w2, float r1, float r2, float m1, float m2)
{
    float dwdt;
    //Calculation broken into smaller parts to help debug equations
    float top; 
    float top1; 
    float top2; 
    float top3; 
    float bot; 
    
    top1 = 2.0 * sin(theta1 - theta2);
    
    top2 = pow(w1,2) * r1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + pow(w2,2.0) * r1 * m2 * cos(theta1 - theta2);
    
    bot = r1 * (2.0 * m1 + m2 - m2 * cos(2.0 * theta1 - 2.0 * theta2));
    
    top = top1 * top2;
    
    dwdt = top / bot;
    
    return dwdt;

    // return (2.0 * sin(theta1 - theta2) * (pow(w1,2) * r1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + pow(w2,2.0) * r1 * m2 * cos(theta1 - theta2)))/(r1 * (2.0 * m1 + m2 - m2 * cos(2.0 * theta1 - 2.0 * theta2)));
}


//
//-----------------------------------------------------------------------
//  This function converts radians to degrees, taking into account that 
//  the angle is represented from 0 to 180 and 0 to -180
//  
//-----------------------------------------------------------------------
float Radians2Degrees(float radians)
{
  float angle;
  float numRot;
  
      angle = radians * 57.2957795;
      
      if (angle > 360.0)
      {
          numRot = floor(angle / 360);
          angle = angle - (360.0 * numRot);
        //   std::cout << "angle " << angle << "\n";
      } else if (angle < (-360.0))
      {
          numRot = floor(abs(angle) / 360.0);
          angle = angle + (360.0 * numRot);
        //   std::cout << "angle " << angle << "\n";
      }
      
      if (angle > 180)
      {
          angle = angle - 360.0;
      }
      if (angle < -180)
      {
          angle = angle + 360.0;
      }

  return angle;

}


int main() 
{

  //-----------------------------------------------------------------------
  // 
  // Double pendulum assumes a massless rod and a frictionless pivot
  // 
  //-----------------------------------------------------------------------

  int iStep;
  int previStep;
  int iCol = 1;
  int i;
  int nStep = 10000; 
  
  char ans;
  
  float dt = 0.01;
  float dtheta11, dtheta12, dtheta13, dtheta14, dtheta1;//Elements of the derivation for m1
  float dw1, dw2, dw3, dw4, dw;//Elements of the derivation for m1
  float dtheta21, dtheta22, dtheta23, dtheta24, dtheta2;//Elements of the derivation for m1
  float dw21, dw22, dw23, dw24, dww2;//Elements of the derivation for m2
  float m1; //Mass of first pendulum kg
  float m2; //Mass of second pendulum kg
  float r1; //Length of first pendulum rod m
  float r2; //Length of second pendulum rod m
  float theta1init; //Initial Angle 1 
  float theta2init; //Initial Angle 2
  
  // Arrays for plotting
  float x1[nStep+1], y1[nStep+1]; //Position of mass 1
  float x2[nStep+1], y2[nStep+1];
  float theta1[nStep+1]; //Angle 1 in Radians
  float theta2[nStep+1]; //Angle 2 in Raddians
  float w1[nStep+1];  //Angular velocity 1
  float w2[nStep+1];  //Angular velocity 2
  float a1[nStep+1];  //Acceleration theta 1
  float a2[nStep+1];
  float timeaxis[nStep+1];
  float theta1angle[nStep+1]; //Angle 1 in degrees
  float theta2angle[nStep+1]; //Angle 2 in degrees
  float miny;
  float maxy;
  float bobsize;  
  float plotsize;



    //Default values
    r1 = 1.0;    // meters
    r2 = 1.0; 
    m1 = 1.0;     //kg
    m2 = 1.0;
    theta1init = 90.0; //degrees
    theta2init = 45.0;
  
  
    std::cout << "\nPress D for default initial parameters.\n";
    std::cout << "Any other key to specify your own intial parameters.\n";
    std::cin >>  ans;
    
    if (toupper(ans) == 'D') {
     std::cout << "I'm picking the initial parameters standby ...\n";
     //Pick default above
    }
    else 
    {
      std::cout << "Enter mass (kg) for first pendulum must be > 0\n";
      std::cin >>  m1;
      if (m1 == 0.0) {
        std::cout << "Invalid entry try again!!!!\n";
        exit(0);
      }
      std::cout << "Enter mass (kg) for second pendulum\n";
      std::cin >>  m2;
      
      std::cout << "Enter rod length (m) for first pendulum must be > 0\n";
      std::cin >>  r1;
      if (r1 == 0.0) {
        std::cout << "Invalid entry try again!!!!\n";
        exit(0);
      }
      std::cout << "Enter rod length (m) for second pendulum must be > 0\n";
      std::cin >>  r2;
      if (r2 == 0.0) {
        std::cout << "Invalid entry try again!!!!\n";
        exit(0);
      }
      std::cout << "Enter the angle (dgrees) of the first pendulum rod\n";
      std::cin >>  theta1init;
      std::cout << "Enter the angle (degrees) of the second pendulum rod\n";
      std::cin >>  theta2init;   
    }
    
    if (cpgopen("/XWINDOW") < 0) exit(1);
   
    iStep = 0;
 
    //Initial conditions
    theta1[0] = (theta1init * M_PI) / 180;
    theta2[0] = (theta2init * M_PI) / 180;
    
    theta1angle[0] = theta1init;
    theta2angle[0] = theta2init;
    
    
    x1[0] = r1 * sin(theta1[0]); 
    y1[0] = (-r1) * cos(theta1[0]);
   
    x2[0] = x1[0] + (r2 * sin(theta2[0])); 
    y2[0] = y1[0] - (r2 * cos(theta2[0]));
    
    w1[0] = 0.0;
    w2[0] = 0.0;
    timeaxis[0] = 0.;
    
    for (iStep = 1; iStep <= nStep;iStep++) { 
    
      previStep = iStep - 1;
      timeaxis[iStep] = iStep;
      
      //Runge-Kutta
      
      dw1 = dt * dw1dt(theta1[previStep], theta2[previStep], w1[previStep], w2[previStep], r1, r2, m1, m2);
      dtheta11 = dt * w1[previStep];
      
      //
      dw21 = dt * dw2dt(theta1[previStep], theta2[previStep], w1[previStep], w2[previStep], r1, r2, m1, m2);
      dtheta21 = dt * w2[previStep];
      
      dw2 = dt * dw1dt(theta1[previStep] + (dtheta11/2), theta2[previStep]+ (dtheta21/2), w1[previStep] + (dw1/2), w2[previStep] + (dw21/2), r1, r2, m1, m2);
      dtheta12 = dt * (w1[previStep] + dw1/ 2);
      
      //
      dw22 = dt * dw2dt(theta1[previStep] + (dtheta11/2), theta2[previStep] + (dtheta21/2), w1[previStep] + (dw1/2), w2[previStep] + (dw21/2), r1, r2, m1, m2);
      dtheta22 = dt * (w2[previStep] + dw21/ 2);
      
      
      dw3 = dt * dw1dt(theta1[previStep] + (dtheta12/2), theta2[previStep]+ (dtheta22/2), w1[previStep] + (dw2/2), w2[previStep]+ (dw22/2), r1, r2, m1, m2);
      dtheta13 = dt * (w1[previStep] + dw2/ 2);
      
      //
      dw23 = dt * dw2dt(theta1[previStep]+ (dtheta12/2), theta2[previStep] + (dtheta22/2), w1[previStep]+ (dw2/2), w2[previStep] + (dw22/2), r1, r2, m1, m2);
      dtheta23 = dt * (w2[previStep] + dw22/ 2);
      
      dw4 = dt * dw1dt(theta1[previStep] + (dtheta13), theta2[previStep]+ (dtheta23), w1[previStep] + (dw3), w2[previStep]+ (dw23), r1, r2, m1, m2);
      dtheta14 = dt * (w1[previStep] + dw3);
      
      //
       dw24 = dt * dw2dt(theta1[previStep] + (dtheta13), theta2[previStep] + (dtheta23), w1[previStep]+ (dw3), w2[previStep] + (dw23), r1, r2, m1, m2);
      dtheta24 = dt * (w2[previStep] + dw23);    
      
      dw = (dw1 + 2 * dw2 + 2 * dw3 + dw4) / 6;
      dtheta1 = (dtheta11 + 2 * dtheta12 + 2 * dtheta13 + dtheta14) / 6;
       
      //
      dww2 = (dw21 + 2 * dw22 + 2* dw23 + dw24) / 6;
      dtheta2 = (dtheta21 + 2 * dtheta22 + 2 * dtheta23 + dtheta24) / 6;  
                    
      theta1[iStep] = theta1[previStep] + dtheta1;
      w1[iStep] = w1[previStep] + dw  ;
      
      //
      theta2[iStep] = theta2[previStep] + dtheta2;
      w2[iStep] = w2[previStep] + dww2 ;
             
      x1[iStep] = r1 * sin(theta1[iStep]);
      y1[iStep] = (-r1) * cos(theta1[iStep]);
   
      x2[iStep] = x1[iStep] + (r2 * sin(theta2[iStep])); 
      y2[iStep] = y1[iStep] - (r2 * cos(theta2[iStep]));
      
      
      theta1angle[iStep] = Radians2Degrees(theta1[iStep]);
      theta2angle[iStep] = Radians2Degrees(theta2[iStep]);
                    
     } //END FOR
    
    
    plotsize = (r1 + r2) + ((r1 + r2)*0.1);
    
    //bobsize i.e circle size needs to be relative to the plot size
    bobsize = ((r1 + r2)*0.1)/2;
    // Set-up plot axes;
 
    cpgenv(-plotsize,plotsize,-plotsize,plotsize,1,1);
    cpglab("X", "Y", "DOUBLE PENDULUM");
  
    //INITIAL CONDITION
    //Draw Origin   
    cpgsci(1);
    cpgcirc(0,0,bobsize/2);  
    
    //Draw line between origin and m1
    cpgsci(4);
    DrawStraightLine(0.,0.,x1[0],y1[0],2);

    //Draw initial position of mass1
    //Draw the circle relative the co-ordinate plotsize
    cpgsci(2);
    cpgcirc(x1[0],y1[0],bobsize);
    
    //Draw initial position of mass2
    cpgsci(3);
    cpgcirc(x2[0],y2[0],bobsize);
    
    //cpgline(2,lx,ly) ie line between m1 and m2
    cpgsci(2);
    DrawStraightLine(x1[0],y1[0],x2[0],y2[0],1);

    //END CONDITION
    //Draw line between origin and m1
    cpgsci(5);
    DrawStraightLine(0.,0.,x1[nStep],y1[nStep],5);

    //Draw initial position of mass1
    cpgsci(6);
    cpgcirc(x1[nStep],y1[nStep],bobsize);
    
    //Draw initial position of mass2
    cpgsci(7);
    cpgcirc(x2[nStep],y2[nStep],bobsize);
    
    //cpgline(2,lx,ly);
    cpgsci(8);
    DrawStraightLine(x1[nStep],y1[nStep],x2[nStep],y2[nStep],8);

   // cpgeras();
     
    // Plot the trajectory of mass 1 and mass 2
    cpgsci(1); 
    cpgline(nStep,x1,y1);
    
    cpgsci(2); 
    cpgline(nStep,x2,y2);
 
   
    //PLOT the angle in radians
 
  /*  
    cpgsci(1);
    cpgenv(0,nStep,-10,10,0,1);
    cpgsci(1);
    cpglab("X", "Y", "DOUBLE PENDULUM PLOT RADIANS: WHITE MASS1 GREEN MASS2");
    
    cpgsci(1); 
    cpgline(nStep,timeaxis,theta1);
    
    cpgsci(3); 
    cpgline(nStep,timeaxis,theta2);
  */
    
    //PLOT the angle in degrees
    
    miny = *std::min_element(theta2angle, theta2angle+nStep) - 10.0;
    maxy = *std::max_element(theta2angle, theta2angle+nStep) + 10.0;
    cpgsci(1);
    cpgenv(0,nStep,miny,maxy,0,1);
    //cpgsci(1);
    cpglab("X", "Y", "DOUBLE PENDULUM PLOT ANGLE: WHITE MASS1 GREEN MASS2");
    
    cpgsci(1); 
    cpgline(nStep,timeaxis,theta1angle);
    
    cpgsci(3); 
    cpgline(nStep,timeaxis,theta2angle);
    
   std::cout << "\nDo you want a plot the angles in steps? Y or N\n";
   std::cin >>  ans;

    if (toupper(ans) == 'Y') {
      for (i = 4; i > 0;i-=2) { //used to zoom in into the initial values
        cpgsci(1);  
        cpgenv(0,nStep/i,miny,maxy,0,1);
        //cpgsci(1);
        cpglab("X", "Y", "DOUBLE PENDULUM PLOT ANGLE: WHITE MASS1 GREEN MASS2");
        
        cpgsci(1); 
        cpgline(nStep/i,timeaxis,theta1angle);
        
        cpgsci(3); 
        cpgline(nStep/i,timeaxis,theta2angle);
      }
    
    }
  
  /*  
  std::cout << "\nDo you want a plot of the angular velocity? Y or N\n";
  std::cin >>  ans;
  
  if (toupper(ans) == 'Y') {
  
      miny = *std::min_element(w2, w2+nStep);
      maxy = *std::max_element(w2, w2+nStep);
      cpgsci(1); 
      cpgenv(0,nStep/4,miny,maxy,0,1);
      //cpgsci(1);
      cpglab("X", "Y", "DOUBLE PENDULUM PLOT ANGULAR VELOCITY: WHITE MASS1 GREEN MASS2");
      
      cpgsci(1); 
      cpgline(nStep/4,timeaxis,w1);
      
      cpgsci(3); 
      cpgline(nStep/4,timeaxis,w2);
  
  }
  */
  
    
  // Pause and then close plot window
  cpgclos();
}
