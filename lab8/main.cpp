// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Lab: Lab 8
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <GL/glut.h>
#include "color.h"
#include "vector3.h"

using namespace std;


const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const float VIEW_LEFT = 0.0;
const float VIEW_RIGHT = WINDOW_WIDTH;
const float VIEW_BOTTOM = 0.0;
const float VIEW_TOP = WINDOW_HEIGHT;
const float VIEW_NEAR = -400;
const float VIEW_FAR = 400;


#define MAX_ITERS 30

void renderJuliaSet(double xmin, double xmax, double ymin, double ymax );
void GLrender();
void GLInit(int* argc, char** argv);
void renderPixel(int x, int y, Color3d& color, float sz = 1.0);

int main(int argc, char** argv)
{
  GLInit(&argc, argv);
  glutMainLoop();
  
  return 0;
}

void renderJuliaSet(double xmin, double xmax, double ymin, double ymax )
{
  double dx = (xmax-xmin)/WINDOW_WIDTH;
  double dy = (ymax-ymin)/WINDOW_HEIGHT;
  double x = 0, y = 0;

  complex<double> test( -.31415, .671); // -0.7796,0.134 );
 
 for(double xi = xmin; xi < xmax; xi += dx)
 {
    for(double yi = ymin; yi < ymax; yi += dy)
    {
      int iterations = 0;
      complex<double> com(xi, yi);

      // Determine if we diverge or not
      while( (abs(com) < 2) && (iterations < MAX_ITERS))
      {
	// f(x) = x^2 + c
        com = (pow( com, 2) + test );
        ++iterations;
      }

      // Get color      
      Color3d color(.0, .5, .4);
      color.rotateHue( ( (120) / MAX_ITERS) * iterations);

      // Render the pixel with appropriate color
      renderPixel(x, y, color );
      ++y;
    }
    
    ++x;
    y = 0;
  } 
}

void renderPixel(int x, int y, Color3d& color, float sz )
{
  glPointSize(sz);
  glColor3f(color.r, color.g, color.b );
  glBegin(GL_POINTS);
  glVertex2i(x, y);
  glEnd();
}

void GLrender()
{
  glClear(GL_COLOR_BUFFER_BIT);
  
  // Render the Julia Set based on the parameters for mins and maxs
  renderJuliaSet( -2, 2, -2, 2 );

  glFlush();	
  glutSwapBuffers();
}

void GLInit(int* argc, char** argv)
{
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutCreateWindow("Lab 8 - Sergio Morales");
  glutDisplayFunc(GLrender);
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  glOrtho(VIEW_LEFT, VIEW_RIGHT, VIEW_BOTTOM, VIEW_TOP, VIEW_NEAR, VIEW_FAR);
}


