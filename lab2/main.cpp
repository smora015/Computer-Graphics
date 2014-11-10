// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Lab: 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////

#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "point2d.h"

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

// Global vars used for constructing Bezier curves
vector< Point2D > points;

// ############################################# //
// ########## Function Declarations ############ //
// ############################################# //

// Factorial function
inline int fact( int x )
{
  // cout << "factorial of " << x << " is ";
  if( x <= 1 )
  {
    //cout << 1 << endl;
    return 1;
  }
  int fact = 1;

  for ( int i = 1; i <= x; ++i )
    fact *= i;

  //cout << fact << endl;
  return fact;
}

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r, float g, float b);

// Renders line based on incremental DDA
void renderLine( int x0, int y0, int x1, int y1);

// Finds B(t)
Point2D Bezier( double t );

// Renders Bezier quadratic curve
void renderBezier();

// Returns appropriate window coordinates
Point2D GLScreenToWindowCoordinates( const Point2D & p );

// Output functon to OpenGL Buffer
void GL_render();

// Mouse even handler
void GL_mouse(int button, int state, int x, int y);

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv);



int main( int argc, char** argv )
{
  GLInit(&argc, argv);
  glutMainLoop();
}


// ############################################# //
// ########## Function Definitions ############# //
// ############################################# //
void GLInit(int* argc, char** argv)
{
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

  // Create window
  glutCreateWindow("CS 130 - Lab 2 <Sergio Morales>");

  // Register rendering function
  glutDisplayFunc( GL_render );
	
  // NEW: Register mouse event function
  glutMouseFunc( GL_mouse );


  // The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
  // For the purposes of this lab, this is set to the number of pixels
  // in each dimension.
  glMatrixMode(GL_PROJECTION_MATRIX);
  glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
}

void GL_render()
{
  glClear( GL_COLOR_BUFFER_BIT );

  // Render Bezier curve after 2 mouse clicks on screen
  if( points.size() > 1 )
    renderBezier();

  glutSwapBuffers();
}

void GL_mouse(int button, int state, int x, int y)
{
  // Right click clears the screen and the vector of points
  if( button == GLUT_RIGHT_BUTTON && state == GLUT_UP )
  {
    points.clear();
    glutPostRedisplay();
  }

  // Left click gets points
  if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
  {
    // Get first/second coordinate
    Point2D temp = GLScreenToWindowCoordinates( Point2D(x, y) );
    points.push_back( temp );

    // After second click, form Bezier curve
    if( points.size() > 1 )
      glutPostRedisplay();
  }
  
  return;
} 

// Bezier generalization for nth-order 
Point2D Bezier( double t )
{
  Point2D b(0,0);

  int n = points.size();
  for ( int i = 0; i < n; ++i)
  {
    // Summation from 0 -> # of points.  ( n-1 C i ) * t^i * ( 1 - t ) ^ (n-1 - 1)
    b += points[i] * static_cast<double>((fact(n-1) / (fact(i) * fact( n-1 - i ) )) * pow( t, i ) * pow( 1 - t, n-1 - i )); 
  }
  
  return b;

}
void renderBezier()
{
  Point2D a = points[0], b = points[0]; // Initialize endpoints

  int n = points.size(); // The degree of the Bezier curve

  // Use small t increments
  for( double t = 0; t <= 1; t += 0.01 )
  {
    b = Bezier( t ); // Find Bezier curve end point for specific t
    renderLine( a.x, a.y, b.x, b.y ); // Start from old b or first point, up until new b to draw a line
    a = b;          // a gets the old Bezier end point, to connect all the lines to form the curve

    renderPixel( b.x, b.y, 1, 0, 0);
  }

}

void renderPixel(int x, int y, float r, float g, float b)
{
  // Set glBegin flags
  glBegin( GL_POINTS );

  // Change color to RGB parameter
  glColor3f( r, g, b );

  // Render a point based off a vertex
  glVertex2f( x, y);

  // Signal finish rendering
  glEnd();
}

void renderLine( int x0, int y0, int x1, int y1)
{
  // Swap endpoints to make rendering go from left -> right
  if( x0 > x1 )
  {
  	float temp = x0;
  	x0 = x1;
  	x1 = temp;
  	
  	temp = y0;
  	y0 = y1;
  	y1 = temp;
  }
  

  // Declare variables
  float dy = (y1 - y0); // Change in y
  float dx = (x1 - x0); // Change in x
  float m = 0;

  // Find slope
  if( dx != 0 )
      m = dy / dx;

  float x = x0, y = y0; // Set initial x and y positions
 
  if( m > 1 ) // Case 1: Slope > 1 (steep lines)
  {
    while( abs( (int) y) != abs( (int)y1))
    {
      renderPixel(  (int) round( x ),  (int) y, 1, 1, 1 );
      x +=  1/m; // Since the slope is greater than 1, we move x at a rate of 1/slope
      y += (dy < 0) ? -1 : 1; // Increment or decrement y based on direction      
    }
  }
  else if( m < -1 ) // Case 2: Slope < -1 (steep lines)
  {
    while( abs( (int) y) != abs( (int) y1))
    {
      renderPixel(  (int) round( x ),  (int) y, 1, 1, 1 );
      x -=  1/m; // Since the slope is less than 1, we move x at a rate of -1/slope
      y += (dy < 0) ? -1 : 1; // Increment or decrement y based on direction
    }
  }
  else if ( m == 0 ) // Case 3: Slope is 0. Horizontal or Vertical lines
  {
    while( abs( (int) x) != abs( (int) x1) || abs( (int) y) != abs( (int)y1) )
    {
      renderPixel(  (int) round( x ),  (int) round( y ), 1, 1, 1);
      if( y != y1 )
		y += (dy < 0) ? -1 : 1;
      if( x != x1 )
		x += (dx < 0) ? -1 : 1;
    }
  }
  else if( m > -1 && m <= 1) // Case 4:  -1 < Slope <= 1. Normal diagonal lines
  {
    while( abs( (int) x) != abs( (int) x1) )
    {
      renderPixel(  (int) x,  (int) round( y ), 1, 1, 1 );
      y += m;

      x += (dx < 0) ? -1 : 1;
    }
  }
}

Point2D GLScreenToWindowCoordinates( const Point2D & p )
{
  GLdouble model[16];
  glGetDoublev( GL_MODELVIEW_MATRIX, model );

  GLdouble projection[16];
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
	
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
	
  Point2D converted;
  GLdouble temp;
	
  gluUnProject( p.x, viewport[3] - p.y, 0, model, projection, viewport, &converted.x, &converted.y, &temp );
	
  return converted;
}
