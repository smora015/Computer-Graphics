// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Lab: 3
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
vector< Point2D > points; // Holds all clicked points
typedef struct clr        // Holds different triangle colors
{
  double R;
  double G;
  double B;

  clr() : R(0), G(0), B(0) {}
  clr( const double & nR, const double & nG, const double & nB ) : R(nR), G(nG), B(nB) {}

} clr;
clr colors[5] = { clr(.5, .1, .1), clr(.1, .5, .1), clr(.1, .1, .5), clr(.5, .75, 1), clr(.3, .2, .6) };


// ################ //
// Helper Functions //
// ################ //

int ymax( int number ) // Determines max y value out of 3 points
{
    int max = number*3;
    
    for(int i = number*3; i < number*3 + 3; ++i)
    {
      if( points[i].y >= points[max].y )
	max = i;
    }
    
    return max;
}

int ymin( int number ) // Determines minimum x value out of 3 points
{
    int min = number*3;

    for(int i = number*3; i < number*3 + 3; ++i)
    {
        if( points[i].y <= points[min].y )
            min = i;
    }
    
    return min;
}

// ############################################# //
// ########## Function Declarations ############ //
// ############################################# //
// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r, float g, float b);

// Renders line based on incremental DDA
void renderLine( int x0, int y0, int x1, int y1);

// Renders a triangular polygon
void renderTPolygon( int number );

// Returns appropriate window coordinates
Point2D GLScreenToWindowCoordinates( const Point2D & p );

// Output functon to OpenGL Buffer
void GL_render();

// Mouse event handler
void GL_mouse(int button, int state, int x, int y);

// Keyboard event handler
void GL_key( unsigned char key, int x, int y );

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv);



// ############################################ //
// ##############    Main    ################# //
// ############################################ //
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
  glutCreateWindow("CS 130 - Lab 3 <Sergio Morales>");

  // Register rendering function
  glutDisplayFunc( GL_render );
	
  // NEW: Register mouse event function
  glutMouseFunc( GL_mouse );

  // NEW: Register keyboard event function
  glutKeyboardFunc( GL_key );

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
  if( ( (points.size() % 3) == 0) && points.size() > 0)
  {
    for( int i = 0; i < points.size()/3; ++i)
      renderTPolygon(i); // Render triangle polygon
      
  }
  glutSwapBuffers();
}


// Finds mid y value out of 3 points, or returns false if not found
bool find_ymid( Point2D & mid, Point2D y_min, Point2D y_max, int number)
{
  for( int i = number*3; i < number*3 + 3; ++i )
  {
    if( (points[i].y > y_min.y) && (points[i].y < y_max.y) )
    {
      mid = points[i];
      return true;
    }
  }

  return true;
}

// Renders a triangle using 3 
void renderTPolygon( int number )
{
    // Find min and max values for polygon points
    int min = ymin( number );
    int max = ymax( number );
  
    Point2D y_min = points[min];
    Point2D y_max = points[max];
  
    // Triangle check 
    if( y_min.y == y_max.y )
    {
        cout << "Not a real triangle!" << endl;
        return;
    }

    Point2D mid;
    if( find_ymid( mid, y_min, y_max, number ) ) // Case 1: Split triangle into 2 
    {
      // Fetch Pnew. Confirmed to work.
      double s = (y_max.x - y_min.x) / ( y_max.y - y_min.y );
      double c = (y_min.x - s*y_min.y );
      Point2D p_new( s*mid.y + c, mid.y );

      // Triangle one (using p_new, p_mid, and y_max )
      double m_left = ( p_new.y - y_max.y ) / ( p_new.x - y_max.x);
      double m_right = ( mid.y - y_max.y ) / (mid.x - y_max.x ); 
      for ( int y = (int)mid.y; y < y_max.y; ++y)
      {
	// x = (1/m)*(y - b ) + x0
	double x0 =  (y - y_max.y) / m_right + y_max.x;
	double x1 =  (y - y_max.y) / m_left + y_max.x;

	// Swap x values if going in descending order
	if( x0 > x1 )
	{
	  double temp = x0;
	  x0 = x1;
	  x1 = temp;
	}

	// Render pixels, starting from left slope of the triangle, to the right slope.
	for ( int x = (int) x0; x < x1; ++x )
	{
	  renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
	}
      }

      // Triangle two ( using pnew, pmid, pmin )
      m_left = ( p_new.y - y_min.y ) / ( p_new.x - y_min.x );
      m_right = ( mid.y - y_min.y ) / ( mid.x - y_min.x );
      for ( int y = (int)y_min.y; y < mid.y; ++y)
      {
	// x = (1/m)*(y - b ) + x0
	double x0 =  (y - y_min.y) / m_left + y_min.x;
	double x1 =  (y - y_min.y) / m_right + y_min.x;
	
	// Swap x values if going in descending order
	if( x0 > x1 )
	{
	  double temp = x0;
	  x0 = x1;
	  x1 = temp;
	}

	// Render pixels
	for ( int x = (int) x0; x < x1; ++x )
	{
	  renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
	}
      }
      
    }
    else // Case 2 and 3: Triangle has 2 ymins or 2 ymaxs
    {
      cout << "No midpoint!" << endl;

      int p_newn = number*3;
      for(int i = number*3; i < number*3 + 3; ++i)
      {
	if( ((i != min)) && ((i != max)) )
	  p_newn = i;
      }

     Point2D p_new = points[p_newn];

      // Triangle one (using p_new, p_mid, and y_max )
 
     double m_left = ( p_new.y - y_min.y ) / ( p_new.x - y_min.x);
     double m_right = ( y_max.y - y_min.y ) / (y_max.x - y_min.x ); 

     for ( int y = (int)y_min.y; y < y_max.y; ++y)
      {
	// x = (1/m)*(y - b ) + x0
	double x0 =  (y - y_max.y) / m_right + y_max.x;
	double x1 =  (y - y_max.y) / m_left + y_max.x;

	// Swap x values if going in descending order
	if( x0 > x1 )
	{
	  double temp = x0;
	  x0 = x1;
	  x1 = temp;
	}

	// Render pixels, starting from left slope of the triangle, to the right slope.
	for ( int x = (int) x0; x < x1; ++x )
	{
	  renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
	}
      }
      
    }
    

}

void GL_mouse(int button, int state, int x, int y)
{
  // Left click gets points
  if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
  {
    // Get first/second coordinate
    Point2D temp = GLScreenToWindowCoordinates( Point2D(x, y) );
    points.push_back( temp );

    // After second click, form Bezier curve
    if( points.size() % 3 == 0 && points.size() > 0)
      glutPostRedisplay();
  }
  
  return;
} 

void GL_key( unsigned char key, int x, int y )
{
  if( key == '0' )
  {
    points.resize(0);
    glutPostRedisplay();
  }

  return;
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
    int temp = x0;
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
