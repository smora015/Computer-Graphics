// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Lab: 1
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <iostream>
#include <math.h>

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r = 1.0, float g = 1.0, float b = 1.0)
{
  glBegin( GL_POINTS );

  // Render a point based off a vertex
  glVertex2f( x, y);

  glEnd();
}

// Renders line based on incremental DDA
void render_line( float x0, float y0, float x1, float y1)
{
  // Declare variables
  float dy = (y1 - y0); // Change in y
  float dx = (x1 - x0); // Change in x
  float m = 0;

  // Find slope
  if( dx != 0 )
      m = dy / dx;

  float x = x0, y = y0; // Set initial x and y positions
 
  // Use 3 cases for rendering a line.
  if( abs( (int) m) > 1 ) // Case 1: Slope > 1 || Slope < -1 (steep lines)
  {
    while( y != y1)
    {
      renderPixel(  (int) round( x ),  (int) y );
      x +=  1/m; // Since the slope is greater than 1, we move x at a rate of 1/slope
      y += (dy < 0) ? -1 : 1; // Increment or decrement y based on direction
    }
  }
  else if ( m == 0 ) // Case 2: Slope is 0. Horizontal or Vertical lines
  {
    while( x != x1 || y != y1 )
    {
      renderPixel(  (int) round( x ),  (int) round( y ) );
      if( y != y1 )
	y += (dy < 0) ? -1 : 1;
      if( x != x1 )
	x += (dx < 0) ? -1 : 1;
    }
  }
  else				 // Case 3:  0 > Slope <= 1. Normal diagonal lines
  {
    while( x != x1 )
    {
      renderPixel(  (int) x,  (int) round( y ) );
      y += m;

      x += (dx < 0) ? -1 : 1;
    }
  }
}

float f2( float x, float cx, float y, float cy, float r )
{
  return ( (x - cx)*(x - cx) + (y - cy)*(y - cy) - r*r );
}

void render_circle ( float cx, float cy, float r )
{
  float x = cx; // x starts at the center
  float y = cy + round(r); // y starts at the center plus the radius

  // Render the first octant and perform matrix transpose
  while ( (x-cx) <= (y - cy) )
  {    
    // Render original pixel, then transpose on to the remaining 7 octants
    renderPixel(  (int) x,  (int) round(y) );
    renderPixel(  (int) (-x + 2*cx),  (int) (y ) ); // Reflect along x
    renderPixel(  (int) (x),  (int) (-y + 2*cy ) ); // Reflect along y
    renderPixel(  (int) (-x + 2*cx),  (int) (-y + 2*cy ) ); // Reflect along x and y
	renderPixel( (int) (cx + (y-cy)), (int) ( cy + (x-cx)) ); // Rotate 45 degrees
	renderPixel( (int) (cx + (y-cy)), (int) ( cy - (x-cx)) ); // Reflect along rotations
	renderPixel( (int) (cx - (y-cy)), (int) ( cy + (x-cx)) ); // "
	renderPixel( (int) (cx - (y-cy)), (int) ( cy - (x-cx)) ); // "
    
    // Incremement x, then decrement y based on condition
    ++x;
    if( (r*r - (x - cx)*(x - cx)- ((y - cy) - 1/2)*((y - cy) - 1/2)) < 0 ) // If point falls inside circle
      --y; // Then we move to lower pixel
  }
}


// Output functon to OpenGL Buffer
void GL_render()
{
  glClear( GL_COLOR_BUFFER_BIT );

  // NEW: Call the render function
  render_line( 400, 400, 400, 600 );
  render_line( 400, 400, 600, 600 );
  render_line( 400, 400, 600, 400 ); 
  render_line( 400, 400, 600, 200 ); 
  render_line( 400, 400, 400, 200 );
  render_line( 400, 400, 200, 200 ); 
  render_line( 400, 400, 200, 400 );
  render_line( 400, 400, 200, 600 ); 


  render_circle( 400, 400, 50 );
  render_circle( 300, 300, 100 );
  render_circle( 500, 500, 100 );
  render_circle( 400, 250, 50 );

  glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

  // Create window
  glutCreateWindow("CS 130 - <Sergio Morales>");

  // NEW: Register rendering function
  glutDisplayFunc( GL_render );

  // The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
  // For the purposes of this lab, this is set to the number of pixels
  // in each dimension.
  glMatrixMode(GL_PROJECTION_MATRIX);
  glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
}

int main(int argc, char** argv)
{
  GLInit(&argc, argv);
  glutMainLoop();

  return 0;
}


