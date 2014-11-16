// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Project: 1
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////

#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include "point2d.h"
#include "point3d.h"

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

/// ------------------------- Struct/Type Definitions ---------------------///
typedef struct clr // Holds different triangle colors
{
  double R;
  double G;
  double B;

  clr() : R(0), G(0), B(0) {}
  clr( const double & nR, const double & nG, const double & nB ) 
        : R(nR), G(nG), B(nB) {}

} clr;

typedef struct Triangle
{
  int p1; // Holds the index to the actual 3-Dimensional point
  int p2;
  int p3;
  clr Color;
  
  Triangle() : p1(0), p2(0), p3(0) {}
  Triangle( const int & p11, const int & p22, const int & p33 )
        : p1(p11), p2(p22), p3(p33) {}
  
} Triangle;

/// --------------------------- Global Variables -------------------------///
const int DEPTH_INIT = -999; // This is our initial z-value
clr colors[5] = { clr(.5, .1, .1), clr(.1, .5, .1), clr(.1, .1, .5), clr(.5, .75, 1), clr(.3, .2, .6) }; // Array of 5 colors to cycle through
vector< Point2D > points; // Container for 2D points. Used for scan conversion
vector< Point3D > cartesian_points; // Container of 3D points. Raw point data
vector< Triangle > triangles; // Container of 3 sets of 3D points, or a triangle. Used to get points array data
vector< vector< double > > zbuffer( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) ); // Initialize our zbuffer to 800x800 points at DEPTH_INIT


/// -------------------------- Function Declarations --------------------///
// Grab all model information
void parse_model( string file );

// Output functon to OpenGL Buffer
void GL_render();

// Mouse event handler
void GL_mouse(int button, int state, int x, int y);

// Keyboard event handler
void GL_key( unsigned char key, int x, int y );

// Keyboard event handler
void GL_spec_key( int key, int x, int y );

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv);

// Returns appropriate window coordinates
Point2D GLScreenToWindowCoordinates( const Point2D & p );

// Gets the average of all vertices in our points from model
Point3D average_vertices();

// Converts degrees to radians
double get_rad( double deg );

// Rotates the screen by the specified angle	
void rotate( double q );

// Determines max y value out of 3 points.
// @param: number   The starting point in our points array for the 3 vertices.
int ymax( int number );

// Determines minimum y value out of 3 points
// @param: number   The starting point in our points array for the 3 vertices
int ymin( int number );

// Renders a triangular polygon
// @param: number   The starting point in our points array for the 3 vertices
void renderTPolygon( int number );

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r, float g, float b);

// Renders line based on incremental DDA
void renderLine( int x0, int y0, int x1, int y1);


/// ----------------------- Main --------------------------------------///
int main( int argc, char** argv )
{
  if( argc < 2)
  {
    cerr << "Error: File not specified!" << endl; return -1;
  }

  // Parse through model file
  parse_model( argv[1] );
  
  GLInit(&argc, argv);
  glutMainLoop();
}


/// -------------------------- Function Definitions -------------------///

void parse_model( string file )
{
    // Open model file
    fstream model( file.c_str() , std::fstream::in | std::fstream::out );
    
    short num_points, num_triangles;
    
    // Get number of points and number of triangles
    model >> num_points;
    model >> num_triangles;
    
    int data;
    int x,y,z;
    for( int i = 0; i < num_points; ++i )
    {
				// Populate the cartesian points container
        model >> x;
        model >> y;
        model >> z;
        cartesian_points.push_back( Point3D( x, y, z) );
    }    
    
    for( int i = 0; i < num_triangles; ++i )
    {
				// Populate the triangles container with cartesian_points indices
        model >> x;
        model >> y;
        model >> z;
        triangles.push_back( Triangle( x, y, z) );
    }
    
    model.close();
}


void GLInit(int* argc, char** argv)
{
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

  // Create window
  glutCreateWindow("CS 130 - Project 1 <Sergio Morales>");

  // Register rendering function
  glutDisplayFunc( GL_render );
	
  // Register special keyboard event function
  glutSpecialFunc( GL_spec_key );

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

  // Traverse through points in all polygons
  for( int i = 0; i < triangles.size(); ++i )
  {
			// Grab the 2D points and push them to points container
      Point3D temp;
    
      temp = cartesian_points[ triangles[i].p1 ];
      points.push_back( Point2D( temp.x, temp.y ) );
      cout << "p1: " << temp.x << " " << temp.y << " " << temp.z << endl;
      

			temp = cartesian_points[ triangles[i].p2 ];
      points.push_back( Point2D( temp.x, temp.y ) );
			cout << "p2: " << temp.x << " " << temp.y << " " << temp.z << endl;


      temp = cartesian_points[ triangles[i].p3 ];
      points.push_back( Point2D( temp.x, temp.y ) );
      cout << "p3: " << temp.x << " " << temp.y << " " << temp.z << endl << endl;
  }

	// Go through points container and output pixels, based off zbuffer as well
  for( int i = 0; i < points.size()/3; ++i)
  	renderTPolygon(i); // Render triangle polygon
      
	cout << "Triangles size: " << triangles.size() <<  " points size: " << points.size() << endl;
	points.clear();
  glutSwapBuffers();
}


void GL_spec_key( int key, int x, int y )
{
  if(key == GLUT_KEY_RIGHT) // Rotate along y-axis -->
	{
		cout << "Rotated left!" << endl;
		rotate( -25 );
    points.clear();
		glutPostRedisplay();
	}
	else if(key == GLUT_KEY_LEFT) // Rotate along y-axis <--
	{
		cout << "Rotated right!" << endl;
		rotate( 25 );
		points.clear();
		glutPostRedisplay();
	}
}

void GL_key( unsigned char key, int x, int y )
{
  if( key == '0' )
  {
		//cartesian_points.clear();
    points.clear();
		//triangles.clear();
		zbuffer.clear();
		zbuffer = vector< vector<double> >( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) );

    glutPostRedisplay();
  }
	else if( key == '+' )
  {
		//cartesian_points.clear();
    points.clear();
		//triangles.clear();
		//zbuffer.clear();
		//zbuffer = vector< vector<double> >( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) );

		cout << "Zoomed in!" << endl;
    glutPostRedisplay();
  }
	else if( key == '-' )
  {
		//cartesian_points.clear();
    points.clear();
		//triangles.clear();
		//zbuffer.clear();
		zbuffer = vector< vector<double> >( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) );

    glutPostRedisplay();

		cout << "Zoomed out!" << endl;
	}

  return;
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


Point3D average_vertices()
{
	Point3D average;
	for( int i = 0; i < cartesian_points.size(); ++i)
		average += cartesian_points[i];

	float div = cartesian_points.size();
	average /= div;
	return average;

}

double degrees_to_radians( double deg )
{
	return deg*3.14/180;
}
	
void rotate( double q )
{
	Point3D avg = average_vertices(); // Get average of all vertices to translate polygon to origin

	for( int i = 0; i < cartesian_points.size(); ++i ) // Displace vertices to axis
		cartesian_points[i] -= avg;

	for( int i = 0; i < cartesian_points.size(); ++i ) // Rotate
	{
		Point3D rotate;
		rotate.x = cartesian_points[i].x*cos( degrees_to_radians(q)) + cartesian_points[i].z*sin( degrees_to_radians(q));
		rotate.y = cartesian_points[i].y;
		rotate.z = -cartesian_points[i].x*sin( degrees_to_radians(q)) + cartesian_points[i].z*cos( degrees_to_radians(q));
		cartesian_points[i] = rotate;
	}

	for( int i = 0; i < cartesian_points.size(); ++i ) // Translate back to original position 
		cartesian_points[i] += avg;
	
		
}


bool update_zbuffer( double x, double y, int number )
{
	//first calculate plane of current triangle
	//do this by finding a normal vector, or the cross product of two vectors

	Point3D p1v = cartesian_points.at( triangles.at(number).p1 );//stores the values of each vertex in the triangle
	Point3D p2v = cartesian_points.at( triangles.at(number).p2 );
	Point3D p3v = cartesian_points.at( triangles.at(number).p3 );
	
	Point3D vector1 = p1v-p2v; //finds two vectors in the plane equal to the triangle
	Point3D vector2 = p1v-p3v;
	
	Point3D norm = vector1.cross( vector2 ); // the normal vector to the plane
	double D = norm.x*p1v.x + norm.y*p1v.y + norm.z*p1v.z;	
	double z = ( D -norm.x*x - norm.y*y )/norm.z; //will find the z val of the passed in x&y
	
	if( z > zbuffer.at( (x >= 800) ? 799 : x ).at( (y >= 800) ? 799 : y ) ||
			zbuffer.at( (x >= 800) ? 799 : x ).at( (y >= 800) ? 799 : y ) == DEPTH_INIT  )//point is only or closest point
	{
		//cout << "x: " << x << " y: " << y << " z: " << z << endl;
		zbuffer.at( (x >= 800) ? 799 : x ).at( (y >= 800) ? 799 : y ) = z; // we update the zbuff
		return true;// we want to render this point
	}


	return false; //there is already a closer point at this value
}


// ------ Triangle functions ------ //
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

int ymin( int number ) // Determines minimum y value out of 3 points
{
    int min = number*3;

    for(int i = number*3; i < number*3 + 3; ++i)
    {
        if( points[i].y <= points[min].y )
            min = i;
    }
    
    return min;
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

  return false;
}

// Renders a triangle using 3 points
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
      // Fetch Pnew.
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
						if( update_zbuffer(x, y, number) )
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
					if( update_zbuffer(x, y, number) )
     				renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
        }
      }

    }
    else // Case 2 and 3: Triangle has 2 ymins or 2 ymaxs
    {
      // Grab second y-min/y-max
      int p_newn = number*3;
      for(int i = number*3; i < number*3 + 3; ++i)
      {
        if( ((i != min)) && ((i != max)) )
        p_newn = i;
      }
      Point2D p_new = points[p_newn];
      
      // Two possibilities in this case:
      if( p_new.y == y_max.y)
      {
          // Triangle's horizontal line is on top
            
          // Swap if p_new.x is greater than y_max.x
          if( p_new.x > y_max.x )
          {
              Point2D temp = p_new;
              p_new = y_max;
              y_max = temp;
          }
          
          
          double m_left = ( p_new.y - y_min.y ) / ( p_new.x - y_min.x );
          double m_right = ( y_max.y - y_min.y ) / ( y_max.x - y_min.x );
          for ( int y = (int)y_min.y; y < y_max.y; ++y)
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
						if( update_zbuffer(x, y, number) )
 	           renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
          }
        }
      }
      else
      {
          // Triangle's horizontal line is on bottom

          // Swap if p_new.x is greater than y_max.x
          if( p_new.x > y_min.x )
          {
              Point2D temp = p_new;
              p_new = y_min;
              y_min = temp;
          }

          double m_left = ( p_new.y - y_max.y ) / ( p_new.x - y_max.x);
          double m_right = ( y_min.y - y_max.y ) / (y_min.x - y_max.x ); 
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
							if( update_zbuffer(x, y, number) )
  	            renderPixel(x, y, colors[number % 5].R, colors[number % 5].G, colors[number % 5].B);
            } 
              
          }
      }

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
