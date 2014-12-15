// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Project: 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>
#include "point2d.h"
#include "point3d.h"

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

// Pan and Zoom for interface
#define PAN_LEFT_RIGHT 0   // Parameter values to specify type of panning
#define PAN_UP_DOWN 1      // ""
#define DEPTH_INIT -10000  // This is our initial z-value

// Phong Illumination Model
#define Ka .1//.1  // Coefficients
#define Kd .3//.7
#define Ks .7// 1

#define La .4//.2  // Intensities
#define Ld .7//.8
#define Ls .9//.1

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
  Triangle() : p1(0), p2(0), p3(0), Color( clr(1, 1, 1) ) {}
  Triangle( const int & p11, const int & p22, const int & p33, clr colorr )
             : p1(p11), p2(p22), p3(p33), Color( colorr ) {}
} Triangle;

typedef struct Ray
{
    Point3D origin;
    Point3D direction;
    Ray() {Point3D d(0.0, 0.0, 1.0);  direction = d;}
    Ray(const Point3D& o, const Point3D& dir)
    {
        origin = o;
        Point3D d(0.0, 0.0, 1.0);
        float mag = dir.magnitude();
        if (mag > 0.0) {d = dir;}
        direction = d;
    }
} Ray;

/// --------------------------- Global Variables -------------------------///
vector< Point2D > points;           // Container for 2D points. Used for scan conversion
vector< Point3D > cartesian_points; // Container of 3D points. Raw point data
vector< Triangle > triangles;       // Container of 3 sets of 3D points, or a triangle. Used to get points array data
Point3D lights[5]  = {Point3D( 200, 200, -5/*-.5, -.5, -.5*/ ), Point3D( 1, -5, -2 ), Point3D( .3, .5, -1), Point3D( 1, 1, -1),
		      Point3D( 1, -1, 1) };


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
double degrees_radians( double deg );

// Rotates the screen by the specified angle
void rotate_model( double t );

// Moves the screen left-right, or up-down
void pan_model( int mode, double amount );

// Zooms the model in or out
void zoom_model( double amount );

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r, float g, float b);


/// ----------------------- Main --------------------------------------///
int main( int argc, char** argv )
{	
  if( argc < 2)
  {
    cerr << "Error: File not specified!" << endl; 
    return -1;
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
        
        triangles.push_back( Triangle( x, y, z, clr(1,1,1)/*colors[ (color_index++) % 10 ]*/) );
    }
    
    // Display usage information
    cout << "============== Flat/Smooth Shaded Rendering ==============\n"
         << "                   ..Controls.. \n"
         << " Rotate:    'right' -> RIGHT  | 'left' -> LEFT\n"
         << " Pan:       'w'     -> UP     | 'a' -> LEFT \n"
         << "            's'     -> DOWN   | 'd' -> RIGHT \n"
         << " Zoom:      '-'     -> OUT    | '+' -> IN \n"
         << "=========================================================="
         << endl;
         
    model.close();
}

void GLInit(int* argc, char** argv)
{
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    
    // Create window
    glutCreateWindow("CS 130 - Project 2 <Sergio Morales>");
    
    // Register rendering function
    glutDisplayFunc( GL_render );
    
    // Register special keyboard event function
    glutSpecialFunc( GL_spec_key );
    
    // Register keyboard event function
    glutKeyboardFunc( GL_key );
    
    // The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
    // For the purposes of this lab, this is set to the number of pixels
    // in each dimension.
    glClearColor( 1.0, 1.0, 1.0, 0 ); 
    glMatrixMode(GL_PROJECTION_MATRIX);
    glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
}


// Finds the determinant of a matrix made up of 3 vectors (3x3)
float tri_det( const Point3D &v1, const Point3D &v2, const Point3D &v3 )
{
  return ( ((v1.x * v2.y * v3.z) + (v2.x * v3.y * v1.z) + (v3.x * v1.y * v2.z )) - // Lower diaganoal
	   ((v1.z * v2.y * v3.x) + (v2.z * v3.y * v1.x) + (v3.z * v1.y * v2.x)));  // Upper diagonal
}
	
// Returns true if the ray intersects the triangle.
bool rayIntersectsTriangle(const Ray& ray, const Triangle& triangle, float* t )
{
  if( t != NULL ) 
    *t = -1;

  // Set up our points in triangle
  Point3D a = cartesian_points[ triangle.p1 ]; // First point        
  Point3D b = cartesian_points[ triangle.p2 ]; // Second point
  Point3D c = cartesian_points[ triangle.p3 ]; // Third point
  

  Point3D Eb = ( Point3D(b.x - a.x, b.y - a.y, b.z - a.z) );
  Point3D Ec = ( Point3D(c.x - a.x, c.y - a.y, c.z - a.z) );

  // Find t, beta, and gamma
  // t = det(T) / det(A), where A, B, and C are matrices populated with Eb, Ec, and ray.
  // beta = det( B ) / det( A ) and gamma = det( C ) / det( A )
  float det_a = tri_det( Point3D(-Eb.x, -Eb.y, -Eb.z), 
			 Point3D(-Ec.x, -Ec.y, -Ec.z), 
			 (ray.direction) );

  *t = tri_det( Point3D(-Eb.x, -Eb.y, -Eb.z), 
		Point3D(-Ec.x, -Ec.y, -Ec.z), 
		Point3D(a.x - ray.origin.x, a.y - ray.origin.y, a.z - ray.origin.z ) ) / det_a;            // det(T) / det(a)

  // Determine if we intersect using all 3 values
  if( *t > 0 )
  {
    float beta = tri_det( Point3D(a.x - ray.origin.x, a.y - ray.origin.y, a.z - ray.origin.z), 
			  Point3D(-Ec.x, -Ec.y, -Ec.z), 
			  ray.direction ) / det_a;       // det(B) / dat(a)
    float gamma = tri_det( Point3D(-Eb.x, -Eb.y, -Eb.z), 
			   Point3D(a.x - ray.origin.x, a.y - ray.origin.y, a.z - ray.origin.z), 
			   ray.direction ) / det_a;  // det(C) / det(a)

    if( (beta >= 0) && (gamma >= 0) )
    {
      if( (beta + gamma) <= 1 )
	return true;
    }
  }


  return false;
}

//Returns -1 if no sphere should appear in pixel, otherwise returns index of sphere.
int pixelOn( Ray current )
{
    int num_triangles = (int) triangles.size();
    float min = -1.0;
    int min_i = -1;

    // Two t's used to check for intersection
    float t = 0;

    for (int i = 0; i < num_triangles; i++)
    {
      if( rayIntersectsTriangle( current, triangles.at( i ), &t) )
      {
	// Determine the closest intersection
	if( t >= 0 )
	{
	  // If we found an intersection, update index
	  if( min == -1 || (t <= min) )
	  {
	    min = t;
	    min_i = i;
	  }
	}
      }
    }
    return min_i;
}


void GL_render()
{
    glClear( GL_COLOR_BUFFER_BIT );
    
    // For each pixel,
    for( int i = 0; i < WINDOW_WIDTH; ++i )
    {
      for( int j = 0; j < WINDOW_HEIGHT; ++j )
      {
		// Check if there is an intersection between ray and a triangle
		int which_triangle = pixelOn( Ray( Point3D( i, 800 - j, 0), Point3D(0, 0, 1) ) );
 		Point3D viewer(i, 800 - j, 1.0);
	
		if( which_triangle >= 0 )
		{	  
		  // We now find the illumination
		  float sum = 0;
	  
		  // For all light sources
		  for( int x = 0; x < 1; ++x )
		  {
		    // We determine if there is a shadow. If so, only add ambient light
			int shadow = pixelOn( Ray( Point3D( i, 800 - j, 0), Point3D(lights[x].x - i, lights[x].y, lights[x].z - 0) ) );
		    
		    float diffuse = 0, specular = 0, ambient = 0;
		    // If there isn't, add diffuse and specular as well
		    if( shadow < 0 )
		    {
		      // Phong Illumination Model
		      // L = kd * Ld * cos(theta) + ks * Ls * cos( phi ) + ka * La
		      // cos( theta ) = N * l / (||N|| || L || ), where N = normal vector, L = light
		    
  		      // First we find the normal vector of that triangle
		      Point3D p1 = cartesian_points.at( triangles.at( which_triangle ).p1 );  // Get x, y, z
		      Point3D p2 = cartesian_points.at( triangles.at( which_triangle ).p2 );
		      Point3D p3 = cartesian_points.at( triangles.at( which_triangle ).p3 );
		    
		      Point3D vector1 = Point3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);       // Get 2 vectors from pts
		      Point3D vector2 = Point3D(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
		    
		      Point3D normal = vector1.cross( vector2 );                        // Find normal from cross
		    
		      // Using the normal, we calculate the diffuse component as Kd * Ld * ( N * l )/(|N||L|)
		      diffuse = (normal.dot( lights[x] ) / ( normal.magnitude() * lights[x].magnitude() ))*Kd*Ld;
		    
		      // We then find our specular component using the reflected vector from the light source
		      Point3D reflector = normal * 2 * ( normal.dot( (lights[x] / lights[x].magnitude()) ) ) - lights[x];
		      specular = (reflector.dot( viewer ) / ( reflector.magnitude() * viewer.magnitude() ))*Ks*Ls;

		      // Finally, the ambient component
		      ambient = Ka * La;

		    }
		    else
		    {
		      // Finally, the ambient component
		      ambient = Ka * La;
		    }
		      // Find the total illumination
		      sum += diffuse + specular + ambient;
		    
		  }
	  
	  	// After looping through all light sources, we then render the pixel
	    clr colour = triangles.at( which_triangle ).Color;	  	  
	    renderPixel( i, 800 - j,/*colour.R**/sum, /*colour.G**/sum, /*colour.B**/sum);

		}
      
      }
    }
 
    glutSwapBuffers();
}

void GL_spec_key( int key, int x, int y )
{
    // Clear points after rendering entire model
    points.clear();
    
    if(key == GLUT_KEY_RIGHT) // Rotate along y-axis -->
    {
        //cout << "Rotated right!" << endl;
        rotate_model( 25 );
    }
    else if(key == GLUT_KEY_LEFT) // Rotate along y-axis <--
    {
        //cout << "Rotated left!" << endl;
        rotate_model( -25 );
    }
    
    // Refresh screen
    glutPostRedisplay();
}

void GL_key( unsigned char key, int x, int y )
{
    // Clear after rendering entire model
    points.clear();
     
    if( key == '+' )
    {
        zoom_model( 1.2 );
        //cout << "Zooming in!" << endl;
    }
    else if( key == '-' )
    {
        zoom_model( .8 );
        //cout << "Zooming out!" << endl;
    }
    else if( key == 'a' )
    {
        pan_model( PAN_LEFT_RIGHT, -50 );
        //cout << "Panning left!" << endl;
    }
    else if( key == 'd' )
    {
        pan_model( PAN_LEFT_RIGHT, 50 );
        //cout << "Panning right!" << endl;
    }
    else if( key == 'w' )
    {
        pan_model( PAN_UP_DOWN, 50 );
        //cout << "Panning up!" << endl;
    }
    else if( key == 's' )
    {
        pan_model( PAN_UP_DOWN, -50 );
        //cout << "Panning down!" << endl;
    }
    
    // Refresh screen
    glutPostRedisplay();
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
    return deg*3.14/180; // Convert by multiplying degrees by (PI / 180 degrees)
}

void rotate_model( double t )
{
    Point3D avg = average_vertices(); // Get average of all vertices to translate polygon to origin
    for( int i = 0; i < cartesian_points.size(); ++i ) // Displace vertices to axis
        cartesian_points[i] -= avg;

    for( int i = 0; i < cartesian_points.size(); ++i ) // Rotate vertices through translation
    {
        Point3D rotate;
        rotate.x = cartesian_points[i].x*cos( degrees_to_radians(t)) + cartesian_points[i].z*sin( degrees_to_radians(t));
        rotate.y = cartesian_points[i].y;
        rotate.z = -cartesian_points[i].x*sin( degrees_to_radians(t)) + cartesian_points[i].z*cos( degrees_to_radians(t));
        cartesian_points[i] = rotate;
    }

    for( int i = 0; i < cartesian_points.size(); ++i ) // Translate back to original position
        cartesian_points[i] += avg;
}

void pan_model( int mode, double amount )
{
    Point3D avg = average_vertices(); // Get average of all vertices to translate polygon to origin
    for( int i = 0; i < cartesian_points.size(); ++i ) // Displace vertices to axis
        cartesian_points[i] -= avg;
    
    for( int i = 0; i < cartesian_points.size(); ++i ) // Pan by amount using translation
    {
        if( mode == PAN_LEFT_RIGHT )
            cartesian_points[i].x += amount;
            
        else if( mode == PAN_UP_DOWN )
            cartesian_points[i].y += amount;
    }

    for( int i = 0; i < cartesian_points.size(); ++i ) // Translate back to original position
        cartesian_points[i] += avg;
}

void zoom_model( double amount )
{
    Point3D avg = average_vertices(); // Get average of all vertices to translate polygon to origin
    for( int i = 0; i < cartesian_points.size(); ++i ) // Displace vertices to axis
        cartesian_points[i] -= avg;
    
    for( int i = 0; i < cartesian_points.size(); ++i ) // Rotate by scaling everything
        cartesian_points[i] *= amount;
    
    for( int i = 0; i < cartesian_points.size(); ++i ) // Translate back to original position
        cartesian_points[i] += avg;
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



