// Name: Sergio Morales
// Quarter, Year: Fall 2014
// Project: 1
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
#define Ka .3  // Coefficients
#define Kd .6
#define Ks .3

#define La .7  // Intensities
#define Ld .6
#define Ls .4

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
clr colors[10] = { clr(.5, .1, .1), clr(.1, .5, .1), clr(.1, .1, .5), 
		   clr(.5, .75, 1), clr(.3, .2, .6), clr(.3, .5, .6), 
		   clr(.5, .2, .5), clr(.5, .75, 0), clr(.5, .8, .1), 
		   clr(0, .75, 3) }; // Array of 10 colors to cycle through
vector< vector< double > > zbuffer( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) ); // Initialize our zbuffer to 800x800 points at DEPTH_INIT

bool WIREFRAME_TOGGLE = false;      // Toggles wireframe view on/off by pressing '0'
vector< Point2D > points;           // Container for 2D points. Used for scan conversion
vector< Point3D > cartesian_points; // Container of 3D points. Raw point data
vector< Triangle > triangles;       // Container of 3 sets of 3D points, or a triangle. Used to get points array data
unsigned char color_index = 0;          // Helps cycle through different colors for triangles
Point3D lights[5] = {Point3D( 0, 0, -.5 ), Point3D( 1, 5, -2 ), Point3D( 0, 0, -1), Point3D( 1, 1, -1),
		     Point3D( 0, -1, 1) };


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

// Determines if rendering should occur for x,y at a certain z value
bool update_zbuffer( double x, double y, int number );

// Determines max y value out of 3 points.
// @param: number The starting point in our points array for the 3 vertices.
int ymax( int number );

// Determines minimum y value out of 3 points
// @param: number The starting point in our points array for the 3 vertices
int ymin( int number );

// Renders a triangular polygon
// @param: number The starting point in our points array for the 3 vertices
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
        
        triangles.push_back( Triangle( x, y, z, colors[ (color_index++) % 10 ]) );
    }
    
    // Display usage information
    cout << "============== Z-buffer Model Rendering ==============\n"
         << "                   ..Controls.. \n"
         << " Rotate:    'right' -> RIGHT  | 'left' -> LEFT\n"
         << " Pan:       'w'     -> UP     | 'a' -> LEFT \n"
         << "            's'     -> DOWN   | 'd' -> RIGHT \n"
         << " Zoom:      '-'     -> OUT    | '+' -> IN \n"
         << " Wireframe: '0'     -> ON/OFF \n"
         << "======================================================"
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
int pixelOn(float pixel_x, float pixel_y)
{
    Ray current(Point3D(pixel_x, pixel_y, 0.0), Point3D(0.0, 0.0, 1.0));
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
	//cout << i << " " << j << " ";
	// Check if there is an intersection between ray and a triangle
	int which_triangle = pixelOn( i, 800 - j );
	
	if( which_triangle >= 0 )
	{
	  Point3D viewer(i, j, 1.0);
	  
	  // We now find the illumination
	  float sum = 0;

	  // For all light sources
	  for( int x = 0; x < 5; ++x )
	  {
	    // We determine if there is a shadow. If so, only add ambient light
	    
	    // If there isn't, add diffuse and specular as well
	    
	    // L = kd * Ld * cos(theta) + ks * Ls * cos( phi ) + ka * La
	    // cos( theta ) = N * l / (||N|| || L || ), where N = normal vector, L = light
	    Point3D p1 = cartesian_points.at( triangles.at( which_triangle ).p1 );  // Get x, y, z
	    Point3D p2 = cartesian_points.at( triangles.at( which_triangle ).p2 );
	    Point3D p3 = cartesian_points.at( triangles.at( which_triangle ).p3 );
	    
	    Point3D vector1 = Point3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);       // Get 2 vectors from pts
	    Point3D vector2 = Point3D(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
	    
	    Point3D normal = vector1.cross( vector2 );                        // Find normal from cross
	    
	    float diffuse = (normal.dot( lights[x] ) / ( normal.magnitude() * lights[x].magnitude() ))*Kd*Ld;
	    
	    Point3D reflector = normal * 2 * ( normal.dot( (lights[x] / lights[x].magnitude()) ) ) - lights[x];
	    float specular = (reflector.dot( viewer ) / ( reflector.magnitude() * viewer.magnitude() ))*Ks*Ls;
	    
	    float ambient = Ka * La;
	    
	    sum += diffuse + specular + ambient;

	  }
	    clr colour = triangles.at( which_triangle ).Color;	  	  
	    renderPixel( i, 800 - j, colour.R*sum, colour.G*sum, colour.B*sum);

	}
      
      }
    }

    // Project 1 way of rendering triangles
    /*// Traverse through points in all polygons
    for( int i = 0; i < triangles.size(); ++i )
    {
        // Grab the 2D projected points and push them to points container
        Point3D temp;
        temp = cartesian_points[ triangles[i].p1 ]; // First point
        points.push_back( Point2D( temp.x, temp.y ) );
        
        temp = cartesian_points[ triangles[i].p2 ]; // Second point
        points.push_back( Point2D( temp.x, temp.y ) );
        
        temp = cartesian_points[ triangles[i].p3 ]; // Third point
        points.push_back( Point2D( temp.x, temp.y ) );
    }
    
    // Go through points container and output pixels
    for( int i = 0; i < points.size()/3; ++i)
    renderTPolygon(i); // Render triangle polygons, based on zbuffer depth as well*/
    
    glutSwapBuffers();
}

void GL_spec_key( int key, int x, int y )
{
    // Clear points and zbuffer after rendering entire model
    points.clear();
    zbuffer.clear();
    zbuffer = vector< vector<double> >( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) );
    
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
    // Clear zbuffer after rendering entire model
    points.clear();
    zbuffer.clear();
    zbuffer = vector< vector<double> >( WINDOW_WIDTH, vector<double>( WINDOW_HEIGHT, DEPTH_INIT ) );
    
    if( key == '0' )
    {
        WIREFRAME_TOGGLE = (WIREFRAME_TOGGLE ? false : true );
    }
    else if( key == '+' )
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

bool update_zbuffer( double x, double y, int number )
{
    // Get values of each 3D point
    Point3D p1 = cartesian_points.at( triangles.at(number).p1 );
    Point3D p2 = cartesian_points.at( triangles.at(number).p2 );
    Point3D p3 = cartesian_points.at( triangles.at(number).p3 );
    
    // Get vectors in the plane
    Point3D vector1 = p1-p2;
    Point3D vector2 = p1-p3;
    
    // Find the vector normal to the plane
    Point3D normal_plane = vector1.cross( vector2 );
    double value_D = p1.x*normal_plane.x + p1.y*normal_plane.y + p1.z*normal_plane.z;
    
    // Find value of Z using calculated vectors in the plane
    double z = ( value_D - x*normal_plane.x - y*normal_plane.y ) / normal_plane.z;
    
    // Update zbuffer if the point is closer than current
    if( z > zbuffer.at( (unsigned int)abs(x) % 800 ).at( (unsigned int)abs(y) % 800 ) ||
        zbuffer.at( (unsigned int)abs(x) % 800 ).at( (unsigned int)abs(y) % 800 ) == DEPTH_INIT )
    {
        // Update the point if its closer, then render it in our polygon render function
        zbuffer.at( (unsigned int)abs(x) % 800 ).at( (unsigned int)abs(y) % 800 ) = z;
        return true;
    }

    // Depth of z at this point is not closer than current value in buffer. Do not render
    return false;
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
        //cout << "Not a real triangle!" << endl;
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
            double x0 = (y - y_max.y) / m_right + y_max.x;
            double x1 = (y - y_max.y) / m_left + y_max.x;
            
            // Swap x values if going in descending order
            if( x0 > x1 )
            {
                double temp = x0;
                x0 = x1;
                x1 = temp;
            }
            
            if( WIREFRAME_TOGGLE )
            {
                if( y == mid.y )
                {
                    // Draw horizontal line on triangle
                    for ( int x = (int) x0; x < x1; ++x )
                    {
                        if( update_zbuffer(x, y, number) )
                            renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
                }
                else
                {
                    // Draw only the outside edges
                    if( update_zbuffer((int)x0, y, number) )
                        renderPixel((int)x0, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    if( update_zbuffer((int)x1, y, number) )
                        renderPixel((int)x1, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                }
            }
            else
            {
                // Render pixels, starting from left slope of the triangle, to the right slope.
                for ( int x = (int) x0; x < x1; ++x )
                {
                    if( update_zbuffer(x, y, number) )
                        renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                }
            }
        }
        
        
        // Triangle two ( using pnew, pmid, pmin )
        m_left = ( p_new.y - y_min.y ) / ( p_new.x - y_min.x );
        m_right = ( mid.y - y_min.y ) / ( mid.x - y_min.x );
        for ( int y = (int)y_min.y; y < mid.y; ++y)
        {
            // x = (1/m)*(y - b ) + x0
            double x0 = (y - y_min.y) / m_left + y_min.x;
            double x1 = (y - y_min.y) / m_right + y_min.x;

            // Swap x values if going in descending order
            if( x0 > x1 )
            {
                double temp = x0;
                x0 = x1;
                x1 = temp;
            }
            
            if( WIREFRAME_TOGGLE )
            {
                if( y == y_min.y )
                {
                    // Draw horizontal line on triangle
                    for ( int x = (int) x0; x < x1; ++x )
                    {
                        if( update_zbuffer(x, y, number) )
                            renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
                }
                else
                {
                    // Draw only the outside edges
                    if( update_zbuffer(x0, y, number) )
                        renderPixel((int)x0, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    if( update_zbuffer(x1, y, number) )
                        renderPixel((int)x1, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                }
            }
            else
            {
                // Render pixels
                for ( int x = (int) x0; x < x1; ++x )
                {
                    if( update_zbuffer(x, y, number) )
                        renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                }
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
                double x0 = (y - y_min.y) / m_left + y_min.x;
                double x1 = (y - y_min.y) / m_right + y_min.x;
                
                // Swap x values if going in descending order
                if( x0 > x1 )
                {
                    double temp = x0;
                    x0 = x1;
                    x1 = temp;
                }
                
                if( WIREFRAME_TOGGLE )
                {
                    if( y == y_min.y )
                    {
                        // Draw horizontal line on triangle
                        for ( int x = (int) x0; x < x1; ++x )
                        {
                            if( update_zbuffer(x, y, number) )
                                renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                        }
                    }
                    else
                    {
                        // Draw only the outside edges
                        if( update_zbuffer(x0, y, number) )
                            renderPixel((int)x0, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                        if( update_zbuffer(x1, y, number) )
                            renderPixel((int)x1, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
                }
                else
                {
                    // Render pixels
                    for ( int x = (int) x0; x < x1; ++x )
                    {
                        if( update_zbuffer(x, y, number) )
                            renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
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
                double x0 = (y - y_max.y) / m_right + y_max.x;
                double x1 = (y - y_max.y) / m_left + y_max.x;
                
                // Swap x values if going in descending order
                if( x0 > x1 )
                {
                    double temp = x0;
                    x0 = x1;
                    x1 = temp;
                }
                
                if( WIREFRAME_TOGGLE )
                {
                    if( y == y_min.y )
                    {
                        // Draw horizontal line on triangle
                        for ( int x = (int) x0; x < x1; ++x )
                        {
                        if( update_zbuffer(x, y, number) )
                        renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                        }
                    }
                    else
                    {
                        // Draw only the outside edges
                        if( update_zbuffer(x0, y, number) )
                            renderPixel((int)x0, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                        if( update_zbuffer(x1, y, number) )
                            renderPixel((int)x1, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
                }
                else
                {
                    // Render pixels, starting from left slope of the triangle, to the right slope.
                    for ( int x = (int) x0; x < x1; ++x )
                    {
                        if( update_zbuffer(x, y, number) )
                            renderPixel(x, y, colors[number % 10].R, colors[number % 10].G, colors[number % 10].B);
                    }
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
