// Name: Alyssa Gaudioso  
// Quarter, Year: Fall, 2014
// Proj 1
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <GL/glut.h>
#include <cmath>
#include "point2d.h"
#include "point3d.h"
#include "triangle.h"
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

#define round(val) (int)(val+0.5)

using namespace std;

const int WINDOW_WIDTH = 900;
const int WINDOW_HEIGHT = 900;

string file_name;
ifstream myfile;
vector<Point2D> pts;
vector<Point3D> points;
vector<Triangle> triangles;
vector< vector<double> > zbuff(WINDOW_WIDTH, vector<double>(WINDOW_HEIGHT, -100000000) );
int clicks = 0;
float red[5] = {1.0,0.5,0.2,0.0,0.8};
float green[5] = {0.0, 0.4, 0.0, 0.8, 1.0};
float blue[5] = {0.2, 1.0, 0.4, 0.8, 0.0};
int i=0;

Point2D GLScreenToWindowCoordinates(const Point2D &p);
void GLKeyboardPress(unsigned char key, int x, int y);
void GLProcessArrow( int key, int x, int y);
void renderPixel(int x, int y, float r, float g, float b);
void GLMousePress(int button, int state, int x, int y);
bool update_z( Triangle tri, double x, double y );
void rotate( double q );
void GL_triangle();
void GLInit(int* argc, char** argv);
void parse_file();
int main(int argc, char** argv);


Point2D GLScreenToWindowCoordinates(const Point2D &p)
{
    GLdouble model[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    GLdouble projection[16];
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    Point2D converted;
    GLdouble temp;
    gluUnProject(p.x,viewport[3]-p.y, 0, model, projection, 
				 viewport, &converted.x, &converted.y, &temp);
    return converted;
}
void GLKeyboardPress(unsigned char key, int x, int y)
{
	if (key=='0')
	{
		pts.resize(0); 
		glutPostRedisplay();
	}
}
Point3D vert_avg()
{
	double x_total = 0;
	double y_total = 0;
	double z_total = 0;
	for( unsigned int k = 0; k < points.size(); k++ )
	{
		x_total += points[k].x;
		y_total += points[k].y;
		z_total += points[k].z;
	}
	
	double x_avg = (double)x_total/points.size();
	double y_avg = (double)y_total/points.size();
	double z_avg = (double)z_total/points.size();
	
	Point3D avg;
	avg.x = x_avg;
	avg.y = y_avg;
	avg.z = z_avg;
	return avg;		
}
void clear_buff()
{
	for(unsigned int xi = 0; xi < zbuff.size(); xi++)
	{
		for( unsigned int yi = 0; yi < zbuff[xi].size(); yi++ )
		{
			zbuff[xi][yi] = -100000000;
		}
	}
}
void pan( double q )
{
	Point3D temp;
	temp.x = q;
	temp.y = 0;
	temp.z = 0;
	for( unsigned int k = 0; k < points.size(); k++ )
	{
		points[k] += temp;
	}
}
void zoom( double q )
{
    Point3D avg = vert_avg();
    for( unsigned int k = 0; k < points.size(); k++ )//move to center axis
	{
		points[k] -= avg;
	}
    for( unsigned int k = 0; k < points.size(); k++ )
	{
		points[k] *= q;
	}
    for( unsigned int k = 0; k < points.size(); k++ )//move back to poistion 
	{
		points[k] += avg;
	}
}

double get_rad( double deg )
{
	return deg*0.0174532925;
}
	
void rotate( double q )
{
	Point3D avg = vert_avg();

	for( unsigned int k = 0; k < points.size(); k++ )//move to center axis
	{
		points[k] -= avg;
	}
	for( unsigned int k = 0; k < points.size(); k++ )//rotation
	{
		Point3D rotate;
		rotate.x = points[k].x;
		rotate.y = points[k].y*cos(get_rad(q)) - points[k].z*sin(get_rad(q));
		rotate.z = points[k].y*sin(get_rad(q)) + points[k].z*cos(get_rad(q));
		points[k] = rotate;
	}
	for( unsigned int k = 0; k < points.size(); k++ )//move back to poistion 
	{
		points[k] += avg;
	}
	
		
}
void GLProcessArrow( int key, int x, int y)
{
	if(key == GLUT_KEY_LEFT)//pan left
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		pan(-50);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();

	}
	if(key == GLUT_KEY_RIGHT)//pan right
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		pan(50);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();
	}
	if(key == GLUT_KEY_UP)//zoom out
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		zoom(1.2);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();
	}
	if(key == GLUT_KEY_DOWN)//zoom in 
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		zoom(0.5);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();
	}
	if(key == GLUT_KEY_PAGE_UP)//rotate up
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		rotate(-5);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();
	}
	else if(key == GLUT_KEY_PAGE_DOWN)//rotate down
	{
		pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
		clear_buff();
		rotate(5);		
		glutDisplayFunc(GL_triangle);
		glutPostRedisplay();
	}
	        
}
void renderPixel(int x, int y, float r = 1.0, float g = 1.0, float b = 1.0)
{
    if( x <= 900 && y <= 900 )
    {
        //~ cout << "render" << endl;
        //~ cout << x << " " << y << endl;
        glBegin(GL_POINTS);
        glVertex2f(x,y);
        glEnd();
    }
   
}
void GLMousePress(int button, int state, int x, int y)
{
    if (state==GLUT_DOWN)
    {
        pts.push_back(GLScreenToWindowCoordinates(Point2D(x,y)));
        glutPostRedisplay();
        glBegin(GL_POINTS);
        clicks++;
    }
}
bool update_z( Triangle tri, double x, double y )
{
	//first calculate plane of current triangle
	//do this by finding a normal vector, or the cross product of two vectors
	Point3D p1v = points[ tri.p1 ];//stores the values of each vertex in the triangle
	Point3D p2v = points[ tri.p2 ];
	Point3D p3v = points[ tri.p3 ];
	
	Point3D vector1 = p1v-p2v; //finds two vectors in the plane equal to the triangle
	Point3D vector2 = p1v-p3v;
	
	Point3D norm = vector1.cross_product( vector2 ); // the normal vector to the plane
	
	double D = norm.x*p1v.x + norm.y*p1v.y + norm.z*p1v.z;
	
	double z = ( D -norm.x*x - norm.y*y )/norm.z; //will find the z val of the passed in x&y
	
	//~ cout << "x: " << x << "    y: " << y << "    z: " << z << "    zbuff: " << zbuff[x][y];
	
	if( z > zbuff[x][y] || zbuff[x][y] == -100000000)//point is only or closest point
	{
		//~ cout << "    update buff" << endl;
		zbuff[x][y] = z; // we update the zbuff
		return true;// we want to render this point
	}
	else
	{
		//~ cout << "    dont update buff" << endl;
		return false; //there is already a closer point at this value
	}
}

void GL_triangle()
{	

	glClear(GL_COLOR_BUFFER_BIT);	
	if (!pts.empty())
	{
		renderPixel(pts[pts.size()-1].x, pts[pts.size()-1].y);
	}
		
	for( unsigned int k = 0; k < triangles.size(); k++ )
	{
		Point3D P0 = points[ triangles[k].p1 ]; // the three verts of the triangle
		Point3D P1 = points[ triangles[k].p2 ];
		Point3D P2 = points[ triangles[k].p3 ];
		
		//~ cout << "P0: ";
		//~ P0.print();
		//~ cout << "P1: ";
		//~ P1.print();
		//~ cout << "P2: ";
		//~ P2.print();
		
		Point3D Pymin; //in reference to the y values 
		Point3D Pymid;
		Point3D Pymax;

		double m_left; //slope of side connected to Pymax that is more left
		double m_right;//slope of side connected to Pymax that is more right 
		bool horz = true;

		if (i <5) //will change color for every five triangles
		{
			glColor3f(red[i], green[i], blue[i]);
			i++;	
		}	
		if (i>=5) //when reached 5 colors, start over with same 5 colors
		{
			i=0;
		}
		
		//next 6 ifs to test for 6 different cases
		//does not work is there is a horizontal line
		if (P0.y < P1.y && P1.y < P2.y)
		{
			Pymin=P0;
			Pymid=P1;
			Pymax=P2;
			horz = false;
		}
		if (P1.y < P0.y && P0.y < P2.y) 
		{
			Pymin=P1;
			Pymid=P0;
			Pymax=P2;
			horz = false;
		}
		if (P2.y < P0.y && P0.y < P1.y)
		{
			Pymin=P2;
			Pymid=P0;
			Pymax=P1;
			horz = false;
		}
		if (P1.y < P2.y && P2.y< P0.y)
		{
			Pymin=P1;
			Pymid=P2;
			Pymax=P0;
			horz = false;
		}
		if (P0.y < P2.y && P2.y < P1.y)
		{
			Pymin=P0;
			Pymid=P2;
			Pymax=P1;
			horz = false;
		}
		if (P2.y < P1.y && P1.y < P0.y)
		{
			Pymin=P2;
			Pymid=P1;
			Pymax=P0;
			horz = false;
		}

		//check to see if there is a horzontal line
		if( horz == true )
		{
			if( P0.y == P1.y && P0.y == P2.y )//if all 3 pnts on horz line ***i think this is pointless***
			{	
				double x0;
				double x1;
				if(P0.x < P1.x && P0.x < P2.x )//find xmin
					x0 = P0.x;
				if(P1.x < P0.x && P1.x < P2.x )
					x0 = P0.x;
				if(P2.x < P0.x && P2.x < P1.x )
					x0 = P0.x;
					
				if(P0.x > P1.x && P0.x > P2.x )//find ymin
					x1 = P0.x;
				if(P1.x > P0.x && P1.x > P2.x )
					x1 = P0.x;
				if(P2.x > P0.x && P2.x > P1.x )
					x1 = P0.x;
				
				for(int x=x0; x <= x1; x++)
				{
					if( update_z( triangles[k], x, P0.y ) == true )
					{
						renderPixel(x,P0.y);//*************************************************
					}
				}	
			}
			else //six different cases if exactly one horizontal line
			{
				bool top_horz = true; 
				if( P1.y == P2.y )
				{
					if( P0.y > P1.y )// horz on bottom
					{
						Pymin=P1;
						Pymid=P2;
						Pymax=P0;
						top_horz = false; 
					}
					else // horz on top
					{
						Pymin=P0;
						Pymid=P1;
						Pymax=P2;
					}
				}					
				if( P0.y == P2.y )
				{
					if( P1.y > P0.y )// horz on bottom
					{
						Pymin=P0;
						Pymid=P2;
						Pymax=P1;
						top_horz = false;
					}
					else // horz on top
					{
						Pymin=P1;
						Pymid=P0;
						Pymax=P2;
					}
				}
				if( P0.y == P1.y )
				{
					if( P2.y > P0.y )// horz on bottom
					{
						Pymin=P0;
						Pymid=P1;
						Pymax=P2;
						top_horz = false;
					}
					else // horz on top
					{
						Pymin=P2;
						Pymid=P0;
						Pymax=P1;
					}
				}
				//finish checking 6 cases
				if( top_horz == true )
				{
					//~ cout << "top horizontal" << endl;
					for(int y=Pymin.y; y <= Pymid.y; y++)
					{
						m_left=((Pymax.y-Pymin.y)/(Pymax.x-Pymin.x));
						m_right=((Pymid.y-Pymin.y)/(Pymid.x-Pymin.x));
						
						double x0 = (y - Pymin.y) / m_left + Pymin.x;
						double x1 = (y - Pymid.y) / m_right + Pymid.x;
						
						if(x0 > x1)
						{
							swap(x0, x1);
						}
						for(int x=x0; x <= x1; x++)
						{
							if( update_z( triangles[k], x, y ) == true )
							{
								renderPixel(x,y);//********************************************
							}
						}
					}					
				}
				else// horz on bottom
				{
					for(int y=Pymid.y; y <= Pymax.y; y++)
					{
						m_left=((Pymax.y-Pymin.y)/(Pymax.x-Pymin.x));
						m_right=((Pymax.y-Pymid.y)/(Pymax.x-Pymid.x));
						
						double x0 = (y - Pymax.y) / m_left + Pymax.x;
						double x1 = (y - Pymax.y) / m_right + Pymax.x;
						
						if(x0 > x1)//if right is not right of left, switch
						{
							swap(x0, x1);
						}

						for(int x=x0; x <= x1; x++)//scan line
						{
							if( update_z( triangles[k], x, y ) == true )
							{
								renderPixel(x,y);//******************************************
							}
						}
					}
				}
				
			}
		}
		//end check for case
		else
		{
			for(int y=Pymid.y; y <= Pymax.y; y++)
			{
				m_left=((Pymax.y-Pymin.y)/(Pymax.x-Pymin.x));
				m_right=((Pymax.y-Pymid.y)/(Pymax.x-Pymid.x));
				
				double x0 = (y - Pymax.y) / m_left + Pymax.x;
				double x1 = (y - Pymax.y) / m_right + Pymax.x;
				
				if(x0 > x1)//if right is not right of left, switch
				{
					swap(x0, x1);
				}
				for(int x=x0; x <= x1; x++)//scan line
				{
					if( update_z( triangles[k], x, y ) == true )
					{
						renderPixel(x,y);//*****************************************************
					}
				}
			}

			for(int y=Pymin.y; y <= Pymid.y; y++)
			{
				m_left=((Pymax.y-Pymin.y)/(Pymax.x-Pymin.x));
				m_right=((Pymid.y-Pymin.y)/(Pymid.x-Pymin.x));
				
				double x0 = (y - Pymin.y) / m_left + Pymin.x;
				double x1 = (y - Pymid.y) / m_right + Pymid.x;
				
				if(x0 > x1)
				{
					swap(x0, x1);
				}
				for(int x=x0; x <= x1; x++)
				{
					if( update_z( triangles[k], x, y ) == true )
					{
						renderPixel(x,y);//*****************************************************
					}
				}
			}
		}
		pts.resize(0);

	}

    glutSwapBuffers();
}

void GLInit(int* argc, char** argv)
{
    glutInit(argc, argv);//
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);//
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);//
    glutCreateWindow("CS 130 - Alyssa Gaudioso");//
    glMatrixMode(GL_PROJECTION_MATRIX);//
    glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);//
    glutDisplayFunc(GL_triangle);
    glutKeyboardFunc(GLKeyboardPress);//
    glutSpecialFunc(GLProcessArrow);
}
void parse_file()
{
    points.clear();
    triangles.clear();
     
    int num_pts;
    int num_tri;
    
    myfile.clear();
    myfile.open( file_name.c_str() );

    myfile >> num_pts;
    myfile >> num_tri;
    //~ cout << num_pts << " " << num_tri << endl;
    
    for( int i = 0; i < num_pts; i++ )
    {	
		double tmp1, tmp2, tmp3;
		Point3D tmp;
		
		myfile >> tmp1;
		myfile >> tmp2;
		myfile >> tmp3;
		
		//~ cout << tmp1 << " " << tmp2 << " " << tmp3 << endl;
		tmp.x = tmp1;
		tmp.y = tmp2;
		tmp.z = tmp3;
		points.push_back(tmp);
	}
	
    for( int i = 0; i < num_tri; i++ )
    {	
		double tmp1, tmp2, tmp3;
		Triangle tmp;
		
		myfile >> tmp1;
		myfile >> tmp2;
		myfile >> tmp3;
		
		//~ cout << tmp1 << " " << tmp2 << " " << tmp3 << endl;
		tmp.p1 = tmp1;
		tmp.p2 = tmp2;
		tmp.p3 = tmp3;
		triangles.push_back(tmp);
	}
	myfile.close();
	
}
int main(int argc, char** argv)
{   
    if( argc != 2 )
    {	cerr << "Usage: File name" << endl;
		return -1;
	}
	file_name = argv[1];
	parse_file(); 
    cout << "Pan: left arrow/right arrow" << endl
         << "Zoom: up arrow/down arrow" << endl
         << "Rotate: page up/page down" << endl;
	
	GLInit(&argc, argv);
	glutSpecialFunc(GLProcessArrow);
    glutMainLoop();
    
	return 0;
}

