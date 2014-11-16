// Name: Jordon Ng 860977164
// Quarter, Year: Fall 2014
// Lab: 3
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include "point2d.h"
#include <vector>
#include <cmath> //for sin and cos
#include <fstream> //used to open and close files
#include <sstream> //used to parse inputs
#include <iostream>
#include <algorithm> //for fill n

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const float PI = 3.14159265;

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(float x, float y, float r = 1.0, float g = 1.0, float b = 1.0)
{
	glColor3f(r,g,b); //add in colors if specified
	glBegin(GL_POINTS);
	glVertex2f(x,y);
	glEnd();

}

struct Color_RGB //struct used to hold color data
{
	float r;
	float g;
	float b;
        Color_RGB()
                :r(1.0),g(1.0),b(1.0) {}
	Color_RGB(float x, float y, float z)
		:r(x),g(y),b(z) {}     
};

struct my_vertices
{
        int v1;
        int v2;
        int v3;
        my_vertices(int x, int y, int z)
                :v1(x), v2(y), v3(z) {}
};

struct Vertex
{
	float x;
	float y;
	float z;

	Vertex() : x(0.0), y(0.0), z(0.0)
	{}

	Vertex(float x, float y, float z)
		: x(x), y(y), z(z)
	{}

	Vertex(const Vertex & v)
		: x(v.x), y(v.y), z(v.z)
	{}

	Vertex operator+(const Vertex & rhs) const
	{ return Vertex(x + rhs.x, y + rhs.y, z + rhs.z); }
	Vertex operator-(const Vertex & rhs) const
	{ return Vertex(x - rhs.x, y - rhs.y, z - rhs.z); }
	Vertex operator*(float rhs) const
	{ return Vertex(x * rhs, y * rhs, z * rhs); }
	Vertex operator/(float rhs) const
	{ return Vertex(x / rhs, y / rhs, z / rhs); }
	Vertex operator+=(const Vertex & rhs)
	{ x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Vertex operator-=(const Vertex & rhs)
	{ x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Vertex operator*=(float rhs)
	{ x *= rhs; y *= rhs; z *= rhs; return *this; }
	Vertex operator/=(float rhs)
	{ x /= rhs; y /= rhs; z /= rhs; return *this; }

	float magnitude() const
	{ return sqrt(x * x + y * y + z * z); }
	void normalize()
	{ *this /= magnitude(); }
	float dot(const Vertex & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}
	Vertex cross(const Vertex& lhs, const Vertex & rhs) const
	{
		return Vertex(lhs.y * rhs.z - lhs.z * rhs.y,
		              lhs.z * rhs.x - lhs.x * rhs.z,
			      lhs.x * rhs.y - lhs.y * rhs.x);
	}
};

struct my_tri
{
	Vertex vertex_one;
	Vertex vertex_two;
	Vertex vertex_three;
	Color_RGB tri_color;
	my_tri(Vertex one, Vertex two, Vertex three, Color_RGB c)
		:vertex_one(one), vertex_two(two), vertex_three(three), tri_color(c) {}
};

//--------------------------------GLOBAL DATA CONTAINERS------------------------
vector<Point2D> pts;
vector<Vertex> poly_pts;
vector<my_vertices> poly_vertices;
float zbuf[WINDOW_WIDTH][WINDOW_HEIGHT];
Color_RGB frame_buf[WINDOW_WIDTH][WINDOW_HEIGHT]; //holds color information of each point
vector<my_tri> triangles; //holds number of triangles each poly has
Color_RGB output_colors[7] = {Color_RGB(1.0,0.0,0.0),Color_RGB(0.8,0.4,0.0), Color_RGB(1.0,1.0,0.0),
			      Color_RGB(0.0,1.0,0.0),Color_RGB(0.2,0.4,1.0), Color_RGB(0.4,0.0,1.0), Color_RGB(1.0,0.0,0.8)};
float angle_rotate = 0.0;
//--------------------------------GLOBAL DATA CONTAINERS------------------------


Point2D GLscreenToWindow(const Point2D &p)
{
	GLdouble model[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,model);
	GLdouble projection[16];
	glGetDoublev(GL_PROJECTION_MATRIX,projection);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	Point2D converted;
	GLdouble temp;
	gluUnProject(p.x,viewport[3]- p.y,0,model,projection,viewport,&converted.x,&converted.y,&temp);
	return converted;
}

//button = what button, state = up or down, (x,y) is location of the mouse
void GlMousePress(int button, int state, int x, int y)
{
	int get_butt;
	switch(button)
	{
		case GLUT_LEFT_BUTTON: get_butt = 0; break;
		case GLUT_MIDDLE_BUTTON: get_butt = 1; break;
		case GLUT_RIGHT_BUTTON: get_butt = 2; break;
	}
	if(state == GLUT_DOWN && get_butt == 0) //left button is pressed
	{
		pts.push_back(GLscreenToWindow(Point2D(x,y)));
	}

	glutPostRedisplay();
}

void swap(float& x, float& y)
{
	float temp = x;
	x = y;
	y = temp;
}

void FindEdge(Vertex P0, Vertex P1, int* left_edge, int* right_edge)
{
	float x, mx;

	if(P0.y > P1.y) //get y_max
    	{
		swap(P0.x,P1.x);
		swap(P0.y,P1.y);
  	}

	mx = (P0.y == P1.y) ? (P1.x - P0.x) : (P1.x - P0.x) / (P1.y - P0.y);
    	x = P0.x;

    	for(int i = int(P0.y); i <= (int)P1.y; i++) //scan through vertically
    	{
        	if(x < (float)left_edge[i]) left_edge[i] = (int)x;
        	if(x > (float)right_edge[i]) right_edge[i] = (int)x;
        	x += mx;
    	}
}

float get_z(const int& x, const int& y, const Vertex& one, const Vertex& two, const Vertex& three)
{
        Vertex first_sub = three - one;
        Vertex sec_sub = two - one;
        Vertex c_prod = c_prod.cross(first_sub, sec_sub); 

        float a_coeff = c_prod.x, b_coeff = c_prod.y, c_coeff = c_prod.z;
        float d = -((a_coeff * one.x) + (b_coeff * one.y) + (c_coeff * one.z));
        if(c_coeff == 0) return 0;
        return ((-a_coeff * x - b_coeff * y - d) / c_coeff); //return z  
}

float get_next_z(const int& x, const int& y, const Vertex& one, const Vertex& two, const Vertex& three)
{
        Vertex first_sub = three - one;

        Vertex sec_sub = two - one;
        Vertex c_prod = c_prod.cross(first_sub, sec_sub);
        
        float a_coeff = c_prod.x, b_coeff = c_prod.y, c_coeff = c_prod.z;
        float prev_z = get_z(x,y,one,two,three);
        if(c_coeff == 0) return 0;
        return prev_z - (a_coeff / c_coeff); //return next z (x+1)
}

Vertex get_avg_ver() //compute average vertex from every poly
{
        float avg_x = 0.0, avg_y = 0.0, avg_z = 0.0;
        int num_tri = triangles.size();

        for(unsigned i = 0; i < num_tri; ++i)
        {
                avg_x += triangles[i].vertex_one.x;
                avg_x += triangles[i].vertex_two.x;
                avg_x += triangles[i].vertex_three.x;

                avg_y += triangles[i].vertex_one.y;
                avg_y += triangles[i].vertex_two.y;
                avg_y += triangles[i].vertex_three.y;
                
                avg_z += triangles[i].vertex_one.z;
                avg_z += triangles[i].vertex_two.z;
                avg_z += triangles[i].vertex_three.z;
        }
        avg_x /= (float)num_tri * 3.0; 
        avg_y /= (float)num_tri * 3.0;
        avg_z /= (float)num_tri * 3.0;
        
        return Vertex(avg_x, avg_y, avg_z); 
}

//multiplies a 4x4 to a 4x4 matrix and puts into sol_mat
void in_n_out_four_by_four(float (&lhs_mat)[4][4], float (&rhs_mat)[4][4], float (&sol_mat)[4][4])
{
        for(int x = 0; x < 4; ++x)
        {
                for(int y = 0; y < 4; ++y)
                {
                        sol_mat[x][y] = 0.0;
                        for(int z = 0; z < 4; ++z)
                        {
                                sol_mat[x][y] += lhs_mat[x][z] * rhs_mat[z][y];
                        }
                }
        }               
}

//multiplies a 4x4 to a 4x1 matrix and puts into sol_mat
void in_n_out_four_by_one(float (&lhs_mat)[4][4], float (&rhs_mat)[4], float (&sol_mat)[4])
{
        for(int x = 0; x < 4; ++x)
        {
                sol_mat[x] = 0;
                for(int y = 0; y < 4; ++y)
                {
                        sol_mat[x] += lhs_mat[x][y] * rhs_mat[y];
                }
        }
}

float get_rad(float deg)
{
        return (deg * PI) / 180;
}

void Scan_Conv(const my_tri& t, float deg_rotate, Vertex avg)
{	
	int left_edge[800], right_edge[800]; //set edges to scan equal to window size
	for(unsigned i = 0; i < 800; ++i)
	{
		left_edge[i] = 800;
		right_edge[i] = 0;
	}

	
	Vertex one = t.vertex_one;
	Vertex two = t.vertex_two;
	Vertex three = t.vertex_three; 
        
        float trn_mat[4][4] = {{1,0,0,avg.x},{0,1,0,avg.y},{0,0,1,avg.z},{0,0,0,1}};
        float in_trn_mat[4][4] = {{1,0,0,-avg.x},{0,1,0,-avg.y},{0,0,1,-avg.z},{0,0,0,1}};
        float roty_mat[4][4] = {{cos(get_rad(deg_rotate)),0,-sin(get_rad(deg_rotate)),0}, {0,1,0,0},
                                {sin(get_rad(deg_rotate)),0,cos(get_rad(deg_rotate)),0},{0,0,0,1}};
        float m_mat1[4][4];
        float m_mat2[4][4];
        in_n_out_four_by_four(trn_mat, roty_mat, m_mat1);
        in_n_out_four_by_four(m_mat1, in_trn_mat, m_mat2);
        
        float update_mat[4] = {0,0,0,1};
        float old_vert_one[4] = {one.x, one.y, one.z, 1};
        float old_vert_two[4] = {two.x, two.y, two.z, 1};
        float old_vert_three[4] = {three.x, three.y, three.z, 1};
        
        float sin_t = sin(get_rad(deg_rotate));
        float cos_t = cos(get_rad(deg_rotate));

        in_n_out_four_by_one(m_mat2, old_vert_one, update_mat);      
        one.x = update_mat[0]; one.y = update_mat[1]; one.z = update_mat[2];

        in_n_out_four_by_one(m_mat2, old_vert_two, update_mat); 
        two.x = update_mat[0]; two.y = update_mat[1]; two.z = update_mat[2];


        in_n_out_four_by_one(m_mat2, old_vert_three, update_mat);  
        three.x = update_mat[0]; three.y = update_mat[1]; three.z = update_mat[2];     

	FindEdge(one,two,left_edge,right_edge); //find edges between all points in triangle
	FindEdge(two,three,left_edge,right_edge);
	FindEdge(three,one,left_edge,right_edge);       

        float curr_z = 0.0;
        float z_prime = 0.0;        
        
        for(int i = 0; i < 800; ++i)
                for(int j = 0; j < 800; ++j)
                        zbuf[i][j] = -1000;

	for(unsigned j = 0; j < WINDOW_HEIGHT; ++j) //for number of scan lines in the WINDOW y values
    	{
        	if(left_edge[j] <= right_edge[j])
		{
            		for(unsigned i = left_edge[j]; i < right_edge[j]; ++i) //horizontal scan x values
			{
                                if(i < 800 && i >= 0)
                                { 
                                      curr_z = get_z(i,j,one,two,three);
                                      //z_prime = get_next_z(i,j,one,two,three);
                           
                                      //populate z buffer
                                      if(curr_z > zbuf[i][j]) //update the z value 
                                      {                                    
                                              zbuf[i][j] = curr_z; //set new depth
                                              frame_buf[i][j] = t.tri_color; //set new color
                                      }
                                }
                                	
			}
			
		}
    	}
}

void GLKeyboardPress(unsigned char key, int x, int y) //check if 0 is pressed to clear screen
{
	if(key == '0')
	{
		pts.resize(0);
		triangles.clear(); //empty triangle container
		glutPostRedisplay();
	}
}

void GLProcessArrow(int key, int x, int y)
{
        if(key == GLUT_KEY_LEFT)
        {
                cout << "you pressed left\n";
                angle_rotate += 10.0; //rotate to the left
                for(unsigned i = 0; i < WINDOW_HEIGHT; ++i)
                        for(unsigned j = 0; j < WINDOW_WIDTH; ++j)
                                frame_buf[i][j] = Color_RGB();
                glutPostRedisplay();
        }
        else if(key == GLUT_KEY_RIGHT)
        {
                cout << "you pressed right\n";
                angle_rotate -= 10.0; //rotate to the right
                for(unsigned i = 0; i < WINDOW_HEIGHT; ++i)
                        for(unsigned j = 0; j < WINDOW_WIDTH; ++j)
                                frame_buf[i][j] = Color_RGB();
                glutPostRedisplay();
        }
}

void Pts_To_Tri(const vector<my_vertices>& V, const vector<Vertex>& P) //reads in from poly_vertices and connects poly_pts
{
	int num_ver = V.size();
        for(unsigned i = 0, k = 0; i < num_ver; ++i,++k)
        {
                if(k == 7) k = 0;
                triangles.push_back(my_tri(P[V[i].v1], P[V[i].v2], P[V[i].v3], output_colors[k]));
        }

}

void GL_render()
{
	glClear(GL_COLOR_BUFFER_BIT);
        
	Pts_To_Tri(poly_vertices,poly_pts); //convert every 3 pts to a triangle
        Vertex avg = get_avg_ver(); //center of the polygon, that needs to be translated  
	for(unsigned i = 0; i < triangles.size(); ++i) //print out all the triangles
	{
		Scan_Conv(triangles[i], angle_rotate, avg); //fill zbuf with appropriate values for each poly
	}
        for(unsigned i = 0; i < WINDOW_HEIGHT; ++i)
                for(unsigned j = 0; j < WINDOW_WIDTH; ++j)
                        renderPixel(i,j, frame_buf[i][j].r, frame_buf[i][j].g, frame_buf[i][j].b); //output the frame buffer

	glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	// ...
	// Complete this function
	// ...
	
	glutCreateWindow("CS 130 - Jordon Ng - Project 1: z-buffering");

	glutDisplayFunc(GL_render);
	glutKeyboardFunc(GLKeyboardPress);
	
	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// For the purposes of this lab, this is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);


	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
}

/* @param: pass in a filename to be read
 * @return: fills in the elements of three vectors: 1. poly_pts 2. poly_vertices 3. zbuffer array
*/
void ProcessFileInput(string fname)
{
        ifstream inbuf;
        inbuf.open(fname.c_str());
        int num_pts, num_tri;
        inbuf >> num_pts >> num_tri; //get input bounds
        
        for(; num_pts > 0 ; --num_pts)
        {
                int x_in, y_in, z_in;
                inbuf >> x_in >> y_in >> z_in;
                poly_pts.push_back(Vertex(x_in, y_in, z_in)); //put into vertex vector
        }
        for(; num_tri > 0; --num_tri)
        {
                int x_in, y_in, z_in;
                inbuf >> x_in >> y_in >> z_in;
                poly_vertices.push_back(my_vertices(x_in, y_in, z_in));
        }
        inbuf.close();
}

int main(int argc, char** argv)
{
        if(argc <= 1) return 0; //need to specify a file to read from
        string file_name = argv[1];

        ProcessFileInput(file_name);
        	
	GLInit(&argc, argv);
        glutSpecialFunc(GLProcessArrow);
	//glutMouseFunc(GlMousePress);
	glutMainLoop();

	return 0;
}
