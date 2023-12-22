#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
double pi = 3.14159265;
int mode = 1;

/**
 * Class that holds all the points in the polygon. Points are in the x,y,z coordinate system.
 */
class vertex
{
	public:
	double x,y,z,odd;
	vertex(double s=0,double t=0,double u=0,double o=0)
	:x(s),y(t),z(u),odd(o)
	{}
	friend class Polygon_container;
	vertex operator=(vertex& lhs)
    {
      x=lhs.x;
      y=lhs.y;
      z=lhs.z;
      return *this;
    }
};

//overloaded operators for vertex
bool operator==(vertex& rhs, vertex& lhs)
{
  return rhs.x==lhs.x && rhs.y==lhs.y && rhs.z==lhs.z;
};

bool operator!=(vertex& rhs, vertex& lhs)
{
  return !(rhs.x==lhs.x && rhs.y==lhs.y && rhs.z==lhs.z);
};

vertex operator *=(vertex& rhs, double i)
{
	rhs.x*=i;
	rhs.y*=i;
	rhs.z*=i;
	return rhs;
}

vertex operator * ( double i, vertex& rhs)
{
	return vertex(i*rhs.x,i*rhs.y,i*rhs.z);
}

vertex operator += (vertex& rhs, vertex& lhs)
{
	rhs.x+=lhs.x;
	rhs.y+=lhs.y;
	rhs.z+=lhs.z;
	return rhs;
}
vertex operator +(vertex& rhs, vertex& lhs)
{
	return vertex(rhs.x+lhs.x,rhs.y+lhs.y,rhs.z+lhs.z);
}

/**
 * Class that holds all the vertices of a given face inside the polgon.
 * Face is comprised of 3 vertices as polygon is subdivided into triangles.
 */
class face
{
	public:
	vertex one;
	vertex two;
	vertex three;
	face(vertex a, vertex b, vertex c)
	:one(a),two(b),three(c)
	{}
	bool containsVertex(vertex v)
	{
	  if (v==one or v==two or v==three)
	  {
		  return true;
	  }
	  return false;
	}
};

//overloaded operators for face
bool operator==(face& rhs, face& lhs)
{
  if (rhs.one != lhs.one){return false;}
  else if (rhs.one != lhs.two){return false;}
  else if (rhs.one != lhs.three){return false;}
  else if (rhs.two != lhs.one){return false;}
  else if (rhs.two != lhs.two){return false;}
  else if (rhs.two != lhs.three){return false;}
  else if (rhs.three != lhs.one){return false;}
  else if (rhs.three != lhs.two){return false;}
  else if (rhs.three != lhs.three){return false;}
  else{return true;}
};


//comparator for map objects used in polygon container
struct compare_v
{
   bool operator() (const vertex& v, const vertex& w)
   {
          if (v.x==w.x)
          {
                  if (v.y == w.y)
                  return v.z<w.z;
                  else
                  return v.y < w.y;
          }
          return v.x < w.x;
   }
};

/**
 * Class that stores and contains all the data used for the Loop subdivision.
 * Container starts with a set of vertices that are then averaged after subdividing.
 * New faces are also generated off of the original faces after Loop subdivision has
 * been applied.
 */
class Polygon_container
{
	public:
	vector<vertex> verts;
	vector<vertex> scaled_verts;
	vector<face> faces;
	vector <face> new_faces;
	vector<face> scaled_new_faces;
	map< int,vector<int> > face_rel;
    map<vertex,vector<vertex>,compare_v> neighbors;
	Polygon_container()
	{}
	bool isNeighbor(vertex home, vertex neighbor)
	{
		if (neighbors[home].size()==0)
		return false;
		for (int i(0); i < neighbors[home].size(); ++i)
	    {
		   if (neighbors[home][i]==neighbor)
		   return true;
		}
		return false;
	}

	//finds a vertex in the verts vector
	int whereIs(vertex v)
	{
		int counter(0);
		for(; counter < verts.size(); ++counter)
		{
			if (verts[counter]==v)
			{
				return counter;
			}
		}
	}

    //checks whether a vertex a is in vector of vertices
	bool isIn(vertex a, vector<vertex>& v)
	{
		for (int i(0); i < v.size(); ++i)
		{
		   if (v[i]==a)
		   {
		     return true;
		   }
		}
		return false;
	}

	//generates the neighbors for the polygon at its initial state
	void generate_neighbors()
	{
	  for (int i(0); i < faces.size(); i++)
	  {
		  for (int j(0); j < 3; ++j)
		  {
               if (j==0)
               {
				   if(!isNeighbor(faces[i].one,faces[i].two))
				   {
				     neighbors[faces[i].one].push_back(faces[i].two);
				     neighbors[faces[i].two].push_back(faces[i].one);
				   }
				   if (!isNeighbor(faces[i].one,faces[i].three))
				   {
				     neighbors[faces[i].one].push_back(faces[i].three);
				     neighbors[faces[i].three].push_back(faces[i].one);
				   }
			   }
			   else if (j==1)
			   {
			     if(!isNeighbor(faces[i].two,faces[i].one))
				 {
				    neighbors[faces[i].two].push_back(faces[i].one);
				    neighbors[faces[i].one].push_back(faces[i].two);
				 }
				 if (!isNeighbor(faces[i].two,faces[i].three))
				 {
				    neighbors[faces[i].two].push_back(faces[i].three);
				    neighbors[faces[i].three].push_back(faces[i].two);
				 }
			   }
			   else
			   {
			     if(!isNeighbor(faces[i].three,faces[i].one))
				 {
				    neighbors[faces[i].three].push_back(faces[i].one);
				    neighbors[faces[i].one].push_back(faces[i].three);
				 }
				 if (!isNeighbor(faces[i].three,faces[i].two))
				 {
				   neighbors[faces[i].three].push_back(faces[i].two);
				   neighbors[faces[i].two].push_back(faces[i].three);
				 }
			   }
		  }
	  }
	}

	//generates the vertices neighbor relationship amonst the other vertices in the
	//polygon
	void generate_neighbor_rel()
	{
	  for (int i(0); i < verts.size();++i)
	  {
		 for(int j(0); j < neighbors[verts[i]].size(); ++j)
		 {
			 int ind = whereIs(neighbors[verts[i]][j]);
			 face_rel[i].push_back(ind);
		 }
	  }
    }

    //populates the scaled faces vector and is called after Loop subdivision is done
    void generate_scaled_faces()
    {
		for (int i(0); i < new_faces.size(); ++i)
		{
			face f = new_faces[i];
			int one = whereIs(new_faces[i].one);
			int two = whereIs(new_faces[i].two);
			int three = whereIs(new_faces[i].three);
			vertex a = scaled_verts[one];
			vertex b = scaled_verts[two];
			vertex c = scaled_verts[three];
			face scale(a,b,c);
			scaled_new_faces.push_back(scale);
		}

	}

    //Loops subdivision: divides each face into 4 triangular faces and sets up the process for
    //the smoothing/averaging part to happen
	void subdivide()
	{
		for (int i(0); i < faces.size(); ++i)
		{
		   	vertex a = faces[i].one;
		   	vertex b = faces[i].two;
		   	vertex c = faces[i].three;
		   	vertex mid1((a.x+b.x)/2,(a.y+b.y)/2,(a.z+b.z)/2,1);
		   	vertex mid2((a.x+c.x)/2,(a.y+c.y)/2,(a.z+c.z)/2,1);
		   	vertex mid3((c.x+b.x)/2,(c.y+b.y)/2,(c.z+b.z)/2,1);
		   	if (!isIn(mid1,verts))
		   	{
		   	  verts.push_back(mid1);
		   	}
		   	if (!isIn(mid2,verts))
		   	{
		   	  verts.push_back(mid2);
		   	}
		   	if (!isIn(mid3,verts))
		   	{
		   	  verts.push_back(mid3);
		   	}
		   	face face1(a,mid1,mid2);
		   	face face2(b,mid1,mid3);
		   	face face3(c,mid2,mid3);
		   	face face4(mid1,mid2,mid3);
		   	new_faces.push_back(face1);
		   	new_faces.push_back(face2);
		   	new_faces.push_back(face3);
		   	new_faces.push_back(face4);
		   	if (!isIn(mid1,neighbors[a]))
		   	{
				neighbors[a].push_back(mid1);
				neighbors[mid1].push_back(a);
			}
			if (!isIn(mid2,neighbors[a]))
		   	{

				neighbors[a].push_back(mid2);
				neighbors[mid2].push_back(a);
			}
			if (!isIn(mid1,neighbors[b]))
		   	{
				neighbors[b].push_back(mid1);
				neighbors[mid1].push_back(b);
			}
			if (!isIn(mid3,neighbors[b]))
		   	{
				neighbors[b].push_back(mid3);
				neighbors[mid3].push_back(b);
			}
			if (!isIn(mid2,neighbors[c]))
		   	{
				neighbors[c].push_back(mid2);
				neighbors[mid2].push_back(c);
			}
			if (!isIn(mid3,neighbors[c]))
		   	{
				neighbors[c].push_back(mid3);
				neighbors[mid3].push_back(c);
			}

			if (!isIn(mid1,neighbors[mid2]))
			{
				neighbors[mid2].push_back(mid1);
			}
			if (!isIn(mid1,neighbors[mid3]))
			{
				neighbors[mid3].push_back(mid1);
			}
			if (!isIn(mid2,neighbors[mid1]))
			{
				neighbors[mid1].push_back(mid2);
			}
			if (!isIn(mid2,neighbors[mid3]))
			{
				neighbors[mid3].push_back(mid2);
			}
			if (!isIn(mid3,neighbors[mid1]))
			{
				neighbors[mid1].push_back(mid3);
			}
			if (!isIn(mid3,neighbors[mid2]))
			{
				neighbors[mid2].push_back(mid3);
			}

		}
		generate_neighbor_rel();
	}

	//averages the points in the polygon and generates the scaled faces vector
	void average()
	{
	  for (int i(0); i < verts.size(); ++i)
	  {
		  if (verts[i].odd == 1)
		  {
			  verts[i].odd=0;
			  scaled_verts.push_back(verts[i]);
		  }
		  else
		  {

			  double k = neighbors[verts[i]].size();
			  double c= cos(2*pi/k);
			  c*=(1.0/4.0);
			  c+=(3.0/8.0);
			  double d = pow(c,2.0);
			  double B = (1.0/k)*((5.0/8.0)-d );
			  vertex scale = (1-(k*B))*verts[i];
			  for (int j(0); j < neighbors[verts[i]].size(); ++j)
			  {
				  vertex add = neighbors[verts[i]][j];
				  add*=B;
				  scale+=add;
			  }
			  scaled_verts.push_back(scale);
		  }
	  }
	  generate_scaled_faces();
    }

   //used for switching between iterations so that you know when you subdivide the
   //difference between an even vertex and an odd vertex
   void all_even()
   {
	 for (int i(0); i > scaled_verts.size(); ++i)
	 {
	    scaled_verts[i].odd=0;
	 }
	 for ( int j(0); j <  scaled_new_faces.size(); ++j)
	 {
	   	 scaled_new_faces[j].one.odd=0;
	   	 scaled_new_faces[j].two.odd=0;
	   	 scaled_new_faces[j].three.odd=0;
	 }
   }

   //prepares the next iteration's polygon container with the verts and faces that it will use to
   //do Loop's subdivision on
   void transfer_info(Polygon_container& p)
   {
	   	all_even();
	    p.verts.clear();
	    p.faces.clear();
	    copy(scaled_verts.begin(),scaled_verts.end(),back_inserter(p.verts));
	    copy(scaled_new_faces.begin(),scaled_new_faces.end(),back_inserter(p.faces));
   }
};

Polygon_container zero_it;
Polygon_container first_it;
Polygon_container sec_it;
Polygon_container third_it;
// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(double x, double y, float r = 1.0, float g = 1.0, float b = 1.0)
{
glBegin(GL_POINTS);

glVertex2f(x,y);

glEnd();
}

void renderLine(double x0, double y0, double x1, double y1)
{
    glBegin(GL_POINTS);
    if (x0==x1)
    {
		double t = y0;
		if (y0>y1)
		{
		  y0=y1;
		  y1=t;
		}
		for(;y0<=y1;++y0)
		renderPixel(x0,y0);
	}
	else
	{
	  if (x0 > x1)
	  {
		  double temp = x0;
		  x0=x1;
		  x1=temp;
		  double t = y0;
		  t=y0;
		  y0=y1;
		  y1=t;
	  }
      double i = round(x0);
      double m = (y1-y0)/(x1-x0);
      double y = y0+m*(i-x0);
      while(i<=x1)
      {
		 renderPixel(i,y);
		 ++i;
		 y+=m;
      }
     }
    glEnd();
}

//gets the info of the mesh from a text file and stores it in the 1st iterations
//polygon container
void get_info(char** argv, Polygon_container& p)
{
	ifstream read;
	read.open(argv[1]);
	int vert,faces;
	read >> vert >> faces;
	for (int i(0); i < vert; ++i)
	{
	  vertex val;
	  read >> val.x >> val.y >> val.z;
	  val.odd=0;
	  p.verts.push_back(val);
	}
	for (int j(0); j < faces; ++j)
	{
	  int v1,v2,v3;
	  read >> v1 >> v2 >> v3;
	  face fac(p.verts[v1],p.verts[v2],p.verts[v3]);
	  p.faces.push_back(fac);
	}
	p.generate_neighbors();
    read.close();
}

//front_view(0-3) displays the fron front of the polygon at iterations 0-3
void front_view(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+10,4*(p.scaled_verts[i].y)+10,4*(p.scaled_verts[ind].x)+10,4*(p.scaled_verts[ind].y)+10);
		}
	}
	glEnd();
}

void front_view1(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+160,4*(p.scaled_verts[i].y)+160,4*(p.scaled_verts[ind].x)+160,4*(p.scaled_verts[ind].y)+160);
		}
	}
	glEnd();
}

void front_view2(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+310,4*(p.scaled_verts[i].y)+310,4*(p.scaled_verts[ind].x)+310,4*(p.scaled_verts[ind].y)+310);
		}
	}
	glEnd();
}

void front_view3(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+460,4*(p.scaled_verts[i].y)+460,4*(p.scaled_verts[ind].x)+460,4*(p.scaled_verts[ind].y)+460);
		}
	}
	glEnd();
}

//side_view(0-3) displays the fron side of the polygon at iterations 0-3
void side_view(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].z)+10,4*(p.scaled_verts[i].y)+10,4*(p.scaled_verts[ind].z)+10,4*(p.scaled_verts[ind].y)+10);
		}
	}
	glEnd();
}

void side_view1(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].z)+160,4*(p.scaled_verts[i].y)+160,4*(p.scaled_verts[ind].z)+160,4*(p.scaled_verts[ind].y)+160);
		}
	}
	glEnd();
}

void side_view2(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].z)+310,4*(p.scaled_verts[i].y)+310,4*(p.scaled_verts[ind].z)+310,4*(p.scaled_verts[ind].y)+310);
		}
	}
	glEnd();
}
void side_view3(Polygon_container& p)
{
	glBegin(GL_POINTS);
	for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].z)+460,4*(p.scaled_verts[i].y)+460,4*(p.scaled_verts[ind].z)+460,4*(p.scaled_verts[ind].y)+460);
		}
	}
	glEnd();
}

//top_view(0-3) displays the fron top of the polygon at iterations 0-3
void top_view(Polygon_container& p)
{
	glBegin(GL_POINTS);
    for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+10,4*(p.scaled_verts[i].z)+10,4*(p.scaled_verts[ind].x)+10,4*(p.scaled_verts[ind].z)+10);
		}
	}
	glEnd();
}
void top_view1(Polygon_container& p)
{
	glBegin(GL_POINTS);
    for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+160,4*(p.scaled_verts[i].z)+160,4*(p.scaled_verts[ind].x)+160,4*(p.scaled_verts[ind].z)+160);
		}
	}
	glEnd();
}
void top_view2(Polygon_container& p)
{
	glBegin(GL_POINTS);
    for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+310,4*(p.scaled_verts[i].z)+310,4*(p.scaled_verts[ind].x)+310,4*(p.scaled_verts[ind].z)+310);
		}
	}
	glEnd();
}
void top_view3(Polygon_container& p)
{
	glBegin(GL_POINTS);
    for(int i(0); i < p.scaled_verts.size(); ++i)
	{
		for (int j(0); j < p.face_rel[i].size(); ++j)
		{
		  int ind = p.face_rel[i][j];
		  renderLine(4*(p.scaled_verts[i].x)+460,4*(p.scaled_verts[i].z)+460,4*(p.scaled_verts[ind].x)+460,4*(p.scaled_verts[ind].z)+460);
		}
	}
	glEnd();
}

//gets Spacebar to be used to go through face views of the polygon's iterations of Loop subdivision
void specialKey( unsigned char key, int x, int y )
{

       if (key == 32)
       {
		   if(mode==3)
           {
			 mode=1;
		   }
		   else
		   {
			 ++mode;
		   }
	   }
	   glutPostRedisplay();
	   return;
}

void renderCircle(int x, int y, int r)
{
   glBegin(GL_POINTS);
   double i = 0;
   double j = round(r);
   while (i <= j)
   {
	 renderPixel(int(i)+x,int(j)+y);
	 renderPixel(-(int(i))+x,-(int(j))+y);
	 renderPixel(-(int(i))+x,int(j)+y);
	 renderPixel(int(i)+x,-(int(j))+y);
	 renderPixel(int(i)+x,int(j)+y);
	 renderPixel(-(int(j))+x,-(int(i))+y);
	 renderPixel(-(int(j))+x,int(i)+y);
	 renderPixel(int(j)+x,-(int(i))+y);
	 renderPixel(int(j)+x,int(i)+y);
	 ++i;
	 double d = pow(double(r),2.0)-pow(i,2.0)-pow(j,2.0);
	 if(d<0)
	 {
	   --j;
	 }
   }
   glEnd();
}

void GL_render()
{

   glClear(GL_COLOR_BUFFER_BIT);
   if (mode == 0)
   {
   }
   if (mode == 1)
   {
     front_view(zero_it);
     front_view1(first_it);
     front_view2(sec_it);
     front_view3(third_it);
   }
   else if ( mode == 2)
   {
     side_view(zero_it);
     side_view1(first_it);
     side_view2(sec_it);
     side_view3(third_it);
   }
   else
   {
     top_view(zero_it);
     top_view1(first_it);
     top_view2(sec_it);
     top_view3(third_it);
   }
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
	glutCreateWindow("CS 130 - Tyler Wilson");

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// This is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
    glutDisplayFunc(GL_render);
    glutKeyboardFunc(specialKey);
}

int main(int argc, char** argv)
{

	get_info(argv, zero_it);
	copy(zero_it.verts.begin(),zero_it.verts.end(),back_inserter(zero_it.scaled_verts));
     copy(zero_it.faces.begin(),zero_it.faces.end(),back_inserter(zero_it.scaled_new_faces));
     zero_it.generate_neighbor_rel();
     zero_it.transfer_info(first_it);
     first_it.subdivide();
     first_it.average();
     first_it.transfer_info(sec_it);
     sec_it.subdivide();
     sec_it.average();
     sec_it.transfer_info(third_it);
     third_it.subdivide();
     third_it.average();
	GLInit(&argc, argv);
	glutMainLoop();

	return 0;
}
