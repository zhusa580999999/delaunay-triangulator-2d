/* 
=============================================================================== 
 
  FILE:  sp_viewer.cpp 
   
  CONTENTS: 
   
    This little tool visualizes a streaming point cloud. 
   
  PROGRAMMERS: 
   
    martin isenburg@cs.unc.edu 
   
  COPYRIGHT: 
   
    copyright (C) 2005  martin isenburg@cs.unc.edu 
     
    This software is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
   
  CHANGE HISTORY: 
   
    17 February 2005 -- created on the day that i wrote henna's birthday card 
   
=============================================================================== 
*/ 
 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 

#ifdef APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h> 
#endif /* APPLE */

#include "vec3fv.h" 
#include "vec3iv.h" 
 
#include "smreader_sma.h" 
#include "smreader_smb.h" 
#include "smreader_smc.h" 
#include "smreader_ply.h" 
#include "spconverter.h" 
 
#include "spreader_spa.h" 
#include "spreader_spb.h" 
#include "spreader_spc.h" 
#include "spreader_node.h" 
#include "spreader_raw_d.h" 
 
#include "spreadscattered.h" 
 
#include "sscontainer2d.h" 
#include "sscontainer3d.h" 
 
#include "spdelaunay2d.h" 
#include "spdelaunay3d.h" 
 
#include "smwriter_nil.h" 
#include "svwriter_nil.h" 
 
#include "positionquantizer_new.h" 
 
 
#include <ext/hash_map> 
 
// MOUSE INTERFACE 
 
int LeftButtonDown=0; 
int MiddleButtonDown=0; 
int RightButtonDown=0; 
int OldX,OldY,NewX,NewY; 
float Elevation=0; 
float Azimuth=0; 
float DistX=0; 
float DistY=0; 
float DistZ=2; 
 
// VISUALIZATION SETTINGS 
 
float boundingBoxMin[3]; 
float boundingBoxMax[3]; 
float boundingBoxScale=1.0f; 
float boundingBoxTranslateX = 0.0f; 
float boundingBoxTranslateY = 0.0f; 
float boundingBoxTranslateZ = 0.0f; 

// GLOBAL CONTROL VARIABLES 

//int WindowW=1024, WindowH=768; 
int WindowW=800, WindowH=600; 
int InteractionMode=0; 
int AnimationOn=0; 
int WorkingOn=0; 
int PrintToFileOn=0; 
 
// COLORS 
 
float colours_diffuse[10][4]; 
float colours_white[4]; 
float colours_light_blue[4]; 
 
// DATA STORAGE FOR STREAM VISUALIZER OUTPUT 
 
char* file_name = "points.spa"; 
char* file_name_header = 0; 
 
bool delaunay2 = false; 
bool delaunay3 = false; 
bool terrain = true; 
bool flatten = false; 
bool ispa = false; 
bool ispb = false; 
bool inode = false; 
bool isma = false; 
bool ismb = false; 
bool ismc = false; 
int scatter = 0; 
 
SMreader* smreader = 0; 
SPreader* spreader = 0; 
PositionQuantizerNew* pq = 0; 
 
SScontainer2D* sscontainer2d = 0; 
SScontainer3D* sscontainer3d = 0; 
 
SPdelaunay2D* spdelaunay2d = 0; 
SPdelaunay3D* spdelaunay3d = 0; 
 
int npoints; 
 
int EXACTLY_N_STEPS = 100; 
int EVERY_NTH_STEP = -1; 
int NEXT_STEP; 
 
int DIRTY_MESH=1; 
int REPLAY_IT=0; 
int REPLAY_COUNT=0; 
int GRID_PRECISION=8; 
int STREAM_COLORING = 0; 
int RENDER_MODE_UNFINALIZED = 1; 
int RENDER_MODE_ACTIVE = 1; 
int RENDER_MODE_INFINITE = 2; 
int POINT_SIZE = 5;
int RENDER_CIRCUMSPHIRCLES = 1; 
int RENDER_BOUNDINGBOX = 0; 
int LINE_WIDTH = 1; 
 
unsigned char* Framebuffer = 0; 
char* PrintFileName = "frame"; 
int Time=0; 
 
float* SpecialPoint = 0; 
int SpecialElement = -1; 
 
typedef struct GridVertex 
{ 
  float v[3]; 
  int number; 
} GridVertex; 
 
typedef struct RenderVertex 
{ 
  RenderVertex* buffer_next; 
  float v[3]; 
} RenderVertex; 
 
typedef __gnu_cxx::hash_map<int, int> my_grid_hash; 


static my_grid_hash* grid_hash; 
 
// efficient memory allocation 
 
static GridVertex* grid_vertex_buffer = 0; 
static int grid_vertex_buffer_size = 0; 
static int grid_vertex_buffer_alloc = 0; 

static void initGridVertexBuffer(int alloc) 
{ 
  if (grid_vertex_buffer) 
  { 
    if (grid_vertex_buffer_alloc < alloc) 
    { 
      grid_vertex_buffer_alloc = alloc; 
      free(grid_vertex_buffer); 
      grid_vertex_buffer = (GridVertex*)malloc(sizeof(GridVertex)*grid_vertex_buffer_alloc); 
    } 
  } 
  else 
  { 
    grid_vertex_buffer_alloc = alloc; 
    grid_vertex_buffer = (GridVertex*)malloc(sizeof(GridVertex)*grid_vertex_buffer_alloc); 
  } 
  grid_vertex_buffer_size = 0; 
} 
 
static int allocGridVertex() 
{ 
  if (grid_vertex_buffer_size == grid_vertex_buffer_alloc) 
  { 
    grid_vertex_buffer = (GridVertex*)realloc(grid_vertex_buffer,sizeof(GridVertex)*grid_vertex_buffer_alloc*2); 
    if (!grid_vertex_buffer) 
    { 
      fprintf(stderr,"FATAL ERROR: realloc grid_vertex_buffer with %d failed.\n",grid_vertex_buffer_alloc*2); 
      exit(0); 
    } 
    grid_vertex_buffer_alloc *= 2; 
  } 
  int index = grid_vertex_buffer_size; 
  grid_vertex_buffer[index].v[0] = 0.0f; 
  grid_vertex_buffer[index].v[1] = 0.0f; 
  grid_vertex_buffer[index].v[2] = 0.0f; 
  grid_vertex_buffer[index].number = 0; 
  grid_vertex_buffer_size++; 
  return index; 
} 
 
static void destroyGridVertexBuffer() 
{ 
  grid_vertex_buffer_size = 0; 
  grid_vertex_buffer_alloc = 0; 
  if (grid_vertex_buffer) 
  { 
    free(grid_vertex_buffer); 
  } 
  grid_vertex_buffer = 0; 
} 
 
static int* triangle_buffer = 0; 
static int triangle_buffer_size = 0; 
static int triangle_buffer_alloc = 0; 
 
static void initTriangleBuffer(int alloc) 
{ 
  if (triangle_buffer) 
  { 
    if (triangle_buffer_alloc < alloc) 
    { 
      triangle_buffer_alloc = alloc; 
      free(triangle_buffer); 
      triangle_buffer = (int*)malloc(sizeof(int)*triangle_buffer_alloc*3); 
    } 
  } 
  else 
  { 
    triangle_buffer_alloc = alloc; 
    triangle_buffer = (int*)malloc(sizeof(int)*triangle_buffer_alloc*3); 
  } 
  triangle_buffer_size = 0; 
} 
 
static int allocTriangle() 
{ 
  if (triangle_buffer_size == triangle_buffer_alloc) 
  { 
    triangle_buffer = (int*)realloc(triangle_buffer,sizeof(int)*triangle_buffer_alloc*3*2); 
    if (!triangle_buffer) 
    { 
      fprintf(stderr,"FATAL ERROR: realloc triangle_buffer with %d failed.\n",triangle_buffer_alloc*2); 
      exit(0); 
    } 
    triangle_buffer_alloc *= 2; 
  } 
  int index = triangle_buffer_size; 
  triangle_buffer_size++; 
  return index; 
} 
 
static void destroyTriangleBuffer() 
{ 
  triangle_buffer_size = 0; 
  triangle_buffer_alloc = 0; 
  if (triangle_buffer) 
  { 
    free(triangle_buffer); 
  } 
  triangle_buffer = 0; 
} 
 
static int render_vertex_buffer_alloc = 1024; 
static RenderVertex* render_vertex_buffer_next = 0; 
 
static RenderVertex* allocRenderVertex(float* v) 
{ 
  if (render_vertex_buffer_next == 0) 
  { 
    render_vertex_buffer_next = (RenderVertex*)malloc(sizeof(RenderVertex)*render_vertex_buffer_alloc); 
    if (render_vertex_buffer_next == 0) 
    { 
      fprintf(stderr,"malloc for render vertex buffer failed\n"); 
      return 0; 
    } 
    for (int i = 0; i < render_vertex_buffer_alloc; i++) 
    { 
      render_vertex_buffer_next[i].buffer_next = &(render_vertex_buffer_next[i+1]); 
    } 
    render_vertex_buffer_next[render_vertex_buffer_alloc-1].buffer_next = 0; 
    render_vertex_buffer_alloc = 2*render_vertex_buffer_alloc; 
  } 
  // get pointer to next available vertex 
  RenderVertex* vertex = render_vertex_buffer_next; 
  render_vertex_buffer_next = vertex->buffer_next; 
 
  VecCopy3fv(vertex->v, v); 
   
  return vertex; 
} 
 
static void deallocRenderVertex(RenderVertex* vertex) 
{ 
  vertex->buffer_next = render_vertex_buffer_next; 
  render_vertex_buffer_next = vertex; 
} 
 
void SavePPM(char *FileName, unsigned char* Colour, int Width, int Height) 
{ 
  FILE *fp = fopen(FileName, "wb"); 
  fprintf(fp, "P6\n%d %d\n255\n", Width, Height); 
  int NumRowPixels = Width*3; 
  for (int i=(Height-1)*Width*3; i>=0; i-=(Width*3)) 
  { 
    fwrite(&(Colour[i]),1,NumRowPixels,fp); 
  } 
  fclose(fp); 
} 
 
void InitColors() 
{ 
  colours_diffuse[0][0] = 0.0f; colours_diffuse[0][1] = 0.0f; colours_diffuse[0][2] = 0.0f; colours_diffuse[0][3] = 1.0f; // black 
  colours_diffuse[1][0] = 0.6f; colours_diffuse[1][1] = 0.0f; colours_diffuse[1][2] = 0.0f; colours_diffuse[1][3] = 1.0f; // red 
  colours_diffuse[2][0] = 0.0f; colours_diffuse[2][1] = 0.8f; colours_diffuse[2][2] = 0.0f; colours_diffuse[2][3] = 1.0f; // green 
  colours_diffuse[3][0] = 0.0f; colours_diffuse[3][1] = 0.0f; colours_diffuse[3][2] = 0.6f; colours_diffuse[3][3] = 1.0f; // blue 
  colours_diffuse[4][0] = 0.6f; colours_diffuse[4][1] = 0.6f; colours_diffuse[4][2] = 0.0f; colours_diffuse[4][3] = 1.0f; // yellow 
  colours_diffuse[5][0] = 0.6f; colours_diffuse[5][1] = 0.0f; colours_diffuse[5][2] = 0.6f; colours_diffuse[5][3] = 1.0f; // purple 
  colours_diffuse[6][0] = 0.0f; colours_diffuse[6][1] = 0.6f; colours_diffuse[6][2] = 0.6f; colours_diffuse[6][3] = 0.3f; // cyan (tranparent) 
  colours_diffuse[7][0] = 0.7f; colours_diffuse[7][1] = 0.7f; colours_diffuse[7][2] = 0.7f; colours_diffuse[7][3] = 1.0f; // white 
  colours_diffuse[8][0] = 0.2f; colours_diffuse[8][1] = 0.2f; colours_diffuse[8][2] = 0.6f; colours_diffuse[8][3] = 1.0f; // light blue 
  colours_diffuse[9][0] = 0.9f; colours_diffuse[9][1] = 0.4f; colours_diffuse[9][2] = 0.7f; colours_diffuse[9][3] = 0.5f; // violett 
   
  colours_white[0] = 0.7f; colours_white[1] = 0.7f; colours_white[2] = 0.7f; colours_white[3] = 1.0f; // white 
  colours_light_blue[0] = 0.2f; colours_light_blue[1] = 0.2f; colours_light_blue[2] = 0.6f; colours_light_blue[3] = 1.0f; // light blue 
}  
 
void InitLight() 
{ 
  float intensity[] = {1,1,1,1}; 
  float position[] = {1,1,5,0}; // directional behind the viewer 
  glLightfv(GL_LIGHT0,GL_DIFFUSE,intensity); 
  glLightfv(GL_LIGHT0,GL_SPECULAR,intensity); 
  glLightfv(GL_LIGHT0,GL_POSITION,position); 
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE); 
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE); 
} 
 
#ifdef _WIN32 
extern "C" FILE* fopenGzipped(const char* filename, const char* mode); 
#endif 
 
void usage() 
{ 
  fprintf(stderr,"usage:\n"); 
  fprintf(stderr,"sp_viewer -i cloud.spb\n"); 
  fprintf(stderr,"sp_viewer -i cloud.spa.gz -del3\n"); 
  fprintf(stderr,"sp_viewer -i terrain.spb -every 100\n"); 
  fprintf(stderr,"sp_viewer -i terrain.spb.gz -del2\n"); 
  fprintf(stderr,"sp_viewer -i terrain.obj.gz -del2 -flatten\n"); 
  fprintf(stderr,"sp_viewer -i cloud.spb -del3 -steps 500\n"); 
  fprintf(stderr,"sp_viewer -win 1600 1200 -i cloud.obj\n"); 
	fprintf(stderr,"sp_viewer -c -hi A.1.sma -i A.1.spa  (for clarkson)\n");
	fprintf(stderr,"sp_viewer -csm -hi A.1.sma -i A.1.spa  (for clarkson but only displays CT)\n");
  fprintf(stderr,"sp_viewer -h\n"); 
  fprintf(stderr,"\n"); 
  fprintf(stderr,"in addition to our streaming point formats SPA and SPB this tool\n"); 
  fprintf(stderr,"supports SMx, OBJ, SMF, PLY  meshes (optionally gzipped).\n"); 
} 
 
FILE* file = 0; 

//these are settings are used for clarkson
FILE* file_header = 0;

float* clarkson2d_mesh_verts = 0; // CDT vertices
int clarkson2d_mesh_nverts = 0; // # of CDT vertices
int* clarkson2d_mesh_faces = 0; // CDT finalized regions
int clarkson2d_mesh_nfaces = 0; // # of CDT finalized regions

int num_regions_finalized = 0;

bool g_clarkson = false;
bool g_sm_only = false;               // only sm file is provided.

//=========================
// The following functions are for clarkson streaming


/*****************************************************************************/
/*                                                                           */
/*  findcircumcircle()   Find the circumcircle of a triangle.                */
/*                                                                           */
/*  Return the cetner and the radius of the circumcircle of the triangle.    */
/*  This function is copied and modified from triangle.                      */ 
/*                                                                           */
/*****************************************************************************/

void findcircumcircle(float* torg, float* tdest, float* tapex, float* circumcenter, float* radius)
{

  float xdo, ydo, xao, yao;
  float dodist, aodist, dadist;
  float denominator;
  float dx, dy, dxoff, dyoff;

  /* Compute the circumcenter of the triangle. */
  xdo = tdest[0] - torg[0];
  ydo = tdest[1] - torg[1];
  xao = tapex[0] - torg[0];
  yao = tapex[1] - torg[1];
  dodist = xdo * xdo + ydo * ydo;
  aodist = xao * xao + yao * yao;
  dadist = (tdest[0] - tapex[0]) * (tdest[0] - tapex[0]) + (tdest[1] - tapex[1]) * (tdest[1] - tapex[1]);
  denominator = 0.5 / (xdo * yao - xao * ydo); // TODO: watch out for the possible dividing by 0
  dx = (yao * dodist - ydo * aodist) * denominator;
  dy = (xdo * aodist - xao * dodist) * denominator;
  circumcenter[0] = torg[0] + dx;
  circumcenter[1] = torg[1] + dy;
  *radius = sqrt(dx*dx+dy*dy);
}

#define PI 3.14159
void drawCircle(float centerX, float centerY, float radius, float* rgb, int solid) {
  /* draw a circle */
  if(solid)
    glBegin(GL_TRIANGLE_FAN);
  else
    glBegin(GL_LINE_STRIP);

  // Set the drawing color
  glColor3f(rgb[0], rgb[1], rgb[2]);
  for( double angle=0.0; angle<2.1 ; angle+=0.01 ) {
    glVertex2d(centerX +radius*sin(angle*PI), centerY+radius*cos(angle*PI));
  }
  glEnd();
}

// a and b > 0
float floatR(float a, float b) {
  while(a>=b) a-=b;
  return a;
}
float floatQ(float a, float b) {
  return (a - floatR(a,b)) / b;
}

void setRegionColor(float *color) {
  float colorratio = ((float) num_regions_finalized ) / clarkson2d_mesh_nfaces;

  int colorIntScale = (int) (1000*colorratio);
  /* colorration is now 0 ~ 1000 */
  color[0] = ((colorIntScale % 10) / 9.0);
  color[1] = ((colorIntScale/10 % 10) / 9.0);
  color[2] = ((colorIntScale/100 % 10) / 9.0);
  //  fprintf(stderr, "%f\n", colorratio );  
  //  fprintf(stderr, "%f , %f , %f         %d\n", color[0] , color[1], color[2], colorIntScale );
}

// pass initframe=1 for displaying frame
void drawRegion(int index, int initframe) {

  //--------- draw CDT finalization regions -----------------
  glColor3f(0.3f,0.3f,0.3f);

  float r;
  float circumcenter[2];
  float circlecolor[3];
  int ia, ib, ic;
  // note clarkson2d_mesh_faces considers first non-pinf pt at 0
  // but clarkson2d_mesh_verts considers it at 1 where pinf is at 0
  ia = clarkson2d_mesh_faces[index*3  ]+1; // offset due to pinf
  ib = clarkson2d_mesh_faces[index*3+1]+1; // offset due to pinf
  ic = clarkson2d_mesh_faces[index*3+2]+1; // offset due to pinf

  if( initframe ) {// for initial frame
    circlecolor[0]=0.9;
    circlecolor[1]=0.9;
    circlecolor[2]=0.9;
  } else {
    setRegionColor(circlecolor);
  }
  if( ia && ib && ic ) {
    // draw the triangle only if it's not infinite triangle
    glColor3f(0.8f, 0.8f, 0.8f);
    glBegin(GL_LINE_LOOP);

    //    printf("%d\n",clarkson2d_mesh_nverts);//clarkson2d_mesh_v_count);
    /*
    if( ia<0 || clarkson2d_mesh_verts[ia*3]==NULL ) {
      fprintf(stderr, "fatal error! invalid ia %d\n", ia);
      for(int i=ia-1;i<=ia+1;i++)
	fprintf(stderr, "%d->%d  ", i,  clarkson2d_mesh_verts[i*3]);
      exit(0);
      }
    if( ib<0 || clarkson2d_mesh_verts[ib*3]==NULL ) { fprintf(stderr, "fatal error! invalid ib"); exit(0); }
    if( ic<0 || clarkson2d_mesh_verts[ic*3]==NULL ) { fprintf(stderr, "fatal error! invalid ic"); exit(0); }
    */
    glVertex3fv( &(clarkson2d_mesh_verts[ia*3]) );
    glVertex3fv( &(clarkson2d_mesh_verts[ib*3]) );
    glVertex3fv( &(clarkson2d_mesh_verts[ic*3]) );
    glEnd();
    findcircumcircle( &(clarkson2d_mesh_verts[ia*3]), &(clarkson2d_mesh_verts[ib*3]),
		      &(clarkson2d_mesh_verts[ic*3]), circumcenter, &r);
		if (!g_sm_only) {
			drawCircle(circumcenter[0],circumcenter[1],r,circlecolor,0/*!initframe*/);
		}
    //      printf("%d=>(c=%.3f %.3f %.3f, r=%f)\n",i, circumcenter[0],circumcenter[1],circumcenter[2],r );

  } else { // it's in infinite triangle. dont draw it.
    /*
    glBegin(GL_LINE_LOOP);
    glColor3f(0.0f, 1.0f, 0.0f);
    if( ia )  glVertex3fv( &(clarkson2d_mesh_verts[ia*3]) );
    if( ib )  glVertex3fv( &(clarkson2d_mesh_verts[ib*3]) );
    if( ic )  glVertex3fv( &(clarkson2d_mesh_verts[ic*3]) );
    glEnd();
    */
  }

}

//=======================================
// Read clarkson. g_clarkson must be true
void read() {
  bool ispa = false;
  bool ispb = false;

  // open input point file 
	// g_clarkson must be true

  if (file_name || ispa || ispb)
  {
    if (file_name) // if file name is specified
    {
      if (strstr(file_name, ".spa") || ispa) // ascii file
      {
        file = fopen(file_name, "r");
      }
      else if (strstr(file_name, ".spb") || ispb) // binary file
      {
        file = fopen(file_name, "rb");
      }
      else // error occurs!
      {
        if (file_name)
        {
          fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name);
        }
        else
        {
          fprintf(stderr,"ERROR: no ouput format specified\n");
        }
        exit(0);
      }

      if (file == 0) // if files cannot open for write
      {
        fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name);
				exit(0);
      }
    }
    else // file name is from stdin
    {
      file = stdin;
    }

		if(file_name) // if file is not from stdin
		{
			if (strstr(file_name, ".spa") || ispa)
			{
				SPreader_spa* spreader_spa = new SPreader_spa();
				spreader_spa->open(file);
				spreader = spreader_spa;
			}
			else if (strstr(file_name, ".spb") || ispb)
			{
				SPreader_spb* spreader_spb = new SPreader_spb();
				spreader_spb->open(file);
				spreader = spreader_spb;
			}
		}
  }

  // make sure we opened a point cloud with the right type of finalization

  // open output header file ... important ... do not open any earlier to make sure header file is written

  if (file_name_header)
  {
    if (strstr(file_name_header, ".sma"))
    {
      file_header = fopen(file_name_header, "r");
    }
    else if (strstr(file_name_header, ".smb"))
    {
     file_header = fopen(file_name_header, "rb");
    }
    else
    {
      fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name_header);
      exit(0);
    }
    if (file_header == 0)
    {
      fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_header);
      exit(0);
    }
  }
  else
  {
    fprintf(stderr,"ERROR: no file name given for storing header of finalizemethod \n");
    exit(0);
  }

  if (strstr(file_name_header, ".sma"))
  {
    SMreader_sma* smreader_sma = new SMreader_sma();
    smreader_sma->open(file_header);
    smreader = smreader_sma;
  }
  else if (strstr(file_name_header, ".smb"))
  {
    SMreader_smb* smreader_smb = new SMreader_smb();
    smreader_smb->open(file_header);
    smreader = smreader_smb;
  }

  // simple example of reading in a streaming point cloud with SP_CLARKSON_2D finalization

  // *first* we read the header
  
  clarkson2d_mesh_verts = 0;
  clarkson2d_mesh_nverts = smreader->nverts;

  if (clarkson2d_mesh_nverts > 0)
  {
    clarkson2d_mesh_verts = (float*)malloc(clarkson2d_mesh_nverts*sizeof(float)*3);
  }
  else 
  {
    fprintf(stderr,"ERROR: the clarkson finalization mesh does not have any vertices\n");
    exit(0);
  }

  clarkson2d_mesh_faces = 0;
  clarkson2d_mesh_nfaces = smreader->nfaces;

  if (clarkson2d_mesh_nfaces > 0)
  {
    clarkson2d_mesh_faces = (int*)malloc(clarkson2d_mesh_nfaces*sizeof(int)*3);
  }
  else 
  {
    fprintf(stderr,"ERROR: the clarkson finalization mesh does not have any faces\n");
    exit(0);
  }

  SMevent event;

  int clarkson2d_mesh_v_count = 0;
  int clarkson2d_mesh_f_count = 0;

  while ((event = smreader->read_element()))
  {
    if (event == SM_TRIANGLE)
    {
      VecCopy3iv(&(clarkson2d_mesh_faces[clarkson2d_mesh_f_count*3]), smreader->t_idx);
      clarkson2d_mesh_f_count++;
    }
    else if (event == SM_VERTEX)
    {
      VecCopy3fv(&(clarkson2d_mesh_verts[clarkson2d_mesh_v_count*3]), smreader->v_pos_f);
      clarkson2d_mesh_v_count++;
    }
  }

  smreader->close();
  fclose(file_header);
  delete smreader;

  if (clarkson2d_mesh_v_count != clarkson2d_mesh_nverts)
  {
    fprintf(stderr,"WARNING: clarkson2d_mesh_v_count %d is different from clarkson2d_mesh_nverts %d\n", clarkson2d_mesh_v_count, clarkson2d_mesh_nverts);
    exit(0);
  }

  if (clarkson2d_mesh_f_count != clarkson2d_mesh_nfaces)
  {
    fprintf(stderr,"WARNING: clarkson2d_mesh_f_count %d is different from clarkson2d_mesh_nfaces %d\n", clarkson2d_mesh_f_count, clarkson2d_mesh_nfaces);
    exit(0);
  }

  fprintf(stderr,"INFO: read %d vertices and %d triangles from clarkson finalization header\n", clarkson2d_mesh_v_count, clarkson2d_mesh_f_count);


  num_regions_finalized = 0;
  //  last_point_in_region = (int*) malloc(clarkson2d_mesh_f_count*sizeof(int));

  // *then* we read the points... do this in vizContinue... 
  //allocate memory first
  //  int n = 1024;
  //  initGridVertexBuffer(n);
  //  printf("%d vertices are allocated initially\n", n);

  //------------------ find bounding box --------------------

  //----------- we can display the CDT already -------------------

  //============================ display =============================
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

  glDisable(GL_LIGHTING); 
  glDisable(GL_LIGHT0); 
  glDisable(GL_NORMALIZE); 
  glDisable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); // clean the matrix
  glOrtho(-1, 1, -1, 1, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); // clean the matrix
  /*
  for( int i=0; i<clarkson2d_mesh_nfaces; i++) {
    drawRegion(i, 1); // 1 for displaying inital frame
  }
  */

  // get # of pionts from spreader
  npoints = spreader->npoints; 
 
  fprintf(stderr,"npoints %d\n",npoints); 

  if (EVERY_NTH_STEP == -1) 
  { 
    EVERY_NTH_STEP = npoints / EXACTLY_N_STEPS; 
  } 
  if (EVERY_NTH_STEP == 0) 
  { 
    EVERY_NTH_STEP = 1; 
  } 
  NEXT_STEP = EVERY_NTH_STEP; 

  initGridVertexBuffer(1024); 

  grid_hash->clear(); 
  

	// now calculate the bounding box
  SPevent pevent;

  spreader->bb_min_f = new float[3];
  spreader->bb_max_f = new float[3];
  
	// get the first point, and assign min and max to it
  while ((pevent = spreader->read_event()) ) {
    if (pevent == SP_POINT) {
      VecCopy3fv(spreader->bb_min_f, spreader->p_pos_f);
      VecCopy3fv(spreader->bb_max_f, spreader->p_pos_f);
      break;
    }
  }

  fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", spreader->bb_min_f[0], spreader->bb_min_f[1], spreader->bb_min_f[2]); 
  fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", spreader->bb_max_f[0], spreader->bb_max_f[1], spreader->bb_max_f[2]);

  while ((pevent = spreader->read_event()) ) {
    if (pevent == SP_POINT) {
      VecUpdateMinMax3fv(spreader->bb_min_f, spreader->bb_max_f, spreader->p_pos_f);
    }
  }

  fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", spreader->bb_min_f[0], spreader->bb_min_f[1], spreader->bb_min_f[2]); 
  fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", spreader->bb_max_f[0], spreader->bb_max_f[1], spreader->bb_max_f[2]);


  VecCopy3fv(boundingBoxMin, spreader->bb_min_f); 
  VecCopy3fv(boundingBoxMax, spreader->bb_max_f); 
 
  if ((boundingBoxMax[1]-boundingBoxMin[1]) > (boundingBoxMax[0]-boundingBoxMin[0])) 
  { 
    if ((boundingBoxMax[1]-boundingBoxMin[1]) > (boundingBoxMax[2]-boundingBoxMin[2])) 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[1]-boundingBoxMin[1]); 
    } 
    else 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[2]-boundingBoxMin[2]); 
    } 
  } 
  else 
  { 
    if ((boundingBoxMax[0]-boundingBoxMin[0]) > (boundingBoxMax[2]-boundingBoxMin[2])) 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[0]-boundingBoxMin[0]); 
    } 
    else 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[2]-boundingBoxMin[2]); 
    } 
  } 
  boundingBoxTranslateX = - boundingBoxScale * (boundingBoxMin[0] + 0.5f * (boundingBoxMax[0]-boundingBoxMin[0])); 
  boundingBoxTranslateY = - boundingBoxScale * (boundingBoxMin[1] + 0.5f * (boundingBoxMax[1]-boundingBoxMin[1])); 
  boundingBoxTranslateZ = - boundingBoxScale * (boundingBoxMin[2] + 0.5f * (boundingBoxMax[2]-boundingBoxMin[2])); 

  //--- close file and reopen for reading from beginning later ---
  spreader->close();
  fclose(file);

  if (strstr(file_name, ".spa") || ispa)
    {
      file = fopen(file_name, "r");
    }
  else if (strstr(file_name, ".spb") || ispb)
    {
      file = fopen(file_name, "rb");
    }
  if (file == 0)
    {
      fprintf(stderr,"ERROR: cannot open '%s' after closing it once\n", file_name);
      exit(0);
    }
  spreader->open(file);
  //------------------------------------------------------------

}

//========================= Beginning of vizBegin ========================================

void vizBegin() 
{ 
  REPLAY_IT = 0; // just making sure 
  DIRTY_MESH = 1; 
 
  if (file_name == 0 && !isma && !ismb && !ismc) 
  { 
    fprintf(stderr,"ERROR: no input\n"); 
    exit(0); 
  } 
 
  if (file_name) 
  { 
    fprintf(stderr,"loading mesh '%s'...\n",file_name); 
    if (strstr(file_name, ".spa") || ispa) 
    { 
      if (strstr(file_name, ".gz")) 
      { 
        #ifdef _WIN32 
        file = fopenGzipped(file_name, "r"); 
        #else 
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
        exit(1); 
        #endif 
      } 
      else 
      { 
        file = fopen(file_name, "r"); 
      } 
      if (file == 0) 
      { 
        fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
        exit(1); 
      } 
      SPreader_spa* spreader_spa = new SPreader_spa(); 
      spreader_spa->open(file);

	printf("i have %d points\n", spreader_spa->npoints); 
 
      if (spreader_spa->npoints == -1 || spreader_spa->bb_min_f == 0 || spreader_spa->bb_max_f == 0) 
      { 
        if (spreader_spa->bb_min_f && spreader_spa->bb_max_f) 
        { 
          fprintf(stderr,"need additional pass to count the points\n"); 
          while (spreader_spa->read_event()); 
        } 
        else 
        { 
          SPevent event; 
 
          if (spreader_spa->npoints == -1) 
          { 
            fprintf(stderr,"need additional pass to count the points and compute the bounding box\n"); 
          } 
          else 
          { 
            fprintf(stderr,"need additional pass to compute the bounding box\n"); 
          } 
 
          spreader_spa->bb_min_f = new float[3]; 
          spreader_spa->bb_max_f = new float[3]; 
 
          while ((event = spreader_spa->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecCopy3fv(spreader_spa->bb_min_f, spreader_spa->p_pos_f); 
              VecCopy3fv(spreader_spa->bb_max_f, spreader_spa->p_pos_f); 
              break; 
            } 
          } 
          while ((event = spreader_spa->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecUpdateMinMax3fv(spreader_spa->bb_min_f, spreader_spa->bb_max_f, spreader_spa->p_pos_f); 
            } 
          } 
          fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", spreader_spa->bb_min_f[0], spreader_spa->bb_min_f[1], spreader_spa->bb_min_f[2]); 
          fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", spreader_spa->bb_max_f[0], spreader_spa->bb_max_f[1], spreader_spa->bb_max_f[2]); 
        } 
       
        spreader_spa->close(); 
        fclose(file); 
 
        if (strstr(file_name, ".gz")) 
        { 
          #ifdef _WIN32 
          file = fopenGzipped(file_name, "r"); 
          #else 
          fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
          exit(1); 
          #endif 
        } 
        else 
        { 
          file = fopen(file_name, "r"); 
        } 
        if (file == 0) 
        { 
          fprintf(stderr,"ERROR: cannot open %s a second time\n",file_name); 
          exit(1); 
        } 
        spreader_spa->open(file); 
      } 
      spreader = spreader_spa; 
    } 
    else if (strstr(file_name, ".spb") || ispb) 
    { 
      if (strstr(file_name, ".gz")) 
      { 
        #ifdef _WIN32 
        file = fopenGzipped(file_name, "rb"); 
        #else 
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
        exit(1); 
        #endif 
      } 
      else 
      { 
        file = fopen(file_name, "rb"); 
      } 
      if (file == 0) 
      { 
        fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
        exit(1); 
      } 
      SPreader_spb* spreader_spb = new SPreader_spb(); 
      spreader_spb->open(file); 
 
      if (spreader_spb->npoints == -1 || spreader_spb->bb_min_f == 0 || spreader_spb->bb_max_f == 0) 
      { 
        if (spreader_spb->bb_min_f && spreader_spb->bb_max_f)
        { 
          fprintf(stderr,"need additional pass to count the points\n"); 
          while (spreader_spb->read_event()); 
        } 
        else 
        { 
          SPevent event; 
 
          if (spreader_spb->npoints == -1) 
          { 
            fprintf(stderr,"need additional pass to count the points and compute the bounding box\n"); 
          } 
          else 
          { 
            fprintf(stderr,"need additional pass to compute the bounding box\n"); 
          } 
 
          spreader_spb->bb_min_f = new float[3]; 
          spreader_spb->bb_max_f = new float[3]; 
 
          while ((event = spreader_spb->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecCopy3fv(spreader_spb->bb_min_f, spreader_spb->p_pos_f); 
              VecCopy3fv(spreader_spb->bb_max_f, spreader_spb->p_pos_f); 
              break; 
            } 
          } 
          while ((event = spreader_spb->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecUpdateMinMax3fv(spreader_spb->bb_min_f, spreader_spb->bb_max_f, spreader_spb->p_pos_f); 
            } 
          } 
          fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", spreader_spb->bb_min_f[0], spreader_spb->bb_min_f[1], spreader_spb->bb_min_f[2]); 
          fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", spreader_spb->bb_max_f[0], spreader_spb->bb_max_f[1], spreader_spb->bb_max_f[2]); 
        } 
       
        spreader_spb->close(); 
        fclose(file); 
 
        if (strstr(file_name, ".gz")) 
        { 
          #ifdef _WIN32 
          file = fopenGzipped(file_name, "rb"); 
          #else 
          fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
          exit(1); 
          #endif 
        } 
        else 
        { 
          file = fopen(file_name, "rb"); 
        } 
        if (file == 0) 
        { 
          fprintf(stderr,"ERROR: cannot open %s a second time\n",file_name); 
          exit(1); 
        } 
        spreader_spb->open(file); 
      } 
      spreader = spreader_spb; 
    } 
    else if (strstr(file_name, ".node") || inode) 
    { 
      if (strstr(file_name, ".gz")) 
      { 
        #ifdef _WIN32 
        file = fopenGzipped(file_name, "r"); 
        #else 
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
        exit(1); 
        #endif 
      } 
      else 
      { 
        file = fopen(file_name, "r"); 
      } 
      if (file == 0) 
      { 
        fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
        exit(1); 
      } 
      SPreader_node* spreader_node = new SPreader_node(); 
      spreader_node->open(file); 
 
      if (spreader_node->npoints == -1 || spreader_node->bb_min_f == 0 || spreader_node->bb_max_f == 0) 
      { 
        if (spreader_node->bb_min_f && spreader_node->bb_max_f) 
        { 
          fprintf(stderr,"need additional pass to count the points\n"); 
          while (spreader_node->read_event()); 
        } 
        else 
        { 
          SPevent event; 
 
          if (spreader_node->npoints == -1) 
          { 
            fprintf(stderr,"need additional pass to count the points and compute the bounding box\n"); 
          } 
          else 
          { 
            fprintf(stderr,"need additional pass to compute the bounding box\n"); 
          } 
 
          spreader_node->bb_min_f = new float[3]; 
          spreader_node->bb_max_f = new float[3]; 
 
          while ((event = spreader_node->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecCopy3fv(spreader_node->bb_min_f, spreader_node->p_pos_f); 
              VecCopy3fv(spreader_node->bb_max_f, spreader_node->p_pos_f); 
              break; 
            } 
          } 
          while ((event = spreader_node->read_event()) )
          { 
            if (event == SP_POINT) 
            { 
              VecUpdateMinMax3fv(spreader_node->bb_min_f, spreader_node->bb_max_f, spreader_node->p_pos_f); 
            } 
          } 
          fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", spreader_node->bb_min_f[0], spreader_node->bb_min_f[1], spreader_node->bb_min_f[2]); 
          fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", spreader_node->bb_max_f[0], spreader_node->bb_max_f[1], spreader_node->bb_max_f[2]); 
        } 
       
        spreader_node->close(); 
        fclose(file); 
 
        if (strstr(file_name, ".gz")) 
        { 
          #ifdef _WIN32 
          file = fopenGzipped(file_name, "r"); 
          #else 
          fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
          exit(1); 
          #endif 
        } 
        else 
        { 
          file = fopen(file_name, "r"); 
        } 
        if (file == 0) 
        { 
          fprintf(stderr,"ERROR: cannot open %s a second time\n",file_name); 
          exit(1); 
        } 
        spreader_node->open(file); 
      } 
      spreader = spreader_node; 
    } 
    else if (strstr(file_name, ".raw_d")) 
    { 
      if (strstr(file_name, ".gz")) 
      { 
        #ifdef _WIN32 
        file = fopenGzipped(file_name, "rb"); 
        #else 
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
        exit(1); 
        #endif 
      } 
      else 
      { 
        file = fopen(file_name, "rb"); 
      } 
      if (file == 0) 
      { 
        fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
        exit(1); 
      } 
      SPreader_raw_d* spreader_raw_d = new SPreader_raw_d(); 
      spreader_raw_d->open(file); 
      if (file_name_header) 
      { 
        FILE* file_header = fopen(file_name_header, "rb"); 
        spreader_raw_d->load_header(file_header); 
        fclose(file_header); 
      } 
      else 
      { 
        fprintf(stderr,"no header_file_name ...\n"); 
        exit(1); 
      } 
      spreader = spreader_raw_d; 
    } 

    //============================ SMA format Reading starts here ===============================
    else if (strstr(file_name, ".sma") || strstr(file_name, ".obj") || strstr(file_name, ".smf") || isma) 
    { 
      // ---------- open the file ---------------
      if (strstr(file_name, ".gz")) 
      { 
        #ifdef _WIN32 
        file = fopenGzipped(file_name, "r"); 
        #else 
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
        exit(1); 
        #endif 
      } 
      else 
      { 
        file = fopen(file_name, "r"); 
      } 
      if (file == 0) 
      { 
        fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
        exit(1); 
      } 

      //--------- construct smreader_sma with information in file  -------------------
      SMreader_sma* smreader_sma = new SMreader_sma(); 
      smreader_sma->open(file); 
 
      if (smreader_sma->nfaces == -1 || smreader_sma->nverts == -1 || smreader_sma->bb_min_f == 0 || smreader_sma->bb_max_f == 0) 
      { 
        if (smreader_sma->bb_min_f && smreader_sma->bb_max_f) 
        { 
          fprintf(stderr,"need additional pass to count verts and faces\n"); 
          while (smreader_sma->read_element()); 
        } 
        else 
        { 
          SMevent event; 
 
          if (smreader_sma->nfaces == -1 || smreader_sma->nverts == -1) 
          { 
            fprintf(stderr,"need additional pass to count verts and faces and compute bounding box\n"); 
          } 
          else 
          { 
            fprintf(stderr,"need additional pass to compute bounding box\n"); 
          } 
 
          smreader_sma->bb_min_f = new float[3]; 
          smreader_sma->bb_max_f = new float[3]; 

	  // seems like copying the information of bounding box to v_pos_f
	  while ((event = smreader_sma->read_element()) )
          { 
            if (event == SM_VERTEX) 
            { 
              VecCopy3fv(smreader_sma->bb_min_f, smreader_sma->v_pos_f); 
              VecCopy3fv(smreader_sma->bb_max_f, smreader_sma->v_pos_f); 
              break; 
            } 
          } 

          while ((event = smreader_sma->read_element()) )
          { 
            if (event == SM_VERTEX) 
            { 
              VecUpdateMinMax3fv(smreader_sma->bb_min_f, smreader_sma->bb_max_f, smreader_sma->v_pos_f); 
            } 
          } 

          fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", smreader_sma->bb_min_f[0], smreader_sma->bb_min_f[1], smreader_sma->bb_min_f[2]); 
          fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", smreader_sma->bb_max_f[0], smreader_sma->bb_max_f[1], smreader_sma->bb_max_f[2]); 
        } 
       
	//-------- Smreader file is closed because we have all the information in Smreader now -------------
        smreader_sma->close(); 
        fclose(file); 

	//--------- open the file again ---------

        if (strstr(file_name, ".gz")) 
        { 
          #ifdef _WIN32 
          file = fopenGzipped(file_name, "r"); 
          #else 
          fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
          exit(1); 
          #endif 
        } 
        else 
        { 
          file = fopen(file_name, "r"); 
        }
        if (file == 0) 
        { 
          fprintf(stderr,"ERROR: cannot open %s a second time\n",file_name); 
          exit(1); 
        } 

        smreader_sma->open(file); 
      } 
      smreader = smreader_sma; 
    } 
    //======================== SMA format ends here =============================


    else if (strstr(file_name, ".smb") || ismb) 
    {
      if (strstr(file_name, ".gz"))
      {
        #ifdef _WIN32
        file = fopenGzipped(file_name, "rb");
        #else
        fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name);
        exit(1);
        #endif
      }
      else
      {
        file = fopen(file_name, "rb");
      }
      if (file == 0)
      {
        fprintf(stderr,"ERROR: cannot open %s\n",file_name);
        exit(1);
      }
      SMreader_smb* smreader_smb = new SMreader_smb(); 
      smreader_smb->open(file); 
 
      if (smreader_smb->nfaces == -1 || smreader_smb->nverts == -1 || smreader_smb->bb_min_f == 0 || smreader_smb->bb_max_f == 0) 
      { 
        if (smreader_smb->bb_min_f && smreader_smb->bb_max_f) 
        { 
          fprintf(stderr,"need additional pass to count verts and faces\n"); 
          while (smreader_smb->read_element()); 
        } 
        else 
        { 
          SMevent event; 
 
          if (smreader_smb->nfaces == -1 || smreader_smb->nverts == -1) 
          { 
            fprintf(stderr,"need additional pass to count verts and faces and compute bounding box\n"); 
          } 
          else 
          { 
            fprintf(stderr,"need additional pass to compute bounding box\n"); 
          } 
 
          smreader_smb->bb_min_f = new float[3]; 
          smreader_smb->bb_max_f = new float[3]; 
 
          while ((event = smreader_smb->read_element()) )
          { 
            if (event == SM_VERTEX) 
            { 
              VecCopy3fv(smreader_smb->bb_min_f, smreader_smb->v_pos_f); 
              VecCopy3fv(smreader_smb->bb_max_f, smreader_smb->v_pos_f); 
              break; 
            } 
          } 
          while ((event = smreader_smb->read_element()) )
          { 
            if (event == SM_VERTEX) 
            { 
              VecUpdateMinMax3fv(smreader_smb->bb_min_f, smreader_smb->bb_max_f, smreader_smb->v_pos_f); 
            } 
          } 
          fprintf(stderr, "bb_min_f[0] = %gf; bb_min_f[1] = %gf; bb_min_f[2] = %gf;\n", smreader_smb->bb_min_f[0], smreader_smb->bb_min_f[1], smreader_smb->bb_min_f[2]); 
          fprintf(stderr, "bb_max_f[0] = %gf; bb_max_f[1] = %gf; bb_max_f[2] = %gf;\n", smreader_smb->bb_max_f[0], smreader_smb->bb_max_f[1], smreader_smb->bb_max_f[2]); 
        } 
 
        smreader_smb->close(); 
        fclose(file); 
 
        if (strstr(file_name, ".gz")) 
        { 
          #ifdef _WIN32 
          file = fopenGzipped(file_name, "rb"); 
          #else 
          fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
          exit(1); 
          #endif 
        } 
        else 
        { 
          file = fopen(file_name, "rb"); 
        } 
        if (file == 0) 
        { 
          fprintf(stderr,"ERROR: cannot open %s a second time\n",file_name); 
          exit(1); 
        } 
        smreader_smb->open(file); 
      } 
      smreader = smreader_smb; 
    } 
    else 
    { 
      fprintf(stderr,"ERROR: cannot guess input type from file name %s \n", file_name); 
      exit(1); 
    } 
  }
  else 
  { 
    fprintf(stderr,"loading mesh from stdin...\n"); 
    file = stdin; 
    if (isma) 
    { 
      SMreader_sma* smreader_sma = new SMreader_sma(); 
      smreader_sma->open(file); 
      smreader = smreader_sma; 
      isma = false; 
    } 
    else if (ismb) 
    { 
      SMreader_smb* smreader_smb = new SMreader_smb(); 
      smreader_smb->open(file); 
      smreader = smreader_smb; 
      ismb = false; 
    } 
  }

  //=============== finished reading the input files =================
 
  //=== if we are reading the mesh, reduce it to streaming points ====
  if (smreader) 
  { 
    SPconverter* spconverter = new SPconverter(); 
    spconverter->open(smreader); 
    spreader = spconverter; 
    smreader = 0; 
  } 
 
  if (spreader->bb_min_f == 0 || spreader->bb_max_f == 0) 
  { 
    fprintf(stderr,"ERROR: no bounding box info ... exiting\n"); 
    exit(0); 
  } 
 
  // maybe need to flatten 
 
  if (flatten) 
  { 
    spreader->bb_max_f[2] = spreader->bb_min_f[2] = 0; 
  } 
 
  // maybe read scattered 
 
  if (scatter) 
  { 
    SPreadScattered* spreadscattered = new SPreadScattered(); 
    spreadscattered->open(spreader, scatter); 
    spreader = spreadscattered; 
  } 
 
  // setup quantization grid 

  //TODO: we are using different setup cuz of diff data structure
 
  pq->SetPrecision(GRID_PRECISION); 
  pq->SetMinMax(spreader->bb_min_f,spreader->bb_max_f); 
  pq->SetupQuantizer(); 
 
  // setup streaming delaunay triangulator or finalization container 
 
  if (delaunay2) 
  { 
    spdelaunay2d = new SPdelaunay2D(); 
    spdelaunay2d->open(new SMwriter_nil()); 
    if (spreader->datatype == SP_DOUBLE) 
      spdelaunay2d->set_boundingbox(spreader->bb_min_d, spreader->bb_max_d); 
    else 
      spdelaunay2d->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f); 
    sscontainer2d = spdelaunay2d->ss2d; 
  } 
  else if (delaunay3) 
  { 
    spdelaunay3d = new SPdelaunay3D(); 
    spdelaunay3d->open(new SVwriter_nil()); 
    if (spreader->datatype == SP_DOUBLE) 
      spdelaunay3d->set_boundingbox(spreader->bb_min_d, spreader->bb_max_d); 
    else 
      spdelaunay3d->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f); 
    sscontainer3d = spdelaunay3d->ss3d; 
  } 
  else // set up the SScontainer2D
  { 
    sscontainer2d = new SScontainer2D(); 
    if (spreader->datatype == SP_DOUBLE) 
      sscontainer2d->open(spreader->bb_min_d, spreader->bb_max_d); 
    else 
      sscontainer2d->open(spreader->bb_min_f, spreader->bb_max_f); 
  } 
 
  // scale and translate bounding box for rendering 
 
  VecCopy3fv(boundingBoxMin, spreader->bb_min_f); 
  VecCopy3fv(boundingBoxMax, spreader->bb_max_f); 
 
  if ((boundingBoxMax[1]-boundingBoxMin[1]) > (boundingBoxMax[0]-boundingBoxMin[0])) 
  { 
    if ((boundingBoxMax[1]-boundingBoxMin[1]) > (boundingBoxMax[2]-boundingBoxMin[2])) 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[1]-boundingBoxMin[1]); 
    } 
    else 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[2]-boundingBoxMin[2]); 
    } 
  } 
  else 
  { 
    if ((boundingBoxMax[0]-boundingBoxMin[0]) > (boundingBoxMax[2]-boundingBoxMin[2])) 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[0]-boundingBoxMin[0]); 
    } 
    else 
    { 
      boundingBoxScale = 1.0f/(boundingBoxMax[2]-boundingBoxMin[2]); 
    } 
  } 
  boundingBoxTranslateX = - boundingBoxScale * (boundingBoxMin[0] + 0.5f * (boundingBoxMax[0]-boundingBoxMin[0])); 
  boundingBoxTranslateY = - boundingBoxScale * (boundingBoxMin[1] + 0.5f * (boundingBoxMax[1]-boundingBoxMin[1])); 
  boundingBoxTranslateZ = - boundingBoxScale * (boundingBoxMin[2] + 0.5f * (boundingBoxMax[2]-boundingBoxMin[2])); 
 
  // get # of pionts from spreader
  npoints = spreader->npoints; 
 
  fprintf(stderr,"npoints %d\n",npoints); 
 
  if (EVERY_NTH_STEP == -1) 
  { 
    EVERY_NTH_STEP = npoints / EXACTLY_N_STEPS; 
  } 
  if (EVERY_NTH_STEP == 0) 
  { 
    EVERY_NTH_STEP = 1; 
  } 
  NEXT_STEP = EVERY_NTH_STEP; 

  initGridVertexBuffer(1024); 

  grid_hash->clear(); 
} 

//=================================== END OF vizBegin =======================================================
 
void vizEnd() 
{ 
  REPLAY_IT = 0; // just making sure 
  REPLAY_COUNT = grid_vertex_buffer_size; 
  DIRTY_MESH = 0; 
//  fprintf(stderr,"grid: %dx%dx%d\n", pq->m_aiQuantRange[0], pq->m_aiQuantRange[1], pq->m_aiQuantRange[2]); 
//  fprintf(stderr,"  total # of grid cells: %d\n", grid_vertex_buffer_size); 
//  fprintf(stderr,"  average # vertex/cell: %.1f \n", 1.0f*spreader->p_count/grid_vertex_buffer_size); 

  glutSwapBuffers();


  if (delaunay2) 
  { 
    spdelaunay2d->close(); 
    delete spdelaunay2d; 
    spdelaunay2d = 0; 
    sscontainer2d = 0; 
  } 
  else if (delaunay3) 
  { 
    spdelaunay3d->close(); 
    delete spdelaunay3d; 
    spdelaunay3d = 0; 
    sscontainer3d = 0; 
  } 
  else if (sscontainer2d) 
  { 
    sscontainer2d->close(); 
    delete sscontainer2d; 
    sscontainer2d = 0; 
  } 
  else if (sscontainer3d) 
  { 
    sscontainer3d->close(); 
    delete sscontainer3d; 
    sscontainer3d = 0; 
  } 
 
  spreader->close(); 
  fclose(file); 
  delete spreader; 
  spreader = 0; 
} 
 
int vizContinue() 
{ 
  SPevent event = SP_ERROR; 
  my_grid_hash::iterator hash_element; 
 
  REPLAY_IT = 0; // just making sure 
 
  while ((event = spreader->read_event()) )
  {

    if (event == SP_POINT) // if a point is read
    {

			if( g_clarkson ) { // if using 2d clarkson
				//============================== only save the point to buffer ===============================

				//TODO: add this point to the vertex buffer
				//      printf("sp point\n");
				/*
					int vindex = allocGridVertex();
					VecCopy3fv( grid_vertex_buffer[vindex].v, spreader->p_pos_f);
				*/
				// reserve  grid_vertex_buffer[vindex].number for future use

				if( !g_sm_only ) {
					float color[3];
					setRegionColor( color );
					glEnable(GL_POINT_SMOOTH);
					glPointSize(POINT_SIZE);
					glBegin(GL_POINT);
					glColor3f( color[0], color[1], color[2]);      
					glVertex3fv( spreader->p_pos_f );
					glEnd();
					glutSwapBuffers();
				}
				//============================================================================================
			}
			else // other finalization methods
			{			
				int grid_oct; 
				int grid_pos[3]; 
				int grid_idx; 
 
				if (spreader->datatype == SP_DOUBLE) 
					{ 
						VecCopy3fv(spreader->p_pos_f, spreader->p_pos_d); 
					} 
 
				if (flatten) spreader->p_pos_f[2] = 0; 
 
				pq->EnQuantize(spreader->p_pos_f, grid_pos); 

				grid_oct = (grid_pos[0] << (GRID_PRECISION+GRID_PRECISION)) + (grid_pos[1] << GRID_PRECISION) + (grid_pos[2]); 

				// check if a grid cell for this point already exists 
				hash_element = grid_hash->find(grid_oct); 
				if (hash_element == grid_hash->end()) 
					{ 
						printf("allocating a grid because there is no grid for this point yet");
						grid_idx = allocGridVertex(); 
						// all following vertices falling into this grid cell find their grid vertex in the grid_hash 
						grid_hash->insert(my_grid_hash::value_type(grid_oct, grid_idx)); 
					} 
				else 
					{ 
						grid_idx = (*hash_element).second; 
					} 
 
				// all following triangles can find this vertex in the map_hash 

				VecSelfAdd3fv(grid_vertex_buffer[grid_idx].v, spreader->p_pos_f); 
				grid_vertex_buffer[grid_idx].number++; 


				if (rand()%17 == 5 && rand()%21 == 3) {
					VecCopy3fv(grid_vertex_buffer[grid_idx].v, spreader->p_pos_f);
					grid_vertex_buffer[grid_idx].number = 1;
				} 

				// ---------- check datatype, double? ----------------------------------
				if (spreader->datatype == SP_DOUBLE) 
					{ 
						if (delaunay2) 
							{ 
								spdelaunay2d->write_point(spreader->p_pos_d); 
							} 
						else if (delaunay3) 
							{ 
								spdelaunay3d->write_point(spreader->p_pos_d); 
							} 
						else if (sscontainer2d && sscontainer2d->is_finalized(spreader->p_pos_d)) 
							{ 
								int cell_idx = sscontainer2d->get_idx(spreader->p_pos_d, 4); 
								fprintf(stderr,"incoming point %d was %d in finalized space (%f %f %f)\n",  spreader->p_count, cell_idx, spreader->p_pos_d[0], spreader->p_pos_d[1], spreader->p_pos_d[2]); 
								SpecialPoint = spreader->p_pos_f; 
							} 
						else if (sscontainer3d && sscontainer3d->is_point_finalized(spreader->p_pos_d)) 
							{ 
								fprintf(stderr,"incoming point %d was in finalized space (%f %f %f)\n",  spreader->p_count, spreader->p_pos_d[0], spreader->p_pos_d[1], spreader->p_pos_d[2]); 
								SpecialPoint = spreader->p_pos_f; 
							} 
					} 
				else // not SP_DOUBLE
					{ 
						if (delaunay2) 
							{ 
								spdelaunay2d->write_point(spreader->p_pos_f); 
							} 
						else if (delaunay3) 
							{ 
								spdelaunay3d->write_point(spreader->p_pos_f); 
							} 
						else if (sscontainer2d && sscontainer2d->is_finalized(spreader->p_pos_f)) 
							{ 
								fprintf(stderr,"incoming point %d was in finalized space (%g %g %g)\n",  spreader->p_count, spreader->p_pos_f[0], spreader->p_pos_f[1], spreader->p_pos_f[2]); 
								SpecialPoint = spreader->p_pos_f; 
							} 
						else if (sscontainer3d && sscontainer3d->is_point_finalized(spreader->p_pos_f)) 
							{ 
								fprintf(stderr,"incoming point %d was in finalized space (%g %g %g)\n",  spreader->p_count, spreader->p_pos_f[0], spreader->p_pos_f[1], spreader->p_pos_f[2]); 
								SpecialPoint = spreader->p_pos_f; 
							} 
					}
			}
    } 

    // ================ FOR finalized cell event =========================
    else if (event == SP_FINALIZED_CELL) 
      { 
				if ( g_clarkson )
					{
						//-------- record which vertex was last finalized ------------
						//	printf("SP_FINALIZED_CELL ");

						//	last_point_in_region[num_regions_finalized] = grid_vertex_buffer_size-1;

						//	printf("%d\n", spreader->final_idx );

						drawRegion( spreader->final_idx, 0 );
						num_regions_finalized++;
						//------------------------------------------------------------
					}
				else // other finalization methods
					{
						if (delaunay2) 
							{ 
								spdelaunay2d->write_finalize_cell(spreader->final_idx); 
							} 
						else if (delaunay3) 
							{ 
								if (spreader->final_idx >= 0) 
									{ 
										spdelaunay3d->write_finalize_cell(spreader->final_idx); 
										//          fprintf(stderr,"finalizing %d (%.3f %.3f %.3f / %.3f %.3f %.3f) %d\n", spreader->final_idx, sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2], sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2], spreader->p_count); 
									} 
							} 
						else if (sscontainer2d) 
							{ 
								if (spreader->final_idx < 0) 
									{ 
										delete sscontainer2d; 
										sscontainer2d = 0; 
										sscontainer3d = new SScontainer3D(); 
										if (spreader->datatype == SP_DOUBLE) 
											sscontainer3d->open(spreader->bb_min_d, spreader->bb_max_d); 
										else 
											sscontainer3d->open(spreader->bb_min_f, spreader->bb_max_f); 
									} 
								else 
									{ 
										sscontainer2d->finalize_cell(spreader->final_idx);
										fprintf(stderr,"finalizing %d \n", spreader->final_idx);
									} 
							}         
						else if (sscontainer3d) 
							{ 
								sscontainer3d->finalize_cell(spreader->final_idx); 
								//        fprintf(stderr,"finalizing %d (%g %g %g / %g %g %g) %d\n", spreader->final_idx, sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2], sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2], spreader->p_count); 
							}
					} // end of other finalization methods
      } // end of finalized cell event
    
    if (spreader->p_count > NEXT_STEP) 
      { 
	NEXT_STEP += EVERY_NTH_STEP; 
	break; 
      } 
  } 
  
  if (event) 
    { 
    return 1; 
    } 
  else 
    { 
    grid_hash->clear(); 
    return 0; 
    } 
} 
 
void myReshape(int w, int h) 
{ 
  glutReshapeWindow(WindowW,WindowH); 
} 
 
void myIdle() 
{ 
  if (AnimationOn) 
  { 
    AnimationOn = vizContinue(); 
    if (!AnimationOn) 
    { 
      WorkingOn = 0; 
      vizEnd(); 
    } 
    glutPostRedisplay(); 
  } 
  else if (REPLAY_IT) 
  { 
    REPLAY_COUNT += NEXT_STEP; 
    glutPostRedisplay(); 
  } 
} 
 
void full_resolution_rendering() 
{ 
  if (file_name == 0) 
  { 
    fprintf(stderr,"ERROR: no input file\n"); 
  } 
 
  int p_count; 
 
  if (spreader) 
  { 
    p_count = spreader->p_count; 
    fprintf(stderr,"out-of-core rendering of %d points ... \n",p_count); 
    spreader->close(); 
    fclose(file); 
    delete spreader; 
    spreader = 0; 
  } 
  else 
  { 
    p_count = 2000000000; 
    fprintf(stderr,"out-of-core rendering of point cloud ... \n"); 
  } 
 
  if (strstr(file_name, ".sma") || strstr(file_name, ".obj") || strstr(file_name, ".smf")) 
  { 
    if (strstr(file_name, ".gz")) 
    { 
      #ifdef _WIN32 
      file = fopenGzipped(file_name, "r"); 
      #else 
      fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
      exit(1); 
      #endif 
    } 
    else 
    { 
      file = fopen(file_name, "r"); 
    } 
    if (file == 0) 
    { 
      fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
      exit(1); 
    } 
    SMreader_sma* smreader_sma = new SMreader_sma(); 
    smreader_sma->open(file); 
    smreader = smreader_sma; 
  } 
  else if (strstr(file_name, ".smb")) 
  { 
    if (strstr(file_name, ".gz")) 
    { 
      #ifdef _WIN32 
      file = fopenGzipped(file_name, "rb"); 
      #else 
      fprintf(stderr,"ERROR: cannot open gzipped file %s\n",file_name); 
      exit(1); 
      #endif 
    } 
    else 
    { 
      file = fopen(file_name, "rb"); 
    } 
    if (file == 0) 
    { 
      fprintf(stderr,"ERROR: cannot open %s\n",file_name); 
      exit(1); 
    } 
    SMreader_smb* smreader_smb = new SMreader_smb(); 
    smreader_smb->open(file); 
    smreader = smreader_smb; 
  } 
  else 
  { 
    fprintf(stderr,"ERROR: cannot guess input type from file name %s \n", file_name); 
    exit(1); 
  } 
 
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
   
  glViewport(0,0,WindowW,WindowH); 
 
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  gluPerspective(30.0f,(float)WindowW/WindowH,0.0625f,5.0f); 
 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0); 
 
  glRotatef(Elevation,1,0,0); 
  glRotatef(Azimuth,0,1,0); 
 
  glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ); 
  glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale); 
 
  glEnable(GL_DEPTH_TEST); 
  glEnable(GL_LIGHTING); 
  glEnable(GL_LIGHT0); 
  glEnable(GL_NORMALIZE); 
 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
 
  if (p_count == 0) 
  { 
    p_count = spreader->npoints; 
  } 
 
  glBegin(GL_POINTS); 
  while (spreader->p_count < p_count) 
  { 
    switch( (spreader->read_event()) )
    { 
    case SP_POINT: 
      glVertex3fv(spreader->p_pos_f); 
      break; 
    case SM_EOF: 
      glEnd(); 
      spreader->close(); 
      fclose(file); 
      delete spreader; 
      spreader = 0; 
      glutSwapBuffers(); 
      return; 
    default:
      break;

    } 
  } 
  glEnd(); 
  glutSwapBuffers(); 
} 
 
int pick_triangle(int x, int y) 
{ 
	// small viewport for picking 
 
  y = WindowH-y; 
  glViewport(x-5, y-5, 10, 10); 
 
	// set projection matrix 
 
	GLint vp[4]; 
	vp[0] = 0; 
	vp[1] = 0; 
	vp[2] = WindowW; 
	vp[3] = WindowH; 
 
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
	gluPickMatrix(x, y, 10, 10, vp);  
  gluPerspective(30.0f,(float)WindowW/WindowH,0.125f,10.0f); 
 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0); 
  glRotatef(Elevation,1,0,0); 
  glRotatef(Azimuth,0,1,0); 
  glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ); 
  glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale); 
 
  glClearColor(0,0,0,0); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
 
  float v0[3]; 
  float v1[3]; 
  float v2[3]; 
 
  glEnable(GL_DEPTH_TEST); 
  glDisable(GL_LIGHTING); 
  glDisable(GL_NORMALIZE); 
 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glBegin(GL_TRIANGLES); 
  spdelaunay2d->getActiveTrianglesInit(); 
  int idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
  while(idx != -1) 
  { 
    idx += 1; 
		glColor4ubv((const GLubyte *)&idx); 
    // draw 
    glVertex3fv(v0); 
    glVertex3fv(v1); 
    glVertex3fv(v2); 
    idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
  } 
  glEnd(); 
 
  glDisable(GL_DEPTH_TEST); 
 
  unsigned int id; 
 
	glReadPixels(x,y,1,1,GL_RGB,GL_UNSIGNED_BYTE,&id); 
 
  id = 0x00FFFFFF & id; 
 
//  fprintf(stderr, "at %d %d picked color %u\n", x, y, id); 
 
  if (id == 0) 
	{ 
		return -1; 
	} 
	else  
  { 
		return id - 1; 
	} 
} 
 
int pick_tetrahedron(int x, int y) 
{ 
	// small viewport for picking 
 
  y = WindowH-y; 
  glViewport(x-5, y-5, 10, 10); 
 
	// set projection matrix 
 
	GLint vp[4]; 
	vp[0] = 0; 
	vp[1] = 0; 
	vp[2] = WindowW; 
	vp[3] = WindowH; 
 
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
	gluPickMatrix(x, y, 10, 10, vp);  
  gluPerspective(30.0f,(float)WindowW/WindowH,0.125f,10.0f); 
 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0); 
  glRotatef(Elevation,1,0,0); 
  glRotatef(Azimuth,0,1,0); 
  glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ); 
  glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale); 
 
  glClearColor(0,0,0,0); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
 
  float v0[3]; 
  float v1[3]; 
  float v2[3]; 
  float v3[3]; 
 
  glEnable(GL_DEPTH_TEST); 
  glDisable(GL_LIGHTING); 
  glDisable(GL_NORMALIZE); 
 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glBegin(GL_TRIANGLES); 
  spdelaunay3d->getActiveTetrahedraInit(); 
  int idx = spdelaunay3d->getActiveTetrahedraNext(v0, v1, v2, v3); 
  while(idx != -1) 
  { 
    idx += 1; 
		glColor4ubv((const GLubyte *)&idx); 
    // draw 
    glVertex3fv(v0); 
    glVertex3fv(v2); 
    glVertex3fv(v1); 
    // draw 
    glVertex3fv(v0); 
    glVertex3fv(v1); 
    glVertex3fv(v3); 
    // draw 
    glVertex3fv(v2); 
    glVertex3fv(v3); 
    glVertex3fv(v1); 
    // draw 
    glVertex3fv(v0); 
    glVertex3fv(v3); 
    glVertex3fv(v2); 
    idx = spdelaunay3d->getActiveTetrahedraNext(v0, v1, v2, v3); 
  } 
  glEnd(); 
 
  glDisable(GL_DEPTH_TEST); 
 
  unsigned int id; 
 
	glReadPixels(x,y,1,1,GL_RGB,GL_UNSIGNED_BYTE,&id); 
 
  id = 0x00FFFFFF & id; 
 
//  fprintf(stderr, "at %d %d picked color %u\n", x, y, id); 
 
  if (id == 0) 
    { 
      return -1; 
    } 
  else  
    { 
      return id - 1; 
    } 
} 
 
void myMouseFunc(int button, int state, int x, int y) 
{ 
  NewX=x; NewY=y; 
     
  if (button == GLUT_LEFT_BUTTON) 
  { 
    LeftButtonDown = !state; 
    MiddleButtonDown = 0; 
    RightButtonDown = 0; 
  } 
  else if (button == GLUT_RIGHT_BUTTON) 
  { 
    LeftButtonDown = 0; 
    MiddleButtonDown = 0; 
    RightButtonDown = !state; 
  } 
  else if (button == GLUT_MIDDLE_BUTTON) 
  { 
    LeftButtonDown = 0; 
    MiddleButtonDown = !state; 
    RightButtonDown = 0; 
  } 
} 
 
void myMotionFunc(int x, int y) 
{ 
  OldX=NewX; OldY=NewY; 
  NewX=x;    NewY=y; 
   
  float RelX = (NewX-OldX) / (float)glutGet((GLenum)GLUT_WINDOW_WIDTH); 
  float RelY = (NewY-OldY) / (float)glutGet((GLenum)GLUT_WINDOW_HEIGHT); 
  if (LeftButtonDown)  
  {  
    if (InteractionMode == 0) 
    { 
      Azimuth += (RelX*180); 
      Elevation += (RelY*180); 
    } 
    else if (InteractionMode == 1) 
    { 
      DistX-=RelX; 
      DistY+=RelY; 
    } 
    else if (InteractionMode == 2) 
    { 
      DistZ-=RelY*DistZ; 
    } 
  } 
  else if (MiddleButtonDown) 
  { 
    DistX-=RelX*1.0f; 
    DistY+=RelY*1.0f; 
  } 
 
  glutPostRedisplay(); 
} 
 
void MyMenuFunc(int value); 
 
void myKeyboard(unsigned char Key, int x, int y) 
{ 
  switch(Key) 
  { 
  case 'Q': 
  case 'q': 
  case 27: 
		// if we are using clarkson
		if( g_clarkson ) {
			// deallocate before exiting
			free(clarkson2d_mesh_verts);
			free(clarkson2d_mesh_faces);
			//  free(last_point_in_region);
			//  destroyGridVertexBuffer();
	  }
    exit(0); 
    break; 
  case ' ': // rotate, translate, or zoom 
    if (InteractionMode == 2) 
    { 
      InteractionMode = 0; 
    } 
    else 
    { 
      InteractionMode += 1; 
    } 
    glutPostRedisplay(); 
    break; 
  case '>': //zoom in 
    DistZ-=0.1f; 
    break; 
  case '<': //zoom out 
    DistZ+=0.1f; 
    break; 
  case '_': 
    LINE_WIDTH -= 1; 
    if (LINE_WIDTH < 0) LINE_WIDTH = 0; 
    fprintf(stderr,"LINE_WIDTH %d\n",LINE_WIDTH); 
    glutPostRedisplay(); 
    break; 
  case '+': 
    LINE_WIDTH += 1; 
    fprintf(stderr,"LINE_WIDTH %d\n",LINE_WIDTH); 
    glutPostRedisplay(); 
    break; 
  case '-': 
    POINT_SIZE -= 1; 
    if (POINT_SIZE < 0) POINT_SIZE = 0; 
    fprintf(stderr,"POINT_SIZE %d\n",POINT_SIZE); 
    glutPostRedisplay(); 
    break; 
  case '=': 
    POINT_SIZE += 1; 
    fprintf(stderr,"POINT_SIZE %d\n",POINT_SIZE); 
    glutPostRedisplay(); 
    break; 
  case 'B': 
  case 'b': 
    RENDER_BOUNDINGBOX = (RENDER_BOUNDINGBOX+1)%3; 
    fprintf(stderr,"RENDER_BOUNDINGBOX %d\n",RENDER_BOUNDINGBOX); 
    glutPostRedisplay(); 
    break; 
  case 'C': 
  case 'c': 
    RENDER_CIRCUMSPHIRCLES = (RENDER_CIRCUMSPHIRCLES+1)%3; 
    fprintf(stderr,"RENDER_CIRCUMSPHIRCLES %d\n",RENDER_CIRCUMSPHIRCLES); 
    glutPostRedisplay(); 
    break; 
  case 'M': 
  case 'm': 
    if (spdelaunay2d)  
    { 
      SpecialElement = pick_triangle(x, y); 
      if (SpecialElement != -1) 
      { 
        float cen_rad[4]; 
        if (spdelaunay2d->getActiveTrianglesIdx(SpecialElement, cen_rad)) 
        { 
          fprintf(stderr,"picked %d with center (%g,%g) and radius %g\n",SpecialElement,cen_rad[0],cen_rad[1],cen_rad[3]); 
        } 
        else 
        { 
          fprintf(stderr,"picked %d which has no circumcircle computed yet\n",SpecialElement); 
        } 
      } 
      glutPostRedisplay(); 
    } 
    else if (spdelaunay3d)  
    { 
      SpecialElement = pick_tetrahedron(x, y); 
      glutPostRedisplay(); 
    } 
    break; 
  case 'O': 
  case 'o': 
    AnimationOn = 0; 
    REPLAY_IT = 0; 
    break; 
  case 'U': 
    RENDER_MODE_UNFINALIZED = 0; 
    fprintf(stderr,"RENDER_MODE_UNFINALIZED %d\n",RENDER_MODE_UNFINALIZED); 
    glutPostRedisplay(); 
    break; 
  case 'u': 
    RENDER_MODE_UNFINALIZED = RENDER_MODE_UNFINALIZED + 1; 
    if (RENDER_MODE_UNFINALIZED > 3) RENDER_MODE_UNFINALIZED = 1; 
    fprintf(stderr,"RENDER_MODE_UNFINALIZED %d\n",RENDER_MODE_UNFINALIZED); 
    glutPostRedisplay(); 
    break; 
  case 'A': 
    RENDER_MODE_ACTIVE = 0; 
    fprintf(stderr,"RENDER_MODE_ACTIVE %d\n",RENDER_MODE_ACTIVE); 
    glutPostRedisplay(); 
    break; 
  case 'a': 
    RENDER_MODE_ACTIVE = RENDER_MODE_ACTIVE + 1; 
    if (RENDER_MODE_ACTIVE > 3) RENDER_MODE_ACTIVE = 1; 
    fprintf(stderr,"RENDER_MODE_ACTIVE %d\n",RENDER_MODE_ACTIVE); 
    glutPostRedisplay(); 
    break; 
  case 'F': 
  case 'f': 
    break; 
  case 'I': 
    RENDER_MODE_INFINITE = 0; 
    fprintf(stderr,"RENDER_MODE_INFINITE %d\n",RENDER_MODE_INFINITE); 
    glutPostRedisplay(); 
    break; 
  case 'i': 
    RENDER_MODE_INFINITE = RENDER_MODE_INFINITE + 1; 
    if (RENDER_MODE_INFINITE > 2) RENDER_MODE_INFINITE = 1; 
    fprintf(stderr,"RENDER_MODE_INFINITE %d\n",RENDER_MODE_INFINITE); 
    glutPostRedisplay(); 
    break; 
  case 'V': 
  case 'v': 
    STREAM_COLORING = STREAM_COLORING + 1; 
    if (STREAM_COLORING > 3) STREAM_COLORING = 1; 
    fprintf(stderr,"STREAM_COLORING %d\n",STREAM_COLORING); 
    glutPostRedisplay(); 
    break; 
  case 'R': 
  case 'r': 
    full_resolution_rendering(); 
    break; 
  case 'W': 
  case 'w': 
    break; 
  case 'D': 
  case 'd': 
    PrintToFileOn=1; 
    if (Framebuffer) delete [] Framebuffer; 
    Framebuffer = new unsigned char[WindowW*WindowH*3]; 
    fprintf(stderr,"print_to_file %d\n",PrintToFileOn); 
    break; 
  case 'T': 
    if (DIRTY_MESH) 
    { 
      // works only in replay mode 
      fprintf(stderr,"tiny steps only work during second play (replay)\n"); 
    } 
    else 
    { 
      REPLAY_COUNT -= 1; 
      if (REPLAY_COUNT < 0) 
      { 
        REPLAY_COUNT = 0; 
      } 
    } 
    glutPostRedisplay(); 
    break; 
  case 't': 
    if (DIRTY_MESH) 
    { 
      // works only in replay mode 
      fprintf(stderr,"tiny steps only work during second play (replay)\n"); 
    } 
    else 
    { 
      if (REPLAY_COUNT >= triangle_buffer_size) 
      { 
        REPLAY_COUNT = 0; 
      } 
      REPLAY_COUNT += 1; 
    } 
    glutPostRedisplay(); 
    break; 
  case 'S': 
    if (DIRTY_MESH) 
    { 
      // works only in replay mode 
      fprintf(stderr,"back stepping only work during second play (replay)\n"); 
    } 
    else 
    { 
      NEXT_STEP = triangle_buffer_size / EXACTLY_N_STEPS; 
      if (NEXT_STEP == 0) NEXT_STEP = 1; 
      REPLAY_COUNT -= NEXT_STEP; 
      if (REPLAY_COUNT < 0) 
      { 
        REPLAY_COUNT = 0; 
      } 
    } 
    glutPostRedisplay(); 
    break; 
  case 'P': 
  case 'p': 
    DIRTY_MESH = 1; 
    if (DIRTY_MESH) 
    { 
      AnimationOn = !AnimationOn; 
    } 
    else 
    { 
      if (REPLAY_IT == 0) 
      { 
        if (REPLAY_COUNT >= grid_vertex_buffer_size) 
        { 
          REPLAY_COUNT = 0; 
        } 
        NEXT_STEP = grid_vertex_buffer_size / EXACTLY_N_STEPS; 
        if (NEXT_STEP == 0) NEXT_STEP = 1; 
        REPLAY_IT = 1; 
      } 
      else 
      { 
        REPLAY_IT = 0; 
      } 
    } 
  case 's': 
    if (DIRTY_MESH) 
    { 
      if (WorkingOn == 0) 
      { 

	// ---- reading sma ---------------
	// instead of converting it to point clouds via SPconverter
	// and for the purpose of only drawing the coarse triangle
	// we should construct a list of triangles with (x, y, z) coordinates
				
				if( g_clarkson ) {
					read(); 
				} else {
					vizBegin();
				}

	// draw fine points
        WorkingOn = vizContinue(); 
      } 
      else 
      { 
	// draw fine points
        WorkingOn = vizContinue(); 
      } 
      if (WorkingOn == 0) 
      { 
	vizEnd();
        AnimationOn = 0; 
        PrintToFileOn = 0; 
      } 
    } 
    else 
    { 
      if (REPLAY_COUNT >= grid_vertex_buffer_size) 
      { 
        REPLAY_COUNT = 0; 
      } 
      NEXT_STEP = grid_vertex_buffer_size / EXACTLY_N_STEPS; 
      if (NEXT_STEP == 0) NEXT_STEP = 1; 
      REPLAY_COUNT += NEXT_STEP; 
    } 
    glutPostRedisplay(); 
    break; 
  case 'K': 
  case 'k': 
    printf("Azimuth = %gf;\n",Azimuth); 
    printf("Elevation = %gf;\n",Elevation); 
    printf("DistX = %gf; DistY = %gf; DistZ = %gf;\n",DistX,DistY,DistZ); 
    break; 
  case 'L': 
  case 'l': 
  case '1': 
    MyMenuFunc(1); 
    break; 
  case '2': 
    MyMenuFunc(2); 
    break; 
  case '3': 
    MyMenuFunc(3); 
    break; 
  case '4': 
    MyMenuFunc(4); 
    break; 
  case '5': 
    MyMenuFunc(5); 
    break; 
  case '6': 
    MyMenuFunc(6); 
    break; 
  case '7': 
    MyMenuFunc(7); 
    break; 
  case '8': 
    MyMenuFunc(8); 
    break; 
  case '9': 
    MyMenuFunc(9); 
    break; 
  case '0': 
    MyMenuFunc(10); 
    break; 
  }; 
} 
 
void MyMenuFunc(int value) 
{ 
  if (value >= 100) 
  { 
    if (value <= 102) 
    { 
      InteractionMode = value - 100; 
      glutPostRedisplay(); 
    } 
    else if (value == 103) 
    { 
      myKeyboard('s',0,0); 
    } 
    else if (value == 104) 
    { 
      myKeyboard('p',0,0); 
    } 
    else if (value == 105) 
    { 
      myKeyboard('o',0,0); 
    } 
    else if (value == 109) 
    { 
      myKeyboard('q',0,0); 
    } 
    else if (value == 150) 
    { 
      myKeyboard('v',0,0); 
    } 
    else if (value == 151) 
    { 
      myKeyboard('l',0,0); 
    } 
    else if (value == 152) 
    { 
      myKeyboard('c',0,0); 
    } 
    else if (value == 153) 
    { 
      myKeyboard('m',0,0); 
    } 
  } 
  else if (value == 40) 
  { 
    EXACTLY_N_STEPS = 5; 
  } 
  else if (value == 41) 
  { 
    EXACTLY_N_STEPS = 10; 
  } 
  else if (value == 42) 
  { 
    EXACTLY_N_STEPS = 25; 
  } 
  else if (value == 43) 
  { 
    EXACTLY_N_STEPS = 50; 
  } 
  else if (value == 44) 
  { 
    EXACTLY_N_STEPS = 100; 
  } 
  else if (value == 45) 
  { 
    EXACTLY_N_STEPS = 250; 
  } 
  else if (value == 46) 
  { 
    EXACTLY_N_STEPS = 500; 
  } 
  else if (value == 47) 
  { 
    EXACTLY_N_STEPS = 1000; 
  } 
  else if (value == 48) 
  { 
    EXACTLY_N_STEPS = 10000; 
  } 
  else if (value == 71) 
  { 
    if (GRID_PRECISION != 3) DIRTY_MESH = 1; 
    GRID_PRECISION = 3; 
  } 
  else if (value == 72) 
  { 
    if (GRID_PRECISION != 4) DIRTY_MESH = 1; 
    GRID_PRECISION = 4; 
  } 
  else if (value == 73) 
  { 
    if (GRID_PRECISION != 5) DIRTY_MESH = 1; 
    GRID_PRECISION = 5; 
  } 
  else if (value == 74) 
  { 
    if (GRID_PRECISION != 6) DIRTY_MESH = 1; 
    GRID_PRECISION = 6; 
  } 
  else if (value == 75) 
  { 
    if (GRID_PRECISION != 7) DIRTY_MESH = 1; 
    GRID_PRECISION = 7; 
  } 
  else if (value == 76) 
  { 
    if (GRID_PRECISION != 8) DIRTY_MESH = 1; 
    GRID_PRECISION = 8; 
  } 
  else if (value == 77) 
  { 
    if (GRID_PRECISION != 9) DIRTY_MESH = 1; 
    GRID_PRECISION = 9; 
  } 
  else if (value == 78) 
  { 
    if (GRID_PRECISION != 10) DIRTY_MESH = 1; 
    GRID_PRECISION = 10; 
  } 
} 
 
void displayMessage( ) 
{ 
  glColor3f( 0.7f, 0.7f, 0.7f );  // Set colour to grey 
  glMatrixMode( GL_PROJECTION ); 
  glPushMatrix(); 
  glLoadIdentity(); 
  gluOrtho2D( 0.0f, 1.0f, 0.0f, 1.0f ); 
  glMatrixMode( GL_MODELVIEW ); 
  glPushMatrix(); 
  glLoadIdentity(); 
  glRasterPos2f( 0.03f, 0.95f ); 
   
  if( InteractionMode == 0 ) 
  { 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'r'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'o'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 't'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'a'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 't'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'e'); 
  } 
  else if( InteractionMode == 1 ) 
  {     
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 't'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'r'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'a'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'n'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 's'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'l'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'a'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 't'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'e'); 
  } 
  else 
  { 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'z'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'o'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'o'); 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'm'); 
  } 
   
  glPopMatrix(); 
  glMatrixMode( GL_PROJECTION ); 
  glPopMatrix(); 
} 
 
void render_triangle(GridVertex* vertex0, GridVertex* vertex1, GridVertex* vertex2) 
{ 
  float n[3]; 
  float v0[3]; 
  float v1[3]; 
  float v2[3]; 
 
  VecScalarDiv3fv(v0, vertex0->v, (float)vertex0->number); 
  VecScalarDiv3fv(v1, vertex1->v, (float)vertex1->number); 
  VecScalarDiv3fv(v2, vertex2->v, (float)vertex2->number); 
  VecCcwNormal3fv(n, v0, v1, v2); 
 
  glNormal3fv(n); 
  glVertex3fv(v0); 
  glVertex3fv(v1); 
  glVertex3fv(v2); 
} 
 
void render_point(GridVertex* vertex) 
{ 
  float p[3]; 
 
  VecScalarDiv3fv(p, vertex->v, (float)vertex->number); 
 
  glVertex3fv(p); 
} 
 
void render_line(GridVertex* vertex0, GridVertex* vertex1) 
{ 
  float v0[3]; 
  float v1[3]; 
 
  VecScalarDiv3fv(v0, vertex0->v, (float)vertex0->number); 
  VecScalarDiv3fv(v1, vertex1->v, (float)vertex1->number); 
 
  glVertex3fv(v0); 
  glVertex3fv(v1); 
} 
 
void myDisplay() 
{ 
	if( !g_clarkson )
		{
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
		}


  glViewport(0,0,WindowW,WindowH); 
 
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  gluPerspective(30.0f,(float)WindowW/WindowH,0.0625f,5.0f); 
 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0); 
 
  glRotatef(Elevation,1,0,0); 
  glRotatef(Azimuth,0,1,0); 
 
  glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ); 
  glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale); 
 
  glEnable(GL_DEPTH_TEST); 
  glEnable(GL_LIGHTING); 
  glEnable(GL_LIGHT0); 
  glEnable(GL_NORMALIZE); 
 
  int rendered_points; 

	if( g_clarkson )
	{
  //  rendered_points = grid_vertex_buffer_size;
  //============================ display =============================


  /*
  glDisable(GL_LIGHTING); 
  glDisable(GL_LIGHT0); 
  glDisable(GL_NORMALIZE); 
  glDisable(GL_DEPTH_TEST);


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); // clean the matrix
  glOrtho(-1, 1, -1, 1, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); // clean the matrix
  */
  //  printf("f=%d",clarkson2d_mesh_nfaces );


  //  glPointSize(POINT_SIZE); 

  //  glBegin(GL_POINTS); 
  //  float color = 0;
  // note that i = 0 includes pinf
  /*
  printf("num_regions_finalized = %d\n", num_regions_finalized);
  for(int i=0;i< num_regions_finalized; i++)
    printf("%d ",last_point_in_region[i]);
  printf("\n");
  */
  /*
  for (int i = 0, region=0; i < rendered_points; i++) 
    { 

      //      color = ((float)region + 1) / (num_regions_finalized+1);
      //      if( color>0.5 )
      //printf("%.3f ", color);
      glColor3f(color, 0, 0);
      glVertex3fv( grid_vertex_buffer[i].v );
      // if this pt is the last one finalized by this region
      // increment the color so that next point has diff color
      while(i==last_point_in_region[region]) {
	region++;
	color+=inc;
      }
    } 
  glEnd();
  */
  //  printf("\n");

  //============  end of clarson2D display ==============================
	}
 
  if (DIRTY_MESH) 
  { 
    rendered_points = grid_vertex_buffer_size; 
  } 
  else 
  { 
    if (REPLAY_COUNT > grid_vertex_buffer_size) 
    { 
      rendered_points = grid_vertex_buffer_size; 
      REPLAY_IT = 0; 
    } 
    else 
    { 
      rendered_points = REPLAY_COUNT; 
    } 
  } 
 
  //  
  glLineWidth(LINE_WIDTH); 
 
  // draw streaming delaunay state 
 
  if (spdelaunay2d) 
  { 
    float v0[3]; 
    float v1[3]; 
    float v2[3]; 
 
    // draw active triangles 
 
    if (RENDER_MODE_ACTIVE) 
    { 
      if (RENDER_MODE_ACTIVE > 1) 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[5]); 
      glBegin(GL_TRIANGLES); 
      glNormal3f(0,0,1); 
      spdelaunay2d->getActiveTrianglesInit(); 
      int idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
      while(idx != -1) 
      { 
        if (idx == SpecialElement) 
        { 
          idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
          continue; 
        } 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v2); 
        idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
      } 
      glEnd(); 
 
 
      // maybe draw active triangle wireframe on top 
 
      if (RENDER_MODE_ACTIVE > 2) 
      { 
        glLineWidth(1.0f); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[0]); 
        glBegin(GL_TRIANGLES); 
        glNormal3f(0,0,1); 
        spdelaunay2d->getActiveTrianglesInit(); 
        int idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
        while(idx != -1) 
        { 
          glVertex3fv(v0); 
          glVertex3fv(v1); 
          glVertex3fv(v2); 
          idx = spdelaunay2d->getActiveTrianglesNext(v0, v1, v2); 
        } 
        glEnd(); 
        glLineWidth(LINE_WIDTH); 
      } 
    } 
 
    // maybe draw special triangle 
 
    if (SpecialElement != -1) 
    { 
      if (spdelaunay2d->getActiveTrianglesIdx(SpecialElement, v0, v1, v2)) 
      { 
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[1]); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
        glBegin(GL_TRIANGLES); 
        // draw 
        glBegin(GL_TRIANGLES); 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v2); 
        glEnd(); 
      } 
    } 
 
    // maybe draw circumcircles 
 
    if (RENDER_CIRCUMSPHIRCLES) //TODO: change back... jack
    { 
      if (SpecialElement != -1) 
      { 
        float cen_rad[4]; 
        if (spdelaunay2d->getActiveTrianglesIdx(SpecialElement, cen_rad)) 
        { 
          glEnable(GL_CULL_FACE); 
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[4]); 
          glTranslatef(cen_rad[0], cen_rad[1], cen_rad[2]); 
          if (RENDER_CIRCUMSPHIRCLES == 1) 
          { 
            glutWireSphere(cen_rad[3], 64, 2); 
          } 
          else 
          { 
            glutSolidSphere(cen_rad[3], 64, 2); 
          } 
          glDisable(GL_CULL_FACE); 
          glTranslatef(-cen_rad[0], -cen_rad[1], -cen_rad[2]); 
        } 
      } 
    } 
 
    // draw infinite triangles 
 
    if (RENDER_MODE_INFINITE) 
    { 
      if (RENDER_MODE_INFINITE == 1) 
      { 
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
        glEnable(GL_BLEND); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[9]); 
      glBegin(GL_TRIANGLES); 
      glNormal3f(0,0,1); 
      spdelaunay2d->getInfiniteTrianglesInit(); 
      while(spdelaunay2d->getInfiniteTrianglesNext(v0, v1, v2)) 
      { 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v2); 
      } 
      glEnd(); 
      glDisable(GL_BLEND); 
    } 
  } 
  else if (spdelaunay3d) 
  { 
    float v0[3]; 
    float v1[3]; 
    float v2[3]; 
    float v3[3]; 
    float normal[3]; 
 
    // draw active tetrahedra 
 
    if (RENDER_MODE_ACTIVE) 
    { 
      if (RENDER_MODE_ACTIVE == 1) 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[5]); 
      glBegin(GL_TRIANGLES); 
      spdelaunay3d->getActiveTetrahedraInit(); 
      int idx = spdelaunay3d->getActiveTetrahedraNext(v0, v1, v2, v3); 
      while (idx != -1) 
      { 
        if (idx == SpecialElement) 
        { 
          idx = spdelaunay3d->getActiveTetrahedraNext(v0, v1, v2, v3); 
          continue; 
        } 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v2,v1); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v2); 
        glVertex3fv(v1); 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v1,v3); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v3); 
        // draw 
        VecCcwNormNormal3fv(normal,v2,v3,v1); 
        glNormal3fv(normal); 
        glVertex3fv(v2); 
        glVertex3fv(v3); 
        glVertex3fv(v1); 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v3,v2); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v3); 
        glVertex3fv(v2); 
        idx = spdelaunay3d->getActiveTetrahedraNext(v0, v1, v2, v3); 
      } 
      glEnd(); 
    } 
 
    // maybe draw special tetrahedron 
 
    if (SpecialElement != -1) 
    { 
      if (spdelaunay3d->getActiveTetrahedraCoords(SpecialElement, v0, v1, v2, v3)) 
      { 
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[1]); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
        glBegin(GL_TRIANGLES); 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v2,v1); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v2); 
        glVertex3fv(v1); 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v1,v3); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v3); 
        // draw 
        VecCcwNormNormal3fv(normal,v2,v3,v1); 
        glNormal3fv(normal); 
        glVertex3fv(v2); 
        glVertex3fv(v3); 
        glVertex3fv(v1); 
        // draw 
        VecCcwNormNormal3fv(normal,v0,v3,v2); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v3); 
        glVertex3fv(v2); 
        glEnd(); 
      } 
    } 
 
    // maybe draw circumcircles 
 
    if (RENDER_CIRCUMSPHIRCLES) 
    { 
      if (SpecialElement != -1) 
      { 
        float cen_rad[4]; 
        if (spdelaunay3d->getActiveTetrahedraSphere(SpecialElement, cen_rad)) 
        { 
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[4]); 
          glTranslatef(cen_rad[0], cen_rad[1], cen_rad[2]); 
          if (RENDER_CIRCUMSPHIRCLES == 1) 
          { 
            glutWireSphere(cen_rad[3], 64, 64); 
          } 
          else 
          { 
            glutSolidSphere(cen_rad[3], 64, 64); 
          } 
          glTranslatef(-cen_rad[0], -cen_rad[1], -cen_rad[2]); 
        } 
      } 
    } 
 
    // maybe draw the grid cell that stores it 
    if (RENDER_BOUNDINGBOX) 
    { 
      if (SpecialElement != -1) 
      { 
        int cell_idx; 
        if (spdelaunay3d->getActiveTetrahedraCellIdx(SpecialElement, &cell_idx)) 
        { 
          if (sscontainer3d->noiterate(cell_idx)) 
          { 
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[6]); 
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
            glBegin(GL_QUADS); 
            glNormal3f(0,0,-1); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
 
            glNormal3f(0,-1,0); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
 
            glNormal3f(-1,0,0); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
 
            glNormal3f(0,0,1); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
 
            glNormal3f(0,1,0); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
 
            glNormal3f(1,0,0); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
            glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
            glEnd(); 
          } 
        } 
      } 
    } 
 
    // draw infinite triangles 
 
    if (RENDER_MODE_INFINITE) 
    { 
      if (RENDER_MODE_INFINITE == 1) 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[9]); 
      glBegin(GL_TRIANGLES); 
      spdelaunay3d->getInfiniteTetrahedraBaseTriangleInit(); 
      while(spdelaunay3d->getInfiniteTetrahedraBaseTriangleNext(v0, v1, v2)) 
      { 
        VecCcwNormNormal3fv(normal,v0,v1,v2); 
        glNormal3fv(normal); 
        glVertex3fv(v0); 
        glVertex3fv(v1); 
        glVertex3fv(v2); 
      } 
      glEnd(); 
    } 
  } 
 
  // draw unfinalized space 
 
  if (sscontainer2d) 
  { 
    glDisable(GL_DEPTH_TEST); 
    if (RENDER_MODE_UNFINALIZED) 
    { 
      if (RENDER_MODE_UNFINALIZED > 1) 
      { 
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
        glEnable(GL_BLEND); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[6]); 
      glBegin(GL_QUADS); 
      glNormal3f(0,0,-1); 
      sscontainer2d->iterateInit(); 
      while(sscontainer2d->iterateNext()) 
      { 
        glVertex3f(sscontainer2d->r_min_f[0], sscontainer2d->r_min_f[1], spreader->bb_min_f[2]); 
        glVertex3f(sscontainer2d->r_min_f[0], sscontainer2d->r_max_f[1], spreader->bb_min_f[2]); 
        glVertex3f(sscontainer2d->r_max_f[0], sscontainer2d->r_max_f[1], spreader->bb_min_f[2]); 
        glVertex3f(sscontainer2d->r_max_f[0], sscontainer2d->r_min_f[1], spreader->bb_min_f[2]); 
      } 
      glEnd(); 
      glDisable(GL_BLEND); 
 
      // maybe draw wireframe on top 
 
      if (RENDER_MODE_UNFINALIZED > 2) 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[6]); 
        glBegin(GL_QUADS); 
        glNormal3f(0,0,-1); 
        sscontainer2d->iterateInit(); 
        while(sscontainer2d->iterateNext()) 
        { 
          glVertex3f(sscontainer2d->r_min_f[0], sscontainer2d->r_min_f[1], spreader->bb_min_f[2]); 
          glVertex3f(sscontainer2d->r_min_f[0], sscontainer2d->r_max_f[1], spreader->bb_min_f[2]); 
          glVertex3f(sscontainer2d->r_max_f[0], sscontainer2d->r_max_f[1], spreader->bb_min_f[2]); 
          glVertex3f(sscontainer2d->r_max_f[0], sscontainer2d->r_min_f[1], spreader->bb_min_f[2]); 
        } 
        glEnd(); 
      } 
    } 
  } 
  else if (sscontainer3d) 
  { 
    if (RENDER_MODE_UNFINALIZED) 
    { 
      if (RENDER_MODE_UNFINALIZED == 1) 
      { 
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
        glEnable(GL_BLEND); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
      } 
      else 
      { 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
      } 
 
      glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_diffuse[6]); 
      glMaterialfv(GL_BACK, GL_DIFFUSE, colours_diffuse[6]); 
      glBegin(GL_QUADS); 
      sscontainer3d->iterateInit(); 
      while(sscontainer3d->iterateNext()) 
      { 
        glNormal3f(0,0,-1); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
 
        glNormal3f(0,-1,0); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
 
        glNormal3f(-1,0,0); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
 
        glNormal3f(0,0,1); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
 
        glNormal3f(0,1,0); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_min_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
 
        glNormal3f(1,0,0); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_min_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_max_f[1], sscontainer3d->r_max_f[2]); 
        glVertex3f(sscontainer3d->r_max_f[0], sscontainer3d->r_min_f[1], sscontainer3d->r_max_f[2]); 
      } 
      glEnd(); 
      glDisable(GL_BLEND); 
    } 
  } 
 
  glDisable(GL_LIGHTING); 
  glDisable(GL_LIGHT0); 
  glDisable(GL_NORMALIZE); 
  glDisable(GL_DEPTH_TEST); 
 
  // draw points 
 
  if (grid_vertex_buffer_size && POINT_SIZE) 
  { 
    glPointSize(POINT_SIZE); 
    glBegin(GL_POINTS); 
    if (STREAM_COLORING == 0) 
    { 
      glColor3f(0,0,0); 
      for (int i = 0; i < rendered_points; i++) 
      { 
        render_point(&(grid_vertex_buffer[i])); 
      } 
    } 
    else if (STREAM_COLORING == 1) 
    { 
      float color; 
      for (int i = 0; i < rendered_points; i++) 
      { 
        color = 0.1f+0.7f*i/grid_vertex_buffer_size; 
        glColor3f(color,color,color); 
        render_point(&(grid_vertex_buffer[i])); 
      } 
    } 
    else if (STREAM_COLORING == 2) 
    { 
      for (int i = 0; i < rendered_points; i++) 
      { 
        if (i < grid_vertex_buffer_size/3) 
        { 
          glColor3f(0.1f+0.7f*i/(grid_vertex_buffer_size/3),0.1f,0.1f); 
        } 
        else if (i < 2*(grid_vertex_buffer_size/3)) 
        { 
          glColor3f(0.8f,0.1f+0.7f*(i-(grid_vertex_buffer_size/3))/(grid_vertex_buffer_size/3),0.1f); 
        } 
        else 
        { 
          glColor3f(0.8f, 0.8f, 0.1f+0.7f*(i-2*(grid_vertex_buffer_size/3))/(grid_vertex_buffer_size/3)); 
        } 
        render_point(&(grid_vertex_buffer[i])); 
      } 
    } 
    else 
    { 
      for (int i = 0; i < rendered_points; i++) 
      { 
        if (i % 3 == 0) 
        { 
          glColor3f(0,0,0); 
        } 
        else 
        { 
          if (i < grid_vertex_buffer_size/3) 
          { 
            glColor3f(0.1f+0.7f*i/(grid_vertex_buffer_size/3),0.1f,0.1f); 
          } 
          else if (i < 2*(grid_vertex_buffer_size/3)) 
          { 
            glColor3f(0.8f,0.1f+0.7f*(i-(grid_vertex_buffer_size/3))/(grid_vertex_buffer_size/3),0.1f); 
          } 
          else 
          { 
            glColor3f(0.8f, 0.8f, 0.1f+0.7f*(i-2*(grid_vertex_buffer_size/3))/(grid_vertex_buffer_size/3)); 
          } 
        } 
        render_point(&(grid_vertex_buffer[i])); 
      } 
    } 
    glEnd(); 
  } 
 
  if (SpecialPoint) 
  { 
    glColor3f(1,0,0); 
    glPointSize(2.0f); 
    glBegin(GL_POINTS); 
    glVertex3fv(SpecialPoint); 
    glEnd(); 
  } 
 
  if (RENDER_BOUNDINGBOX) 
  { 
    glColor3f(0,0,0); 
    glLineWidth(1.0f); 
    glBegin(GL_QUADS); 
    if (RENDER_BOUNDINGBOX == 1) 
    { 
      glVertex3f(boundingBoxMin[0], boundingBoxMin[1], boundingBoxMin[2]); 
      glVertex3f(boundingBoxMin[0], boundingBoxMax[1], boundingBoxMin[2]); 
      glVertex3f(boundingBoxMax[0], boundingBoxMax[1], boundingBoxMin[2]); 
      glVertex3f(boundingBoxMax[0], boundingBoxMin[1], boundingBoxMin[2]); 
    } 
    else 
    { 
      glVertex3f(boundingBoxMin[0], boundingBoxMin[1], boundingBoxMax[2]); 
      glVertex3f(boundingBoxMin[0], boundingBoxMax[1], boundingBoxMax[2]); 
      glVertex3f(boundingBoxMax[0], boundingBoxMax[1], boundingBoxMax[2]); 
      glVertex3f(boundingBoxMax[0], boundingBoxMin[1], boundingBoxMax[2]); 
    } 
    glEnd(); 
    glLineWidth(1.0f); 
  } 
 
  if (!PrintToFileOn) 
  { 
    displayMessage(); 
  } 
  else 
  { 
    char FileName[256], Command[256]; 
    sprintf(FileName, "./temp.ppm"); 
    glPixelStorei(GL_PACK_ALIGNMENT, 1); 
    #ifdef _WIN32 
    glReadPixels( 0, 0, WindowW, WindowH, GL_RGB, GL_UNSIGNED_BYTE, Framebuffer ); 
    #else 
    glReadPixels( (1024-WindowW)/2, (1024-WindowH)/2, WindowW, WindowH, GL_RGB, GL_UNSIGNED_BYTE, Framebuffer ); 
    #endif 
    SavePPM(FileName, Framebuffer, WindowW, WindowH); 
    #ifdef _WIN32 
      sprintf(Command, "i_view32.exe temp.ppm /convert=%s%d%d%d%d.jpg",PrintFileName,((Time/1000)%10),((Time/100)%10),((Time/10)%10),(Time%10)); 
    #else 
      sprintf(Command, "convert ./temp.ppm ./%s%d%d%d%d.jpg",PrintFileName,((Time/1000)%10),((Time/100)%10),((Time/10)%10),(Time%10)); 
//      sprintf(Command, "imgcopy -q%d ./temp.ppm ./pics/%s%d%d%d%d.jpg",PrintQuality,PrintFileName,((Time/1000)%10),((Time/100)%10),((Time/10)%10),(Time%10)); 
    #endif 
    printf("performing: '%s'\n", Command); 
    system(Command); 
    Time++; 
  } 
  //  glutSwapBuffers(); 
} 
 
int main(int argc, char *argv[]) 
{ 
  if (argc == 1 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"-help") == 0)  
  { 
    usage(); 
    exit(0); 
  } 
  else 
  { 
    for (int i = 1; i < argc; i++) 
    { 
      if (strcmp(argv[i],"-win") == 0) 
      { 
        i++; 
        WindowW = atoi(argv[i]); 
        i++; 
        WindowH = atoi(argv[i]); 
      } 
      else if (strcmp(argv[i],"-del2") == 0) 
      { 
        delaunay2 = true; 
      } 
      else if (strcmp(argv[i],"-del3") == 0) 
      { 
        delaunay3 = true; 
      } 
      else if (strcmp(argv[i],"-flatten") == 0) 
      { 
        flatten = true; 
      } 
      else if (strcmp(argv[i],"-steps") == 0) 
      { 
        i++; 
        EXACTLY_N_STEPS = atoi(argv[i]);; 
      } 
      else if (strcmp(argv[i],"-every") == 0) 
      { 
        i++; 
        EVERY_NTH_STEP = atoi(argv[i]);; 
      } 
      else if (strcmp(argv[i],"-scatter") == 0) 
      { 
        i++; 
        scatter = atoi(argv[i]); 
      } 
      else if (strcmp(argv[i],"-grid") == 0) 
      {                               
        i++; 
        GRID_PRECISION = atoi(argv[i]); 
      } 
      else if (strcmp(argv[i],"-ispa") == 0 || strcmp(argv[i],"-spa") == 0) 
      { 
        ispa = true; 
      } 
      else if (strcmp(argv[i],"-ispb") == 0 || strcmp(argv[i],"-spb") == 0) 
      { 
        ispb = true; 
      } 
      else if (strcmp(argv[i],"-inode") == 0 || strcmp(argv[i],"-node") == 0) 
      { 
        inode = true; 
      } 
      else if (strcmp(argv[i],"-isma") == 0 || strcmp(argv[i],"-sma") == 0) 
      { 
        isma = true; 
      } 
      else if (strcmp(argv[i],"-ismb") == 0 || strcmp(argv[i],"-smb") == 0) 
      { 
        ismb = true; 
      } 
      else if (strcmp(argv[i],"-ismc") == 0 || strcmp(argv[i],"-smc") == 0) 
      { 
        ismc = true; 
      } 
      else if (strcmp(argv[i],"-clarkson") == 0 || strcmp(argv[i],"-c") == 0)
      { 
				// this flag requires -hi (.sm file) and -i (.sp file)
				g_clarkson = true;
      } 
      else if (strcmp(argv[i],"-clarksonSmOnly") == 0 || strcmp(argv[i],"-csm") == 0)
      { 
				// this flag requires -hi (.sm file) and -i (.sp file)
				g_clarkson = true;
				g_sm_only = true;
      } 
      else if (strcmp(argv[i],"-i") == 0) 
      { 
        i++; 
        file_name = argv[i]; 
      } 
      else if (strcmp(argv[i],"-hi") == 0) 
      { 
        i++; 
        file_name_header = argv[i]; 
      } 
      else if (i == argc-1) 
      { 
        file_name = argv[argc-1]; 
      } 
    } 
  } 
 
  if (isma || ismb) 
  { 
    file_name = 0; 
  }
 
  pq = new PositionQuantizerNew(); 
  grid_hash = new my_grid_hash; 
 
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 
  glutInitWindowSize(WindowW,WindowH); 

  glutInitWindowPosition(20,20); 

  glutCreateWindow("Streaming Point Cloud Viewer"); 
   
  glShadeModel(GL_FLAT); 
   
  InitColors(); 
  InitLight(); 

	if( g_clarkson ) {
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	}
   
  glutDisplayFunc(myDisplay); 
  glutReshapeFunc(myReshape); 
  glutIdleFunc(myIdle); 
   
  glutMouseFunc(myMouseFunc); 
  glutMotionFunc(myMotionFunc); 
  glutKeyboardFunc(myKeyboard); 
     
  // grid sub menu 
  int menuGrid = glutCreateMenu(MyMenuFunc); 
  glutAddMenuEntry("3x3x3 grid", 71); 
  glutAddMenuEntry("4x4x4 grid", 72); 
  glutAddMenuEntry("5x5x5 grid", 73); 
  glutAddMenuEntry("6x6x6 grid", 74); 
  glutAddMenuEntry("7x7x7 grid", 75); 
  glutAddMenuEntry("8x8x8 grid", 76); 
  glutAddMenuEntry("9x9x9 grid", 77); 
  glutAddMenuEntry("10x10x10 grid", 78); 
 
  // steps sub menu 
  int menuSteps = glutCreateMenu(MyMenuFunc); 
  glutAddMenuEntry("in 5 steps", 40); 
  glutAddMenuEntry("in 10 steps", 41); 
  glutAddMenuEntry("in 25 steps", 42); 
  glutAddMenuEntry("in 50 steps", 43); 
  glutAddMenuEntry("in 100 steps", 44); 
  glutAddMenuEntry("in 250 steps", 45); 
  glutAddMenuEntry("in 500 steps", 46); 
  glutAddMenuEntry("in 1000 steps", 47); 
  glutAddMenuEntry("in 10000 steps", 48); 
 
  // main menu 
  glutCreateMenu(MyMenuFunc); 
  glutAddSubMenu("grid ...", menuGrid); 
  glutAddMenuEntry("", 0); 
  glutAddSubMenu("steps ...", menuSteps); 
  glutAddMenuEntry("", 0); 
  glutAddMenuEntry("rotate <SPACE>", 100); 
  glutAddMenuEntry("translate <SPACE>", 101); 
  glutAddMenuEntry("zoom <SPACE>", 102); 
  glutAddMenuEntry("", 0); 
  glutAddMenuEntry("<s>tep", 103); 
  glutAddMenuEntry("<p>lay", 104); 
  glutAddMenuEntry("st<o>p", 105); 
  glutAddMenuEntry("", 0); 
  glutAddMenuEntry("stream <c>oloring", 152); 
  glutAddMenuEntry("render <m>ode", 153); 
  glutAddMenuEntry("", 0); 
  glutAddMenuEntry("<Q>UIT", 109); 
  glutAttachMenu(GLUT_RIGHT_BUTTON); 
 
  glutMainLoop(); 
 
  return 0; 
} 
