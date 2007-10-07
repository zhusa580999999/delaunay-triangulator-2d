/*
===============================================================================

  FILE:  sp_diagram.cpp
  
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
  
    20 April 2005 -- created on submission day for the viz conference
  
===============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GL/glut.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "smreader_sma.h"
#include "smreader_smb.h"
#include "smreader_smc.h"

#include "spconverter.h"
#include "spreader_spa.h"
#include "spreader_spb.h"
#include "spreader_ply.h"
#include "spreader_raw_d.h"

#include "sscontainer2d.h"
#include "sscontainer3d.h"

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
float DistZ=1.44369f;

// VISUALIZATION SETTINGS

float boundingBoxMin[3];
float boundingBoxMax[3];
float boundingBoxScale=1.0f;
float boundingBoxTranslateX = 0.0f;
float boundingBoxTranslateY = 0.0f;
float boundingBoxTranslateZ = 0.0f;

// GLOBAL CONTROL VARIABLES

int WindowW=1024, WindowH=768;
//int WindowW=800, WindowH=600;
//int WindowW=1280, WindowH=960;
int InteractionMode=0;
int AnimationOn=0;
int WorkingOn=0;
int PrintToFileOn=0;

// COLORS

float colours_diffuse[10][4];
float colours_white[4];
float colours_light_blue[4];
float colours_red_white[4];

// DATA STORAGE FOR STREAM VISUALIZER OUTPUT

char* file_name = "points.spa";
char* file_name_header = 0;

bool ispa = false;
bool ispb = false;

bool isma = false;
bool ismb = false;
bool ismc = false;

bool terrain = true;

bool top = false;

SVreader* svreader = 0;
SMreader* smreader = 0;
SPreader* spreader = 0;

SScontainer2D* ss2d = 0;
SScontainer3D* ss3d = 0;

int npoints;

int EXACTLY_N_STEPS = 10;
int EVERY_NTH_STEP = -1;
int NEXT_STEP;

int EXACTLY_N_POINTS = 10000;
int EVERY_NTH_POINT = -1;
int NEXT_POINT;

int DIRTY_MESH=1;
int REPLAY_IT=0;
int REPLAY_COUNT=0;
int GRID_PRECISION=4;
int STREAM_COLORING = 0;
int RENDER_MODE = 0;
int RENDER_BOUNDINGBOX = 0;

unsigned char* Framebuffer = 0;
char* PrintFileName = "frame";
int Time=0;

typedef struct GridCell
{
  int idx;
  int first;
  int last;
  int number;
} GridCell;

typedef struct RenderPoint
{
  float v[3];
  int cell_idx;
} RenderPoint;

typedef __gnu_cxx::hash_map<int, int> my_grid_hash;
static my_grid_hash* grid_hash;

// efficient memory allocation

static GridCell* grid_cell_buffer = 0;
static int grid_cell_buffer_size = 0;
static int grid_cell_buffer_alloc = 0;

static void initGridCellBuffer(int alloc)
{
  if (grid_cell_buffer)
  {
    if (grid_cell_buffer_alloc < alloc)
    {
      grid_cell_buffer_alloc = alloc;
      free(grid_cell_buffer);
      grid_cell_buffer = (GridCell*)malloc(sizeof(GridCell)*grid_cell_buffer_alloc);
    }
  }
  else
  {
    grid_cell_buffer_alloc = alloc;
    grid_cell_buffer = (GridCell*)malloc(sizeof(GridCell)*grid_cell_buffer_alloc);
  }
  grid_cell_buffer_size = 0;
}

static int allocGridCell()
{
  if (grid_cell_buffer_size == grid_cell_buffer_alloc)
  {
    grid_cell_buffer = (GridCell*)realloc(grid_cell_buffer,sizeof(GridCell)*grid_cell_buffer_alloc*2);
    if (!grid_cell_buffer)
    {
      fprintf(stderr,"FATAL ERROR: realloc grid_cell_buffer with %d failed.\n",grid_cell_buffer_alloc*2);
      exit(0);
    }
    grid_cell_buffer_alloc *= 2;
  }
  int index = grid_cell_buffer_size;
  grid_cell_buffer[index].idx = -1;
  grid_cell_buffer[index].first = -1;
  grid_cell_buffer[index].last = -1;
  grid_cell_buffer_size++;
  return index;
}

static void destroyGridCellBuffer()
{
  grid_cell_buffer_size = 0;
  grid_cell_buffer_alloc = 0;
  if (grid_cell_buffer)
  {
    free(grid_cell_buffer);
  }
  grid_cell_buffer = 0;
}

static RenderPoint* render_point_buffer = 0;
static int render_point_buffer_size = 0;
static int render_point_buffer_alloc = 0;

static void initRenderPointBuffer(int alloc)
{
  if (render_point_buffer)
  {
    if (render_point_buffer_alloc < alloc)
    {
      render_point_buffer_alloc = alloc;
      free(render_point_buffer);
      render_point_buffer = (RenderPoint*)malloc(sizeof(RenderPoint)*render_point_buffer_alloc);
    }
  }
  else
  {
    render_point_buffer_alloc = alloc;
    render_point_buffer = (RenderPoint*)malloc(sizeof(RenderPoint)*render_point_buffer_alloc);
  }
  render_point_buffer_size = 0;
}

static int allocRenderPoint()
{
  if (render_point_buffer_size == render_point_buffer_alloc)
  {
    render_point_buffer = (RenderPoint*)realloc(render_point_buffer,sizeof(RenderPoint)*render_point_buffer_alloc*2);
    if (!render_point_buffer)
    {
      fprintf(stderr,"FATAL ERROR: realloc render_point_buffer with %d failed.\n",render_point_buffer_alloc*2);
      exit(0);
    }
    render_point_buffer_alloc *= 2;
  }
  int index = render_point_buffer_size;
  render_point_buffer_size++;
  return index;
}

static void destroyRenderPointBuffer()
{
  render_point_buffer_size = 0;
  render_point_buffer_alloc = 0;
  if (render_point_buffer)
  {
    free(render_point_buffer);
  }
  render_point_buffer = 0;
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
  colours_diffuse[6][0] = 0.0f; colours_diffuse[6][1] = 0.6f; colours_diffuse[6][2] = 0.6f; colours_diffuse[6][3] = 1.0f; // cyan
  colours_diffuse[7][0] = 0.7f; colours_diffuse[7][1] = 0.7f; colours_diffuse[7][2] = 0.7f; colours_diffuse[7][3] = 1.0f; // white
  colours_diffuse[8][0] = 0.2f; colours_diffuse[8][1] = 0.2f; colours_diffuse[8][2] = 0.6f; colours_diffuse[8][3] = 1.0f; // light blue
  colours_diffuse[9][0] = 0.9f; colours_diffuse[9][1] = 0.4f; colours_diffuse[9][2] = 0.7f; colours_diffuse[9][3] = 1.0f; // violett
  
  colours_white[0] = 0.7f; colours_white[1] = 0.7f; colours_white[2] = 0.7f; colours_white[3] = 1.0f; // white
  colours_light_blue[0] = 0.2f; colours_light_blue[1] = 0.2f; colours_light_blue[2] = 0.6f; colours_light_blue[3] = 1.0f; // light blue
  colours_red_white[0] = 0.7f; colours_red_white[1] = 0.7f; colours_red_white[2] = 0.7f; colours_red_white[3] = 1.0f; // white or red
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
  fprintf(stderr,"sp_diagram cloud.spa\n");
  fprintf(stderr,"sp_diagram -ispb < cloud.spb\n");
  fprintf(stderr,"sp_diagram -win 640 480 \n");
  fprintf(stderr,"sp_diagram -win 1600 1200 cloud.obj\n");
  fprintf(stderr,"sp_diagram -h\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"in addition to our formats SPA, SPB, and RAW this tool\n");
  fprintf(stderr,"supports SMx, OBJ, SMF, PLY meshes (optionally gzipped).\n");
}

FILE* file = 0;

void vizBegin()
{
  REPLAY_IT = 0; // just making sure
  DIRTY_MESH = 1;

  if (file_name == 0 && !isma && !ismb && !ismc && !ispa && !ispb)
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

          while (event = spreader_spa->read_event())
          {
            if (event == SP_POINT)
            {
              VecCopy3fv(spreader_spa->bb_min_f, spreader_spa->p_pos_f);
              VecCopy3fv(spreader_spa->bb_max_f, spreader_spa->p_pos_f);
              break;
            }
          }
          while (event = spreader_spa->read_event())
          {
            if (event == SP_POINT)
            {
              VecUpdateMinMax3fv(spreader_spa->bb_min_f, spreader_spa->bb_max_f, spreader_spa->p_pos_f);
            }
          }
          fprintf(stderr, "bb_min_f[0] = %ff; bb_min_f[1] = %ff; bb_min_f[2] = %ff;\n", spreader_spa->bb_min_f[0], spreader_spa->bb_min_f[1], spreader_spa->bb_min_f[2]);
          fprintf(stderr, "bb_max_f[0] = %ff; bb_max_f[1] = %ff; bb_max_f[2] = %ff;\n", spreader_spa->bb_max_f[0], spreader_spa->bb_max_f[1], spreader_spa->bb_max_f[2]);
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

          while (event = spreader_spb->read_event())
          {
            if (event == SP_POINT)
            {
              VecCopy3fv(spreader_spb->bb_min_f, spreader_spb->p_pos_f);
              VecCopy3fv(spreader_spb->bb_max_f, spreader_spb->p_pos_f);
              break;
            }
          }
          while (event = spreader_spb->read_event())
          {
            if (event == SP_POINT)
            {
              VecUpdateMinMax3fv(spreader_spb->bb_min_f, spreader_spb->bb_max_f, spreader_spb->p_pos_f);
            }
          }
          fprintf(stderr, "bb_min_f[0] = %ff; bb_min_f[1] = %ff; bb_min_f[2] = %ff;\n", spreader_spb->bb_min_f[0], spreader_spb->bb_min_f[1], spreader_spb->bb_min_f[2]);
          fprintf(stderr, "bb_max_f[0] = %ff; bb_max_f[1] = %ff; bb_max_f[2] = %ff;\n", spreader_spb->bb_max_f[0], spreader_spb->bb_max_f[1], spreader_spb->bb_max_f[2]);
        }
      
        spreader_spb->close();
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
        spreader_spb->open(file);
      }
      spreader = spreader_spb;
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
    else if (strstr(file_name, ".sma") || strstr(file_name, ".obj") || strstr(file_name, ".smf") || isma)
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

      if (smreader_sma->nverts == -1 || smreader_sma->bb_min_f == 0 || smreader_sma->bb_max_f == 0)
      {
        if (smreader_sma->bb_min_f && smreader_sma->bb_max_f)
        {
          fprintf(stderr,"need additional pass to count verts\n");
          while (smreader_sma->read_element());
        }
        else
        {
          SMevent event;

          if (smreader_sma->nverts == -1)
          {
            fprintf(stderr,"need additional pass to count verts and compute bounding box\n");
          }
          else
          {
            fprintf(stderr,"need additional pass to compute bounding box\n");
          }

          smreader_sma->bb_min_f = new float[3];
          smreader_sma->bb_max_f = new float[3];

          while (event = smreader_sma->read_element())
          {
            if (event == SM_VERTEX)
            {
              VecCopy3fv(smreader_sma->bb_min_f, smreader_sma->v_pos_f);
              VecCopy3fv(smreader_sma->bb_max_f, smreader_sma->v_pos_f);
              break;
            }
          }
          while (event = smreader_sma->read_element())
          {
            if (event == SM_VERTEX)
            {
              VecUpdateMinMax3fv(smreader_sma->bb_min_f, smreader_sma->bb_max_f, smreader_sma->v_pos_f);
            }
          }
          fprintf(stderr, "bb_min_f[0] = %ff; bb_min_f[1] = %ff; bb_min_f[2] = %ff;\n", smreader_sma->bb_min_f[0], smreader_sma->bb_min_f[1], smreader_sma->bb_min_f[2]);
          fprintf(stderr, "bb_max_f[0] = %ff; bb_max_f[1] = %ff; bb_max_f[2] = %ff;\n", smreader_sma->bb_max_f[0], smreader_sma->bb_max_f[1], smreader_sma->bb_max_f[2]);
        }
      
        smreader_sma->close();
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
        smreader_sma->open(file);
      }
      smreader = smreader_sma;
    }
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

      if (smreader_smb->nverts == -1 || smreader_smb->bb_min_f == 0 || smreader_smb->bb_max_f == 0)
      {
        if (smreader_smb->bb_min_f && smreader_smb->bb_max_f)
        {
          fprintf(stderr,"need additional pass to count verts\n");
          while (smreader_smb->read_element());
        }
        else
        {
          SMevent event;

          if (smreader_smb->nverts == -1)
          {
            fprintf(stderr,"need additional pass to count verts and compute bounding box\n");
          }
          else
          {
            fprintf(stderr,"need additional pass to compute bounding box\n");
          }

          smreader_smb->bb_min_f = new float[3];
          smreader_smb->bb_max_f = new float[3];

          while (event = smreader_smb->read_element())
          {
            if (event == SM_VERTEX)
            {
              VecCopy3fv(smreader_smb->bb_min_f, smreader_smb->v_pos_f);
              VecCopy3fv(smreader_smb->bb_max_f, smreader_smb->v_pos_f);
              break;
            }
          }
          while (event = smreader_smb->read_element())
          {
            if (event == SM_VERTEX)
            {
              VecUpdateMinMax3fv(smreader_smb->bb_min_f, smreader_smb->bb_max_f, smreader_smb->v_pos_f);
            }
          }
          fprintf(stderr, "bb_min_f[0] = %ff; bb_min_f[1] = %ff; bb_min_f[2] = %ff;\n", smreader_smb->bb_min_f[0], smreader_smb->bb_min_f[1], smreader_smb->bb_min_f[2]);
          fprintf(stderr, "bb_max_f[0] = %ff; bb_max_f[1] = %ff; bb_max_f[2] = %ff;\n", smreader_smb->bb_max_f[0], smreader_smb->bb_max_f[1], smreader_smb->bb_max_f[2]);
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
    if (ispa)
    {
      SPreader_spa* spreader_spa = new SPreader_spa();
      spreader_spa->open(file);
      spreader = spreader_spa;
      ispa = false;
    }
    else if (ispb)
    {
      SPreader_spb* spreader_spb = new SPreader_spb();
      spreader_spb->open(file);
      spreader = spreader_spb;
      ispb = false;
    }
    else if (isma)
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
    exit(1);
  }

  if (terrain)
  {
    ss2d = new SScontainer2D();
    ss2d->open(spreader->bb_min_f, spreader->bb_max_f);
  }
  else
  {
    ss3d = new SScontainer3D();
    ss3d->open(spreader->bb_min_f, spreader->bb_max_f);
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

  npoints = spreader->npoints;

  fprintf(stderr,"npoints %d\n",npoints);

  EVERY_NTH_STEP = npoints / EXACTLY_N_STEPS;
  if (EVERY_NTH_STEP == 0)
  {
    EVERY_NTH_STEP = 1;
  }
  NEXT_STEP = EVERY_NTH_STEP;

  EVERY_NTH_POINT = npoints / EXACTLY_N_POINTS;
  if (EVERY_NTH_POINT == 0)
  {
    EVERY_NTH_POINT = 1;
  }
  NEXT_POINT = EVERY_NTH_POINT;

  initGridCellBuffer(1024);
  initRenderPointBuffer(2048);

  grid_hash->clear();
}

void vizEnd()
{
  REPLAY_IT = 0; // just making sure
  REPLAY_COUNT = grid_cell_buffer_size;
  DIRTY_MESH = 0;
  if (terrain)
    fprintf(stderr,"occupation %3.1f%% as %d of %dx%d=%d cells are occupied\n", 100.0f*grid_cell_buffer_size/((1<<GRID_PRECISION)*(1<<GRID_PRECISION)), grid_cell_buffer_size, (1<<GRID_PRECISION), (1<<GRID_PRECISION), (1<<GRID_PRECISION)*(1<<GRID_PRECISION));
  else
    fprintf(stderr,"occupation %3.1f%% as %d of %dx%dx%d=%d cells are occupied\n", 100.0f*grid_cell_buffer_size/((1<<GRID_PRECISION)*(1<<GRID_PRECISION)*(1<<GRID_PRECISION)), grid_cell_buffer_size, (1<<GRID_PRECISION), (1<<GRID_PRECISION), (1<<GRID_PRECISION), (1<<GRID_PRECISION)*(1<<GRID_PRECISION)*(1<<GRID_PRECISION));

  fprintf(stderr,"  total # of grid cells: %d\n", grid_cell_buffer_size);
  fprintf(stderr,"  average # vertex/cell: %.1f \n", 1.0f*spreader->p_count/grid_cell_buffer_size);

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

  int cell_idx;
  int grid_idx;
  int point_idx;

  while (event)
  {
    event = spreader->read_event();

    if (event == SP_POINT)
    {

      if (spreader->datatype == SP_DOUBLE)
      {
        VecCopy3fv(spreader->p_pos_f, spreader->p_pos_d);
      }

      if (terrain)
        cell_idx = ss2d->get_idx(spreader->p_pos_f, GRID_PRECISION);
      else
        cell_idx = ss3d->get_idx(spreader->p_pos_f, GRID_PRECISION);

      // check if a cell for this point already exists
      hash_element = grid_hash->find(cell_idx);
      if (hash_element == grid_hash->end())
      {
        grid_idx = allocGridCell();
        grid_hash->insert(my_grid_hash::value_type(cell_idx, grid_idx));
        grid_cell_buffer[grid_idx].first = spreader->p_count;
        grid_cell_buffer[grid_idx].idx = cell_idx;
        grid_cell_buffer[grid_idx].number = 0;
      }
      else
      {
        grid_idx = (*hash_element).second;
      }
      grid_cell_buffer[grid_idx].last = spreader->p_count;
      grid_cell_buffer[grid_idx].number++;
      if (spreader->p_count >= NEXT_POINT)
      {
        point_idx = allocRenderPoint();
        render_point_buffer[point_idx].cell_idx = grid_idx;
        VecCopy3fv(render_point_buffer[point_idx].v, spreader->p_pos_f);
        NEXT_POINT += EVERY_NTH_POINT;
      }
    }
    else if (event == SP_EOF)
    {
      grid_hash->clear();
    }
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
    switch (spreader->read_event())
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
    }
  }
  glEnd();
  glutSwapBuffers();
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
  case '-':
    break;
  case '=':
    break;
  case 'B':
  case 'b':
    RENDER_BOUNDINGBOX = (RENDER_BOUNDINGBOX+1)%3;
    fprintf(stderr,"RENDER_BOUNDINGBOX %d\n",RENDER_BOUNDINGBOX);
    glutPostRedisplay();
    break;
  case 'C':
  case 'c':
    STREAM_COLORING = (STREAM_COLORING+1)%4;
    fprintf(stderr,"STREAM_COLORING %d\n",STREAM_COLORING);
    glutPostRedisplay();
    break;
  case 'M':
  case 'm':
    RENDER_MODE = (RENDER_MODE+1)%3;
    fprintf(stderr,"RENDER_MODE %d\n",RENDER_MODE);
    glutPostRedisplay();
    break;
  case 'O':
  case 'o':
    AnimationOn = 0;
    REPLAY_IT = 0;
    break;
  case 'V':
  case 'v':
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
  case 't':
    top = !top;
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
      NEXT_STEP = render_point_buffer_size / EXACTLY_N_STEPS;
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
      DIRTY_MESH = 1;
  case 'p':
    if (DIRTY_MESH)
    {
      AnimationOn = !AnimationOn;
    }
    else
    {
      if (REPLAY_IT == 0)
      {
        if (REPLAY_COUNT >= grid_cell_buffer_size)
        {
          REPLAY_COUNT = 0;
        }
        NEXT_STEP = grid_cell_buffer_size / EXACTLY_N_STEPS;
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
        vizBegin();
        WorkingOn = vizContinue();
      }
      else
      {
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
      if (REPLAY_COUNT >= grid_cell_buffer_size)
      {
        REPLAY_COUNT = 0;
      }
      NEXT_STEP = grid_cell_buffer_size / EXACTLY_N_STEPS;
      if (NEXT_STEP == 0) NEXT_STEP = 1;
      REPLAY_COUNT += NEXT_STEP;
    }
    glutPostRedisplay();
    break;
  case 'K':
  case 'k':
    printf("Azimuth = %ff;\n",Azimuth);
    printf("Elevation = %ff;\n",Elevation);
    printf("DistX = %ff; DistY = %ff; DistZ = %ff;\n",DistX,DistY,DistZ);
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
  else if (value == 1)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 2)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 3)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 4)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 5)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 6)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 7)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 8)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 9)
  {
    file_name = "";
    DIRTY_MESH = 1;
  }
  else if (value == 10)
  {
    file_name = "";
    DIRTY_MESH = 1;
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


void myDiagramDisplayNew()
{
  int sample[3];

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  glClear(GL_COLOR_BUFFER_BIT);
  
  glViewport(0,0,WindowW,WindowH);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0f,(float)WindowW/WindowH,0.125f,5.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0);
  glRotatef(Elevation,1,0,0);
  glRotatef(Azimuth,0,1,0);

  glDisable(GL_DEPTH_TEST);

//  glLineWidth(1.0f);
  glLineWidth(3.0f);

  glScalef(0.5f, 0.5f, 1.0f);

  glPushMatrix();

  glTranslatef(-1.0, 1.0f, 0.0f);

  glScalef(2.0f/npoints, -2.0f/npoints, 1.0f);

  // draw the spatial span of each sampled point
  if (render_point_buffer_size)
  {
    // z-coordinate
    sample[2] = 0;
    // color
    glColor3fv(colours_diffuse[7]);
    // line width
    glLineWidth(2.0f);
    // start drawing
    glBegin(GL_LINES);
    for (int i = 0; i < render_point_buffer_size; i++)
    {
      // x-coordinate (the spatial appearance of the point)
      sample[0] = i*EVERY_NTH_POINT;

      // y-coordinate (the spatial appearance of the point)
      sample[1] = i*EVERY_NTH_POINT;

      glVertex3iv(sample);

      // x-coordinate (the moment this vertex gets spatially finalized)
      sample[0] = grid_cell_buffer[render_point_buffer[i].cell_idx].last;

      glVertex3iv(sample);
    }
    glEnd();
  }

/*
  // draw the start and end points of already decoded points
  if (grid_cell_buffer_size)
  {
    // z-coordinate
    sample[2] = 0;
    // color
    glColor3fv(colours_diffuse[5]);
    // line width
    glPointSize(2.0f);
    // start drawing
    glBegin(GL_POINTS);
    for (int i = 0; i < grid_cell_buffer_size; i++)
    {
      // x-coordinate (first spatial appearance)
      sample[0] = grid_cell_buffer[i].first;

      // y-coordinate (first spatial appearance)
      sample[1] = grid_cell_buffer[i].first;

      glVertex3iv(sample);

      // x-coordinate (last spatial appearance)
      sample[0] = grid_cell_buffer[i].last;

      glVertex3iv(sample);
    }
    glEnd();
  }
*/

  glPopMatrix();

  glColor3fv(colours_diffuse[0]);
  glBegin(GL_LINE_LOOP);
    glVertex3f(-1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, -1.0f, 0.0f);
    glVertex3f(-1.0f, -1.0f, 0.0f);
  glEnd();


  glutSwapBuffers();
}

void myDiagramDisplay()
{
  int sample[3];

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  glClear(GL_COLOR_BUFFER_BIT);
  
  glViewport(0,0,WindowW,WindowH);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0f,(float)WindowW/WindowH,0.125f,5.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0);
  glRotatef(Elevation,1,0,0);
  glRotatef(Azimuth,0,1,0);

  glDisable(GL_DEPTH_TEST);

//  glLineWidth(1.0f);
  glLineWidth(3.0f);

  glScalef(0.5f, 0.5f, 1.0f);

  glPushMatrix();

  glTranslatef(-1.0, 1.0f, 0.0f);

  glScalef(2.0f/npoints, -2.0f/npoints, 1.0f);

  // draw the spatial span of already occupied cells
  if (grid_cell_buffer_size)
  {
    // z-coordinate
    sample[2] = 0;
    // color
    glColor3fv(colours_diffuse[7]);
    // line width
    glLineWidth(2.0f);
    // start drawing
    glBegin(GL_LINES);
    for (int i = 0; i < grid_cell_buffer_size; i++)
    {
      // x-coordinate (first spatial reference)
      sample[0] = grid_cell_buffer[i].first;

      // y-coordinate (first spatial appearance)
      sample[1] = grid_cell_buffer[i].first;

      glVertex3iv(sample);

      // x-coordinate (last spatial appearance)
      sample[0] = grid_cell_buffer[i].last;

      glVertex3iv(sample);
    }
    glEnd();
  }

  // draw the start and end points of already decoded points
  if (grid_cell_buffer_size)
  {
    // z-coordinate
    sample[2] = 0;
    // color
    glColor3fv(colours_diffuse[5]);
    // line width
    glPointSize(2.0f);
    // start drawing
    glBegin(GL_POINTS);
    for (int i = 0; i < grid_cell_buffer_size; i++)
    {
      // x-coordinate (first spatial appearance)
      sample[0] = grid_cell_buffer[i].first;

      // y-coordinate (first spatial appearance)
      sample[1] = grid_cell_buffer[i].first;

      glVertex3iv(sample);

      // x-coordinate (last spatial appearance)
      sample[0] = grid_cell_buffer[i].last;

      glVertex3iv(sample);
    }
    glEnd();
  }

  glPopMatrix();

  glColor3fv(colours_diffuse[0]);
  glBegin(GL_LINE_LOOP);
    glVertex3f(-1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, 1.0f, 0.0f);
    glVertex3f(1.0f, -1.0f, 0.0f);
    glVertex3f(-1.0f, -1.0f, 0.0f);
  glEnd();


  glutSwapBuffers();
}

void myDisplay()
{
//  myDiagramDisplayNew();
//  return;

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
  glEnable(GL_NORMALIZE);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  int rendered_points;

  if (DIRTY_MESH)
  {
    rendered_points = grid_cell_buffer_size;
  }
  else
  {
    if (REPLAY_COUNT > grid_cell_buffer_size)
    {
      rendered_points = grid_cell_buffer_size;
      REPLAY_IT = 0;
    }
    else
    {
      rendered_points = REPLAY_COUNT;
    }
  }

  // draw cells

  if (grid_cell_buffer_size)
  {
    float color;
    float cell_min[3];
    float cell_max[3];

    if (top)
    {
      cell_min[2] = boundingBoxMax[2];
      cell_max[2] = boundingBoxMax[2];
    }
    else
    {
      cell_min[2] = boundingBoxMin[2];
      cell_max[2] = boundingBoxMin[2];
    }

    if (STREAM_COLORING == 0)
    {
      for (int i = 0; i < rendered_points; i++)
      {
        ss2d->get_min_max(cell_min, cell_max, grid_cell_buffer[i].idx);
        glBegin(GL_TRIANGLE_FAN);
        color = 0.1f+0.8f*grid_cell_buffer[i].first/npoints;
        glColor3f(color,color,color);
        glVertex3f((cell_min[0]+cell_max[0])/2.0f,(cell_min[1]+cell_max[1])/2.0f,cell_min[2]);
        color = 0.1f+0.8f*grid_cell_buffer[i].last/npoints;
        glColor3f(color,color,color);
        glVertex3f(cell_min[0],cell_min[1],cell_min[2]);
        glVertex3f(cell_min[0],cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],cell_min[1],cell_min[2]);
        glVertex3f(cell_min[0],cell_min[1],cell_min[2]);
        glEnd();
      }
    }
    else if (STREAM_COLORING == 1)
    {
      for (int i = 0; i < rendered_points; i++)
      {
        ss2d->get_min_max(cell_min, cell_max, grid_cell_buffer[i].idx);
        glBegin(GL_TRIANGLE_FAN);
        color = 0.1f+0.8f*grid_cell_buffer[i].first/npoints;
        if (grid_cell_buffer[i].first < npoints/3)
        {
          glColor3f(0.1f+0.8f*grid_cell_buffer[i].first/(npoints/3),0.1f,0.1f);
        }
        else if (grid_cell_buffer[i].first < 2*(npoints/3))
        {
          glColor3f(0.9f,0.1f+0.8f*(grid_cell_buffer[i].first-(npoints/3))/(npoints/3),0.1f);
        }
        else
        {
          glColor3f(0.9f, 0.9f, 0.1f+0.8f*(grid_cell_buffer[i].first-2*(npoints/3))/(npoints/3));
        }
        glVertex3f((cell_min[0]+cell_max[0])/2.0f,(cell_min[1]+cell_max[1])/2.0f,cell_min[2]);
        if (grid_cell_buffer[i].last < npoints/3)
        {
          glColor3f(0.1f+0.8f*grid_cell_buffer[i].last/(npoints/3),0.1f,0.1f);
        }
        else if (grid_cell_buffer[i].last < 2*(npoints/3))
        {
          glColor3f(0.9f,0.1f+0.8f*(grid_cell_buffer[i].last-(npoints/3))/(npoints/3),0.1f);
        }
        else
        {
          glColor3f(0.9f, 0.9f, 0.1f+0.8f*(grid_cell_buffer[i].last-2*(npoints/3))/(npoints/3));
        }
        glVertex3f(cell_min[0],cell_min[1],cell_min[2]);
        glVertex3f(cell_min[0],0.5f*(cell_min[1]+cell_max[1]),cell_min[2]);
        glVertex3f(cell_min[0],cell_max[1],cell_min[2]);
        glVertex3f(0.5f*(cell_min[0]+cell_max[0]),cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],0.5f*(cell_min[1]+cell_max[1]),cell_min[2]);
        glVertex3f(cell_max[0],cell_min[1],cell_min[2]);
        glVertex3f(0.5f*(cell_min[0]+cell_max[0]),cell_min[1],cell_min[2]);
        glVertex3f(cell_min[0],cell_min[1],cell_min[2]);
        glEnd();
      }
    }
    else
    {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_diffuse[6]);
      glBegin(GL_QUADS);
      glNormal3f(0,0,-1);
      for (int i = 0; i < rendered_points; i++)
      {
        ss2d->get_min_max(cell_min, cell_max, grid_cell_buffer[i].idx);
        glVertex3f(cell_min[0],cell_min[1],cell_min[2]);
        glVertex3f(cell_min[0],cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],cell_max[1],cell_min[2]);
        glVertex3f(cell_max[0],cell_min[1],cell_min[2]);
      }
      glEnd();
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
    }
  }

  glDisable(GL_NORMALIZE);
  glDisable(GL_DEPTH_TEST);
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

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
  glutSwapBuffers();
}

int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    usage();
  }
  else
  {
    if (strcmp(argv[1],"-h") == 0) 
    {
      usage();
      exit(0);
    }
    else if (strcmp(argv[1],"-help") == 0)
    {
      usage();
      exit(0);
    }
    for (int i = 1; i < argc; i++)
    {
      if (strcmp(argv[i],"-win") == 0)
      {
        i++;
        WindowW = atoi(argv[i]);
        i++;
        WindowH = atoi(argv[i]);
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
      else if (strcmp(argv[i],"-ispa") == 0 || strcmp(argv[i],"-spa") == 0)
      {
        ispa = true;
      }
      else if (strcmp(argv[i],"-ispb") == 0 || strcmp(argv[i],"-spb") == 0)
      {
        ispb = true;
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
      else if (strcmp(argv[i],"-grid") == 0)
      {
        i++;
        GRID_PRECISION = atoi(argv[i]);
      }
      else if (strcmp(argv[i],"-terrain") == 0)
      {
        terrain = true;
      }
      else if (strcmp(argv[i],"-top") == 0)
      {
        top = true;
      }
      else if (i == argc-1)
      {
        file_name = argv[argc-1];
      }
    }
  }

  grid_hash = new my_grid_hash;

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(WindowW,WindowH);
  glutInitWindowPosition(180,100);
  glutCreateWindow("Streaming Point Cloud Viewer");
  
  glShadeModel(GL_SMOOTH);
  
  InitColors();
  InitLight();
  
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

  // samples sub menu
  int menuLines = glutCreateMenu(MyMenuFunc);
  glutAddMenuEntry("10 lines", 81);
  glutAddMenuEntry("15 lines", 82);
  glutAddMenuEntry("20 lines", 83);
  glutAddMenuEntry("25 lines", 84);
  glutAddMenuEntry("40 lines", 85);
  glutAddMenuEntry("50 lines", 86);
  glutAddMenuEntry("60 lines", 87);
  glutAddMenuEntry("75 lines", 88);

  // main menu
  glutCreateMenu(MyMenuFunc);
  glutAddSubMenu("grid ...", menuGrid);
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
  glutAddMenuEntry("render <m>ode", 153);
  glutAddMenuEntry("", 0);
  glutAddMenuEntry("<Q>UIT", 109);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  glutMainLoop();

  return 0;
}
