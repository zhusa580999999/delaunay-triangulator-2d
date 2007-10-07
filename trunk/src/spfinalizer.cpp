/*
===============================================================================

  FILE:  spfinalize.cpp
  
  CONTENTS:
  
    This spatial finalizer makes two passes over a point cloud. In the
    first pass it creates a "spatial histogram" of the point cloud, which it
    uses in the second pass to add finalization information to the point could.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    16 January 2006 -- merge children with less than minimum_points into parent
    15 January 2006 -- histogram option for cell occupancy
    14 January 2006 -- local (points in cell) & global (point in grid) sprinkle 
    11 January 2006 -- detect & correct wrong npoints in header of input
    15 December 2005 -- early finalization of empty tree cells 
    25 July 2005 -- created after pruning the blackberries at 2611 etna st. 
  
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

//#define OUTPUT_DOUBLE

#ifdef OUTPUT_DOUBLE
#define OUTPUT_DATATYPE double
#else
#define OUTPUT_DATATYPE float
#endif

#include "spreader_spa.h"
#include "spreader_spb.h"
#include "spreader_node.h"
#include "spreader_ply.h"
#include "spreader_raw_d.h"
#include "spreader_raw.h"
#include "spreader_tiles.h"

#include "smreader_sma.h"
#include "smreader_smb.h"
#include "smreader_smc.h"

#include "svreader_sva.h"
#include "svreader_svb.h"
#include "svreader_svc.h"

#include "spconverter.h"

#include "spwriter_spa.h"
#include "spwriter_spb.h"
#include "spwriter_nil.h"

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

#include "sscontainer2d.h"
#include "sscontainer3d.h"

#include <ext/hash_map>

#define GLOBAL_SCATTER

#ifdef _WIN32
extern "C" FILE* fopenGzipped(const char* filename, const char* mode);
extern "C" int gettime_in_msec();
extern "C" int gettime_in_sec();
extern "C" void settime();
#endif

void usage()
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"spfinalize -i terrain.txt -terrain -ospb\n");
  fprintf(stderr,"spfinalize -i cloud.spa -o cloud_finalized.spa\n");
  fprintf(stderr,"spfinalize -i terrain.spb -o terrain_finalized.spb -terrain\n");
  fprintf(stderr,"spfinalize -i mesh.obj.gz -ispa -o cloud_finalized.spa -level 2 -empty 1\n");
  fprintf(stderr,"spfinalize -i mesh.ply -o cloud_finalized.spa -preserve\n");
  fprintf(stderr,"spfinalize -i terrain.txt.gz -o terrain_finalized.spb -sprinkle -level 5 -empty 4\n");
  fprintf(stderr,"spfinalize -h\n");
  exit(0);
}

typedef struct SPFpoint
{
  SPFpoint* buffer_next;    // used for efficient memory management
  OUTPUT_DATATYPE pos[3];
} SPFpoint;

typedef struct SPFcell
{
  SPFcell* buffer_next;     // used for efficient memory management
  SPFpoint* points;
  int num_points;
  SPFcell* parent;
  int num_children;
  int idx;
#ifdef GLOBAL_SCATTER
  SPFpoint* sprinkle_point;
#endif
} SPFcell;

typedef __gnu_cxx::hash_map<int, SPFcell*> my_cell_hash;

// statistics

static int point_buffer_size;
static int cell_buffer_size;

static int point_buffer_maxsize;
static int cell_buffer_maxsize;

// efficient memory management

static int point_buffer_alloc = 16;
static SPFpoint* point_buffer_next = 0;

static int initPointBuffer(int size)
{
  point_buffer_next = (SPFpoint*)malloc(sizeof(SPFpoint)*size);

  if (point_buffer_next == 0)
  {
    fprintf(stderr,"malloc for point buffer failed\n");
    return 0;
  }
  for (int i = 0; i < size; i++)
  {
    point_buffer_next[i].buffer_next = &(point_buffer_next[i+1]);
  }
  point_buffer_next[size-1].buffer_next = 0;
  point_buffer_alloc = size;
  point_buffer_size = 0;
  return 1;
}

static SPFpoint* allocPoint(const float* pos_f)
{
  if (point_buffer_next == 0)
  {
    point_buffer_next = (SPFpoint*)malloc(sizeof(SPFpoint)*point_buffer_alloc);
    if (point_buffer_next == 0)
    {
      fprintf(stderr,"malloc for point buffer failed\n");
      return 0;
    }
    for (int i = 0; i < point_buffer_alloc; i++)
    {
      point_buffer_next[i].buffer_next = &(point_buffer_next[i+1]);
    }
    point_buffer_next[point_buffer_alloc-1].buffer_next = 0;
    point_buffer_alloc = 1.5*point_buffer_alloc;
  }
  // get pointer to next available point
  SPFpoint* point = point_buffer_next;
  point_buffer_next = point->buffer_next;

  // clean point
  point->buffer_next = 0;
#ifdef OUTPUT_DOUBLE
  VecCopy3dv(point->pos, pos_f);
#else
  VecCopy3fv(point->pos, pos_f);
#endif

  point_buffer_size++; if (point_buffer_size > point_buffer_maxsize) point_buffer_maxsize = point_buffer_size;

  return point;
}

static SPFpoint* allocPoint(const double* pos_d)
{
  if (point_buffer_next == 0)
  {
    point_buffer_next = (SPFpoint*)malloc(sizeof(SPFpoint)*point_buffer_alloc);
    if (point_buffer_next == 0)
    {
      fprintf(stderr,"malloc for point buffer failed\n");
      return 0;
    }
    for (int i = 0; i < point_buffer_alloc; i++)
    {
      point_buffer_next[i].buffer_next = &(point_buffer_next[i+1]);
    }
    point_buffer_next[point_buffer_alloc-1].buffer_next = 0;
    point_buffer_alloc = 1.5*point_buffer_alloc;
  }
  // get pointer to next available point
  SPFpoint* point = point_buffer_next;
  point_buffer_next = point->buffer_next;

  // clean point
  point->buffer_next = 0;
#ifdef OUTPUT_DOUBLE
  VecCopy3dv(point->pos, pos_d);
#else
  VecCopy3fv(point->pos, pos_d);
#endif

  point_buffer_size++; if (point_buffer_size > point_buffer_maxsize) point_buffer_maxsize = point_buffer_size;

  return point;
}

static void deallocPoint(SPFpoint* point)
{
  point->buffer_next = point_buffer_next;
  point_buffer_next = point;
  point_buffer_size--;
}

static int cell_buffer_alloc = 16;
static SPFcell* cell_buffer_next = 0;

static int initCellBuffer(int size)
{
  cell_buffer_next = (SPFcell*)malloc(sizeof(SPFcell)*size);

  if (cell_buffer_next == 0)
  {
    fprintf(stderr,"malloc for cell buffer failed\n");
    return 0;
  }
  for (int i = 0; i < size; i++)
  {
    cell_buffer_next[i].buffer_next = &(cell_buffer_next[i+1]);
  }
  cell_buffer_next[size-1].buffer_next = 0;
  cell_buffer_alloc = size;
  cell_buffer_size = 0;
  return 1;
}

static SPFcell* allocCell()
{
  if (cell_buffer_next == 0)
  {
    cell_buffer_next = (SPFcell*)malloc(sizeof(SPFcell)*cell_buffer_alloc);
    if (cell_buffer_next == 0)
    {
      fprintf(stderr,"malloc for cell buffer failed\n");
      return 0;
    }
    for (int i = 0; i < cell_buffer_alloc; i++)
    {
      cell_buffer_next[i].buffer_next = &(cell_buffer_next[i+1]);
    }
    cell_buffer_next[cell_buffer_alloc-1].buffer_next = 0;
    cell_buffer_alloc = 1.5*cell_buffer_alloc;
  }
  // get index of next available cell
  SPFcell* cell = cell_buffer_next;
  cell_buffer_next = cell->buffer_next;
 
  // clean cell
#ifdef GLOBAL_SCATTER
  cell->points = (SPFpoint*)-1;
#else
  cell->points = 0;
#endif

  cell->num_points = 0;
  cell->parent = 0;
  cell->num_children = 0;

  cell_buffer_size++; if (cell_buffer_size > cell_buffer_maxsize) cell_buffer_maxsize = cell_buffer_size;

  return cell;
}

static void deallocCell(SPFcell* cell)
{
  cell->buffer_next = cell_buffer_next;
  cell_buffer_next = cell;
  cell_buffer_size--;
}

#ifdef GLOBAL_SCATTER
static void write_global_sprinkle_points_terrain(SPFcell* cell, SPwriter* spwriter, my_cell_hash* cell_hash)
{
  if (cell->parent)
  {
    int parent_idx = cell->parent->idx;
    // write ancestors
    write_global_sprinkle_points_terrain(cell->parent, spwriter, cell_hash);
    // write siblings (including me)
    for (int s = 1; s < 5; s++)
    {
      // look for siblings in the hash
      my_cell_hash::iterator hash_element = cell_hash->find(4*parent_idx+s);
      // does they exist
      if (hash_element != cell_hash->end())
      {
        cell = (*hash_element).second;
        if (cell->sprinkle_point) // a parent may have taken that sprinkle point
        {
          spwriter->write_point(cell->sprinkle_point->pos);
          deallocPoint(cell->sprinkle_point);
          cell->sprinkle_point = 0;
        }
      }
    }
  }
  else
  {
    assert(cell->idx == 0);
    if (cell->sprinkle_point)
    {
      spwriter->write_point(cell->sprinkle_point->pos);
      deallocPoint(cell->sprinkle_point);
      cell->sprinkle_point = 0;
    }
  }
}
static void write_global_sprinkle_points(SPFcell* cell, SPwriter* spwriter, my_cell_hash* cell_hash)
{
  if (cell->parent)
  {
    int parent_idx = cell->parent->idx;
    // write ancestors
    write_global_sprinkle_points(cell->parent, spwriter, cell_hash);
    // write siblings (including me)
    for (int s = 1; s < 9; s++)
    {
      // look for siblings in the hash
      my_cell_hash::iterator hash_element = cell_hash->find(8*parent_idx+s);
      // does they exist
      if (hash_element != cell_hash->end())
      {
        cell = (*hash_element).second;
        if (cell->sprinkle_point) // a parent may have taken that sprinkle point
        {
          spwriter->write_point(cell->sprinkle_point->pos);
          deallocPoint(cell->sprinkle_point);
          cell->sprinkle_point = 0;
        }
      }
    }
  }
  else
  {
    assert(cell->idx == 0);
    if (cell->sprinkle_point)
    {
      spwriter->write_point(cell->sprinkle_point->pos);
      deallocPoint(cell->sprinkle_point);
      cell->sprinkle_point = 0;
    }
  }
}
#endif

static bool loadHeader(FILE* file, SPreader* spreader)
{
  fread(&(spreader->npoints),sizeof(int),1,file);
  fread(&(spreader->datatype),sizeof(SPdatatype),1,file);
  if (spreader->datatype == SP_DOUBLE)
  {
    if (spreader->bb_min_d == 0) spreader->bb_min_d = new double[3];
    if (spreader->bb_max_d == 0) spreader->bb_max_d = new double[3];
    fread(spreader->bb_min_d,sizeof(double),3,file);
    fread(spreader->bb_max_d,sizeof(double),3,file);
  }
  else if (spreader->datatype == SP_FLOAT)
  {
    if (spreader->bb_min_f == 0) spreader->bb_min_f = new float[3];
    if (spreader->bb_max_f == 0) spreader->bb_max_f = new float[3];
    fread(spreader->bb_min_f,sizeof(float),3,file);
    fread(spreader->bb_max_f,sizeof(float),3,file);
  }
  else if (spreader->datatype == SP_INT)
  {
    if (spreader->bb_min_i == 0) spreader->bb_min_i = new int[3];
    if (spreader->bb_max_i == 0) spreader->bb_max_i = new int[3];
    fread(spreader->bb_min_i,sizeof(int),3,file);
    fread(spreader->bb_max_i,sizeof(int),3,file);
  }
  else
  {
    return false;
  }
  return true;
}

int main(int argc, char *argv[])
{
  int i;
  bool dry = false;
  bool terrain = false;
  bool preserve = false;
  bool sprinkle = true;
  bool histogram = false;
  bool tiles = false;
  int tiles_x = 1;
  int tiles_y = 1;
  int grid_precision = 3;
  int finalize_empty = 0;
  int minimum_points = 0;
  int set_npoints = -1;
  char* file_name_in = 0;
  char* header_file_name_in = 0;
  bool ispa = false;
  bool ispb = false;
  bool ospa = false;
  bool ospb = false;
  bool ospc = false;
  char* file_name_out = 0;
  
  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-dry") == 0)
    {
      dry = true;
    }
    else if (strcmp(argv[i],"-terrain") == 0)
    {
      terrain = true;
    }
    else if (strcmp(argv[i],"-preserve") == 0)
    {
      preserve = true;
    }
    else if (strcmp(argv[i],"-empty") == 0)
    {
      i++;
      finalize_empty = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-min") == 0 || strcmp(argv[i],"-minimum") == 0)
    {
      i++;
      minimum_points = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-nosprinkle") == 0)
    {
      sprinkle = false;
    }
    else if (strcmp(argv[i],"-sprinkle") == 0 || strcmp(argv[i],"-scatter") == 0)
    {
      sprinkle = true;
    }
    else if (strcmp(argv[i],"-grid") == 0 || strcmp(argv[i],"-level") == 0)
    {
      i++;
      grid_precision = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-i") == 0)
    {
      i++;
      file_name_in = argv[i];
    }
    else if (strcmp(argv[i],"-hi") == 0)
    {
      i++;
      header_file_name_in = argv[i];
    }
    else if (strcmp(argv[i],"-o") == 0)
    {
      i++;
      file_name_out = argv[i];
    }
    else if (strcmp(argv[i],"-ispa") == 0)
    {
      ispa = true;
    }
    else if (strcmp(argv[i],"-ispb") == 0)
    {
      ispb = true;
    }
    else if (strcmp(argv[i],"-ospa") == 0)
    {
      ospa = true;
    }
    else if (strcmp(argv[i],"-ospb") == 0)
    {
      ospb = true;
    }
    else if (strcmp(argv[i],"-ospc") == 0)
    {
      ospc = true;
    }
    else if (strcmp(argv[i],"-tiles") == 0)
    {
      tiles = true;
      i++;
      tiles_x = atoi(argv[i]);
      i++;
      tiles_y = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-hist") == 0 || strcmp(argv[i],"-histogram") == 0)
    {
      histogram = true;
    }
    else
    {
      usage();
    }
  }

  if ((file_name_in == 0 && !ispa && ispb) || (file_name_out == 0  && !ospa && !ospb && !ospc && !dry))
  {
    usage();
  }

  if (preserve && sprinkle)
  {
    fprintf(stderr,"WARNING: cannot preserve and sprinkle at the same time. turning off sprinkle.\n");
    sprinkle = false;
  }

  if (finalize_empty > grid_precision)
  {
    fprintf(stderr,"WARNING: lowering empty finalization from %d to %d to match finalization grid level\n", finalize_empty, grid_precision);
    finalize_empty = grid_precision;
  }

  // open file for first pass

  SMreader* smreader = 0;
  SVreader* svreader = 0;
  SPreader* spreader = 0;
  FILE* file_in;

  if (strstr(file_name_in, ".gz"))
  {
#ifdef _WIN32
    if (strstr(file_name_in, ".spa.gz") || strstr(file_name_in, ".node.gz") || strstr(file_name_in, ".sma.gz") || strstr(file_name_in, ".sva.gz") || ispa)
    {
      file_in = fopenGzipped(file_name_in, "r");
    }
    else
    {
      file_in = fopenGzipped(file_name_in, "rb");
    }
#else
    fprintf(stderr,"ERROR: cannot open gzipped file '%s'\n",file_name_in);
    exit(0);
#endif
  }
  else
  {
    if (strstr(file_name_in, ".spa") || strstr(file_name_in, ".node.gz") || strstr(file_name_in, ".sma") || strstr(file_name_in, ".sva") || ispa )
    {
      file_in = fopen(file_name_in, "r");
    }
    else
    {
      file_in = fopen(file_name_in, "rb");
    }
  }

  if (file_in == 0)
  {
    fprintf(stderr,"ERROR: cannot open '%s' for read\n", file_name_in);
    exit(0);
  }
  
  if (strstr(file_name_in, ".spa") || ispa)
  {
    SPreader_spa* spreader_spa = new SPreader_spa();
    spreader_spa->open(file_in);

    if (spreader_spa->bb_min_f == 0 || spreader_spa->bb_max_f == 0)
    {
      SPevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (spreader_spa->bb_min_f == 0) spreader_spa->bb_min_f = new float[3];
      if (spreader_spa->bb_max_f == 0) spreader_spa->bb_max_f = new float[3];

      while ((event = spreader_spa->read_event()) > SM_EOF)
      {
        if (event == SP_POINT)
        {
          VecCopy3fv(spreader_spa->bb_min_f, spreader_spa->p_pos_f);
          VecCopy3fv(spreader_spa->bb_max_f, spreader_spa->p_pos_f);
          break;
        }
      }
      while ((event = spreader_spa->read_event()) > SM_EOF)
      {
        if (event == SP_POINT)
        {
          VecUpdateMinMax3fv(spreader_spa->bb_min_f, spreader_spa->bb_max_f, spreader_spa->p_pos_f);
        }
      }
      // fix for wrongly reported number of points in the header 
      if (spreader_spa->npoints != -1 && spreader_spa->npoints != spreader_spa->p_count) set_npoints = spreader_spa->p_count;
      spreader_spa->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", spreader_spa->bb_min_f[0], spreader_spa->bb_min_f[1], spreader_spa->bb_min_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", spreader_spa->bb_max_f[0], spreader_spa->bb_max_f[1], spreader_spa->bb_max_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "r");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "r");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      spreader_spa->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    spreader = spreader_spa;
  }
  else if (strstr(file_name_in, ".spb") || ispb)
  {
    SPreader_spb* spreader_spb = new SPreader_spb();
    spreader_spb->open(file_in);

    if (spreader_spb->bb_min_f == 0 || spreader_spb->bb_max_f == 0)
    {
      SPevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (spreader_spb->bb_min_f == 0) spreader_spb->bb_min_f = new float[3];
      if (spreader_spb->bb_max_f == 0) spreader_spb->bb_max_f = new float[3];

      while ((event = spreader_spb->read_event()) > SM_EOF)
      {
        if (event == SP_POINT)
        {
          VecCopy3fv(spreader_spb->bb_min_f, spreader_spb->p_pos_f);
          VecCopy3fv(spreader_spb->bb_max_f, spreader_spb->p_pos_f);
          break;
        }
      }
      while ((event = spreader_spb->read_event()) > SM_EOF)
      {
        if (event == SP_POINT)
        {
          VecUpdateMinMax3fv(spreader_spb->bb_min_f, spreader_spb->bb_max_f, spreader_spb->p_pos_f);
        }
      }
      // fix for wrongly reported number of points in the header 
      if (spreader_spb->npoints != -1 && spreader_spb->npoints != spreader_spb->p_count) set_npoints = spreader_spb->p_count;
      spreader_spb->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", spreader_spb->bb_min_f[0], spreader_spb->bb_min_f[1], spreader_spb->bb_min_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", spreader_spb->bb_max_f[0], spreader_spb->bb_max_f[1], spreader_spb->bb_max_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "rb");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "rb");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      spreader_spb->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    spreader = spreader_spb;
  }
  else if (strstr(file_name_in, ".node"))
  {
    SPreader_node* spreader_node = new SPreader_node();
    spreader_node->open(file_in);

    if (spreader_node->bb_min_f == 0 || spreader_node->bb_max_f == 0)
    {
      SPevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (spreader_node->bb_min_f == 0) spreader_node->bb_min_f = new float[3];
      if (spreader_node->bb_max_f == 0) spreader_node->bb_max_f = new float[3];

      while ((event = spreader_node->read_event()) > SP_EOF)
      {
        if (event == SP_POINT)
        {
          VecCopy3fv(spreader_node->bb_min_f, spreader_node->p_pos_f);
          VecCopy3fv(spreader_node->bb_max_f, spreader_node->p_pos_f);
          break;
        }
      }
      while ((event = spreader_node->read_event()) > SP_EOF)
      {
        if (event == SP_POINT)
        {
          VecUpdateMinMax3fv(spreader_node->bb_min_f, spreader_node->bb_max_f, spreader_node->p_pos_f);
        }
      }
      // fix for wrongly reported number of points in the header 
      if (spreader_node->npoints != -1 && spreader_node->npoints != spreader_node->p_count) set_npoints = spreader_node->p_count;
      spreader_node->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", spreader_node->bb_min_f[0], spreader_node->bb_min_f[1], spreader_node->bb_min_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", spreader_node->bb_max_f[0], spreader_node->bb_max_f[1], spreader_node->bb_max_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "r");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "r");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      spreader_node->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    spreader = spreader_node;
  }
  else if (strstr(file_name_in, ".raw_d") && tiles)
  {
    SPreader_tiles* spreader_tiles = new SPreader_tiles();
    if (header_file_name_in)
    {
      spreader_tiles->open(file_name_in,header_file_name_in, tiles_x, tiles_y);
    }
    else
    {
      fprintf(stderr,"no header_file_name_in for SPreader_tiles...\n");
      exit(1);
    }
    spreader = spreader_tiles;
  }
  else if (strstr(file_name_in, ".raw_d"))
  {
    SPreader_raw_d* spreader_raw_d = new SPreader_raw_d();
    spreader_raw_d->open(file_in);
    if (header_file_name_in)
    {
      FILE* file_header = fopen(header_file_name_in, "rb");
      loadHeader(file_header,spreader_raw_d);
      fclose(file_header);
    }
    else
    {
      fprintf(stderr,"no header_file_name_in ...\n");
      exit(1);
    }
    spreader = spreader_raw_d;
  }
  else if (strstr(file_name_in, ".raw"))
  {
    SPreader_raw* spreader_raw = new SPreader_raw();
    spreader_raw->open(file_in);
    if (header_file_name_in)
    {
      FILE* file_header = fopen(header_file_name_in, "rb");
      loadHeader(file_header,spreader_raw);
      fclose(file_header);
    }
    else
    {
      fprintf(stderr,"no header_file_name_in ...\n");
      exit(1);
    }
    spreader = spreader_raw;
  }
  else if (strstr(file_name_in, ".sma"))
  {
    SMreader_sma* smreader_sma = new SMreader_sma();
    smreader_sma->open(file_in);

    if (smreader_sma->bb_min_f == 0 || smreader_sma->bb_max_f == 0)
    {
      SMevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (smreader_sma->bb_min_f == 0) smreader_sma->bb_min_f = new float[3];
      if (smreader_sma->bb_max_f == 0) smreader_sma->bb_max_f = new float[3];

      while ((event = smreader_sma->read_element()) > SM_EOF)
      {
        if (event == SM_VERTEX)
        {
          VecCopy3fv(smreader_sma->bb_min_f, smreader_sma->v_pos_f);
          VecCopy3fv(smreader_sma->bb_max_f, smreader_sma->v_pos_f);
          break;
        }
      }
      while ((event = smreader_sma->read_element()) > SM_EOF)
      {
        if (event == SM_VERTEX)
        {
          VecUpdateMinMax3fv(smreader_sma->bb_min_f, smreader_sma->bb_max_f, smreader_sma->v_pos_f);
        }
      }
      // fix for wrongly reported number of vertices in the header 
      if (smreader_sma->nverts != -1 && smreader_sma->nverts != smreader_sma->v_count) set_npoints = smreader_sma->v_count;
      smreader_sma->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", smreader_sma->bb_min_f[0], smreader_sma->bb_min_f[1], smreader_sma->bb_min_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", smreader_sma->bb_max_f[0], smreader_sma->bb_max_f[1], smreader_sma->bb_max_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "r");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "r");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      smreader_sma->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    smreader = smreader_sma;
  }
  else if (strstr(file_name_in, ".sva"))
  {
    SVreader_sva* svreader_sva = new SVreader_sva();
    svreader_sva->open(file_in);

    if (svreader_sva->v_pmin_f == 0 || svreader_sva->v_pmin_f == 0)
    {
      SVevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (svreader_sva->v_pmin_f == 0) svreader_sva->v_pmin_f = new float[3];
      if (svreader_sva->v_pmax_f == 0) svreader_sva->v_pmax_f = new float[3];

      while ((event = svreader_sva->read_element()) > SV_EOF)
      {
        if (event == SV_VERTEX)
        {
          VecCopy3fv(svreader_sva->v_pmin_f, svreader_sva->v_prop_f);
          VecCopy3fv(svreader_sva->v_pmax_f, svreader_sva->v_prop_f);
          break;
        }
      }
      while ((event = svreader_sva->read_element()) > SV_EOF)
      {
        if (event == SV_VERTEX)
        {
          VecUpdateMinMax3fv(svreader_sva->v_pmin_f, svreader_sva->v_pmax_f, svreader_sva->v_prop_f);
        }
      }
      // fix for wrongly reported number of vertices in the header 
      if (svreader_sva->nverts != -1 && svreader_sva->nverts != svreader_sva->v_count) set_npoints = svreader_sva->v_count;
      svreader_sva->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", svreader_sva->v_pmin_f[0], svreader_sva->v_pmin_f[1], svreader_sva->v_pmin_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", svreader_sva->v_pmax_f[0], svreader_sva->v_pmax_f[1], svreader_sva->v_pmax_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "r");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "r");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      svreader_sva->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    svreader = svreader_sva;
  }
  else if (strstr(file_name_in, ".smb"))
  {
    SMreader_smb* smreader_smb = new SMreader_smb();
    smreader_smb->open(file_in);

    if (smreader_smb->bb_min_f == 0 || smreader_smb->bb_max_f == 0)
    {
      SMevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (smreader_smb->bb_min_f == 0) smreader_smb->bb_min_f = new float[3];
      if (smreader_smb->bb_max_f == 0) smreader_smb->bb_max_f = new float[3];

      while ((event = smreader_smb->read_element()) > SM_EOF)
      {
        if (event == SM_VERTEX)
        {
          VecCopy3fv(smreader_smb->bb_min_f, smreader_smb->v_pos_f);
          VecCopy3fv(smreader_smb->bb_max_f, smreader_smb->v_pos_f);
          break;
        }
      }
      while ((event = smreader_smb->read_element()) > SM_EOF)
      {
        if (event == SM_VERTEX)
        {
          VecUpdateMinMax3fv(smreader_smb->bb_min_f, smreader_smb->bb_max_f, smreader_smb->v_pos_f);
        }
      }
      // fix for wrongly reported number of vertices in the header 
      if (smreader_smb->nverts != -1 && smreader_smb->nverts != smreader_smb->v_count) set_npoints = smreader_smb->v_count;
      smreader_smb->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", smreader_smb->bb_min_f[0], smreader_smb->bb_min_f[1], smreader_smb->bb_min_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", smreader_smb->bb_max_f[0], smreader_smb->bb_max_f[1], smreader_smb->bb_max_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "rb");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "rb");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      smreader_smb->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    smreader = smreader_smb;
  }
  else if (strstr(file_name_in, ".svb"))
  {
    SVreader_svb* svreader_svb = new SVreader_svb();
    svreader_svb->open(file_in);

    if (svreader_svb->v_pmin_f == 0 || svreader_svb->v_pmax_f == 0)
    {
      SVevent event;

      fprintf(stderr,"first pass: compute bounding box\n");
      
#ifdef _WIN32
      settime();
#endif

      if (svreader_svb->v_pmin_f == 0) svreader_svb->v_pmin_f = new float[3];
      if (svreader_svb->v_pmax_f == 0) svreader_svb->v_pmax_f = new float[3];

      while ((event = svreader_svb->read_element()) > SV_EOF)
      {
        if (event == SV_VERTEX)
        {
          VecCopy3fv(svreader_svb->v_pmin_f, svreader_svb->v_prop_f);
          VecCopy3fv(svreader_svb->v_pmax_f, svreader_svb->v_prop_f);
          break;
        }
      }
      while ((event = svreader_svb->read_element()) > SV_EOF)
      {
        if (event == SV_VERTEX)
        {
          VecUpdateMinMax3fv(svreader_svb->v_pmin_f, svreader_svb->v_pmax_f, svreader_svb->v_prop_f);
        }
      }
      // fix for wrongly reported number of vertices in the header 
      if (svreader_svb->nverts != -1 && svreader_svb->nverts != svreader_svb->v_count) set_npoints = svreader_svb->v_count;
      svreader_svb->close();
      fclose(file_in);
      fprintf(stderr, "bb_min %g %g %g\n", svreader_svb->v_pmin_f[0], svreader_svb->v_pmin_f[1], svreader_svb->v_pmin_f[2]);
      fprintf(stderr, "bb_max %g %g %g\n", svreader_svb->v_pmax_f[0], svreader_svb->v_pmax_f[1], svreader_svb->v_pmax_f[2]);
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        file_in = fopenGzipped(file_name_in, "rb");
#endif
      }
      else
      {
        file_in = fopen(file_name_in, "rb");
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' the second time\n", file_name_in);
        exit(0);
      }
      svreader_svb->open(file_in);
#ifdef _WIN32
      fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif
    }
    else
    {
      fprintf(stderr,"have bounding box. first pass not needed.\n");
    }
    svreader = svreader_svb;
  }
  else
  {
    fprintf(stderr,"ERROR: cannot determine which reader to use for '%s'\n", file_name_in);
    exit(0);
  }

  // do we need the converter?

  if (smreader)
  {
    SPconverter* spconverter = new SPconverter();
    spconverter->open(smreader);
    spreader = spconverter;
  }
  else if (svreader)
  {
    SPconverter* spconverter = new SPconverter();
    spconverter->open(svreader);
    spreader = spconverter;
  }

#ifdef OUTPUT_DOUBLE
  if (spreader->datatype == SP_FLOAT)
  {
    fprintf(stderr,"STUPIDITY WARNING: double output wanted from float input? exiting ...\n");
    exit(1);
  }
#else
  if (spreader->datatype == SP_DOUBLE)
  {
    // float output is wanted so turn doubles into floats
    if (spreader->bb_min_f == 0) spreader->bb_min_f = new float[3];
    if (spreader->bb_max_f == 0) spreader->bb_max_f = new float[3];
    VecCopy3fv(spreader->bb_min_f, spreader->bb_min_d);
    VecCopy3fv(spreader->bb_max_f, spreader->bb_max_d);
  }
#endif

  SScontainer2D* ss2d = 0;
  SScontainer3D* ss3d = 0;
  
  if (terrain)
  {
    ss2d = new SScontainer2D();
#ifdef OUTPUT_DOUBLE
    ss2d->open(spreader->bb_min_d,spreader->bb_max_d);
#else
    ss2d->open(spreader->bb_min_f,spreader->bb_max_f);
#endif
  }
  else
  {
    ss3d = new SScontainer3D();
#ifdef OUTPUT_DOUBLE
    ss3d->open(spreader->bb_min_d,spreader->bb_max_d);
#else
    ss3d->open(spreader->bb_min_f,spreader->bb_max_f);
#endif
  }

  // create hash

  my_cell_hash* cell_hash = new my_cell_hash;
  my_cell_hash::iterator hash_element;
  
  SPFcell* cell_list = 0;
  SPFcell* cell;
  int cell_idx;

  int event;

#ifdef _WIN32
  settime();
#endif

  if (terrain)
  {
    fprintf(stderr,"terrain finalization (using 2d quadtree) to level %d with empty %d.\n", grid_precision, finalize_empty);
  }
  else
  {
    fprintf(stderr,"space finalization (using 3d octree) to level %d with empty %d.\n", grid_precision, finalize_empty);
  }

  if (set_npoints != -1)
  {
    fprintf(stderr,"npoints of input file header is corrected from %d to %d\n",spreader->npoints,set_npoints);
    spreader->npoints = set_npoints;
  }

  // start second pass over points
  fprintf(stderr,"second pass: populating finalization counter grid ... \n");

  int sprinkle_count = 0;

  if (spreader->datatype == SP_FLOAT)
  {
    while (event = spreader->read_event())
    {
      switch (event)
      {
      case SP_POINT:
        // get index of the grid cell (a leaf) into which the point falls
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_f, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_f, grid_precision);
        // check if grid cell for the grid position of this point already exists
        hash_element = cell_hash->find(cell_idx);
        if (hash_element == cell_hash->end())
        {
          cell = allocCell();
          // all following points falling into this grid position find their grid cell in the cell_hash
          cell_hash->insert(my_cell_hash::value_type(cell_idx, cell));
          cell->idx = cell_idx;
          cell->buffer_next = cell_list;
          cell_list = cell;
  #ifdef GLOBAL_SCATTER
          if (sprinkle) // always use the first point as the sprinkle point of a leaf
          {
            cell->sprinkle_point = allocPoint(spreader->p_pos_f);
            sprinkle_count++;
          }
  #endif
        }
        else
        {
          cell = (*hash_element).second;
        }
        // increase point counter
        cell->num_points++;
        break;
      default:
        break;
      }
    }
  }
  else
  {
    while (event = spreader->read_event())
    {
      switch (event)
      {
      case SP_POINT:
        // get index of the grid cell (a leaf) into which the point falls
#ifdef OUTPUT_DOUBLE
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_d, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_d, grid_precision);
#else
        VecCopy3fv(spreader->p_pos_f, spreader->p_pos_d);
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_f, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_f, grid_precision);
#endif
        // check if grid cell for the grid position of this point already exists
        hash_element = cell_hash->find(cell_idx);
        if (hash_element == cell_hash->end())
        {
          cell = allocCell();
          // all following points falling into this grid position find their grid cell in the cell_hash
          cell_hash->insert(my_cell_hash::value_type(cell_idx, cell));
          cell->idx = cell_idx;
          cell->buffer_next = cell_list;
          cell_list = cell;
  #ifdef GLOBAL_SCATTER
          if (sprinkle) // always use the first point as the sprinkle point of a leaf
          {
#ifdef OUTPUT_DOUBLE
            cell->sprinkle_point = allocPoint(spreader->p_pos_d);
#else
            cell->sprinkle_point = allocPoint(spreader->p_pos_f);
#endif
            sprinkle_count++;
          }
  #endif
        }
        else
        {
          cell = (*hash_element).second;
        }
        // increase point counter
        cell->num_points++;
        break;
      default:
        break;
      }
    }
  }

  if (sprinkle)
  {
    fprintf(stderr,"second pass: done. created %d leaf grid cells and stored %d sprinkle points\n", cell_buffer_size, sprinkle_count);
  }
  else
  {
    fprintf(stderr,"second pass: done. created %d leaf grid cells\n",cell_buffer_size);
  }

  // close reader

  spreader->close();
  if (file_in && file_name_in) fclose(file_in);

#ifdef _WIN32
  fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif

  // are we supposed to output some info on the occupancy
  
  if (histogram)
  {
    int min_points = spreader->npoints;
    int max_points = 0;
    int logi;
    int logi_max = 0;
    int logi_histogram[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int logi_numpoints[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    cell = cell_list;
    while (cell)
    {
      if (cell->num_points < min_points) min_points = cell->num_points;
      if (cell->num_points > max_points) max_points = cell->num_points;

      i = (cell->num_points)>>1;
      logi = 0;
      while(i)
      {
        i = i >> 1;
        logi++;
      }
      if (logi > logi_max) logi_max = logi;
      logi_histogram[logi]++;
      logi_numpoints[logi] += cell->num_points;
      cell = cell->buffer_next;
    }
    int cum_cell = 0;
    int cum_point = 0;
    for (i = 0; i <= logi_max; i++)
    {
      cum_cell += logi_histogram[i];
      cum_point += logi_numpoints[i];
      if (logi_histogram[i]) fprintf(stderr, "%6d cells (%4.1f%% %5.1f%%) have on average %6d points (%4.1f%% %5.1f%% of total points)\n", logi_histogram[i], ((float)logi_histogram[i]*100.0f/cell_buffer_size), ((float)cum_cell*100.0f/cell_buffer_size), logi_histogram[i] ? logi_numpoints[i]/logi_histogram[i] : 0 ,((float)logi_numpoints[i]*100.0f/spreader->npoints),((float)cum_point*100.0f/spreader->npoints));
    }
    assert(cum_point == spreader->npoints);
    fprintf(stderr, "+---> total of %d cells contain %d points with an average (min/max) of %d (%d/%d)\n",cum_cell,cum_point,cum_point/cum_cell, min_points, max_points);
  }

  // complete finalization grid

  SPFcell* cell_list_other = 0;
  SPFcell* parent = 0;
  SPFcell* sibling = 0;

  if (terrain)
  {
    for (int l = 0; l < grid_precision; l++)
    {
      while (cell = cell_list)
      {
        cell_list = cell->buffer_next;

        if (cell->parent == 0)
        {
          // create parent
          parent = allocCell();
          parent->idx = (cell->idx - 1) / 4;
          // attach parent into tree
          cell->parent = parent;
          parent->num_children = 1;
          // add parent cell into the cell_hash
          cell_hash->insert(my_cell_hash::value_type(parent->idx, parent));
          // add parent to other cell list
          parent->buffer_next = cell_list_other;
          cell_list_other = parent;
  #ifdef GLOBAL_SCATTER
          float mid_f[2]; // for center of parent
          float sprinkle_dist_best;
          SPFcell* sprinkle_best;
          if (sprinkle)
          {
            // get center of parent through position of sprinkle point of child
            ss2d->get_mid(mid_f, cell->sprinkle_point->pos, parent->idx);
            sprinkle_dist_best = (cell->sprinkle_point->pos[0]-mid_f[0])*(cell->sprinkle_point->pos[0]-mid_f[0])+(cell->sprinkle_point->pos[1]-mid_f[1])*(cell->sprinkle_point->pos[1]-mid_f[1]);
            sprinkle_best = cell;
          }
  #endif
          // use i to count the point total of all children cells
          if (cell->num_points > 0)
          {
            i = cell->num_points;
          }
          else
          {
            assert (cell->num_points == 0);
            i = minimum_points; // if the cell has zero points then it's children have enough points already
          }
          // check for siblings
          for (int s = 1; s < 5; s++)
          {
            // compute index of siblings
            cell_idx = 4*parent->idx+s;
            // is it a real sibling (and not the cell itself)
            if (cell_idx != cell->idx)
            {
              // look for it in the hash
              hash_element = cell_hash->find(cell_idx);
              // does it exist
              if (hash_element != cell_hash->end())
              {
                sibling = (*hash_element).second;
                // attach parent into tree
                sibling->parent = parent;
                parent->num_children++;
  #ifdef GLOBAL_SCATTER
                if (sprinkle)
                {
                  float sprinkle_dist = (sibling->sprinkle_point->pos[0]-mid_f[0])*(sibling->sprinkle_point->pos[0]-mid_f[0])+(sibling->sprinkle_point->pos[1]-mid_f[1])*(sibling->sprinkle_point->pos[1]-mid_f[1]);
                  if (sprinkle_dist < sprinkle_dist_best)
                  {
                    sprinkle_best = sibling;
                    sprinkle_dist_best = sprinkle_dist; 
                  }
                }
  #endif
                // use i to count the total points
                if (sibling->num_points > 0)
                {
                  i += sibling->num_points;
                }
                else
                {
                  assert (sibling->num_points == 0);
                  i = minimum_points; // if the sibling has zero points then it's children have enough points
                }
              }
            }
          }
  #ifdef GLOBAL_SCATTER
          if (sprinkle)
          {
            parent->sprinkle_point = sprinkle_best->sprinkle_point;
            sprinkle_best->sprinkle_point = 0; // this point will get sprinkleed by the parent ... not the child
          }
  #endif
          // should we merge the children into one parent node to meet the  requested minimum?
          if (minimum_points && i && (i < minimum_points))
          {
            // the cell should still have the correct parent
            assert(cell->parent == parent);
            // mark the leaf as 'unused' for finalization purposes
            cell->num_points = -cell->num_points;
            // loop over the children again
            for (int s = 1; s < 5; s++)
            {
              // compute index of siblings
              cell_idx = 4*parent->idx + s;
              // is it a real sibling (and not the cell itself)
              if (cell_idx != cell->idx)
              {
                // look for it in the hash
                hash_element = cell_hash->find(cell_idx);
                // does it exist
                if (hash_element != cell_hash->end())
                {
                  sibling = (*hash_element).second;
                  // it should have the correct parent
                  assert(sibling->parent == parent);
                  // mark the leaf as 'unused' for finalization purposes
                  sibling->num_points = -sibling->num_points;
                }
              }
            }
            // the parent should have children 
            assert(parent->num_children > 0);
            // the parent should not yet have points 
            assert(parent->num_points == 0);
            // give the parent the joint point count
            parent->num_points = i;
          }
          // debug assertion
        }
      }
      cell_list =  cell_list_other;
      cell_list_other = 0;
    }
  }
  else
  {
    for (int l = 0; l < grid_precision; l++)
    {
      while (cell = cell_list)
      {
        cell_list = cell->buffer_next;

        if (cell->parent == 0)
        {
          // create parent
          parent = allocCell();
          parent->idx = (cell->idx - 1) / 8;
          // attach parent into tree
          cell->parent = parent;
          parent->num_children = 1;
          // add parent cell into the cell_hash
          cell_hash->insert(my_cell_hash::value_type(parent->idx, parent));
          // add parent to other cell list
          parent->buffer_next = cell_list_other;
          cell_list_other = parent;
  #ifdef GLOBAL_SCATTER
          float mid_f[3]; // for center of parent
          float sprinkle_dist_best;
          SPFcell* sprinkle_best;
          if (sprinkle)
          {
            // get center of parent through position of sprinkle point of child
            ss3d->get_mid(mid_f, cell->sprinkle_point->pos, parent->idx);
            sprinkle_dist_best = VecSquaredDistance3fv(mid_f, cell->sprinkle_point->pos);
            sprinkle_best = cell;
          }
  #endif
          // use i to count the point total of all children cells
          if (cell->num_points > 0)
          {
            i = cell->num_points;
          }
          else
          {
            assert (cell->num_points == 0);
            i = minimum_points; // if the cell has zero points then it's children have enough points already
          }
          // check for siblings
          for (int s = 1; s < 9; s++)
          {
            // compute index of siblings
            cell_idx = 8*parent->idx+s;
            // is it a real sibling (and not the cell itself)
            if (cell_idx != cell->idx)
            {
              // look for it in the hash
              hash_element = cell_hash->find(cell_idx);
              // does it exist
              if (hash_element != cell_hash->end())
              {
                sibling = (*hash_element).second;
                // attach parent into tree
                sibling->parent = parent;
                parent->num_children++;
  #ifdef GLOBAL_SCATTER
                if (sprinkle)
                {
                  float sprinkle_dist = VecSquaredDistance3fv(mid_f, sibling->sprinkle_point->pos);
                  if (sprinkle_dist < sprinkle_dist_best)
                  {
                    sprinkle_best = sibling;
                    sprinkle_dist_best = sprinkle_dist; 
                  }
                }
  #endif
                // use i to count the total points
                if (sibling->num_points > 0)
                {
                  i += sibling->num_points;
                }
                else
                {
                  assert (sibling->num_points == 0);
                  i = minimum_points; // if the sibling has zero points then it's children have enough points
                }
              }
            }
          }
  #ifdef GLOBAL_SCATTER
          if (sprinkle)
          {
            parent->sprinkle_point = sprinkle_best->sprinkle_point;
            sprinkle_best->sprinkle_point = 0; // this point will get sprinkleed by the parent ... not the child
          }
  #endif
          // should we merge the children into one parent node to meet the  requested minimum?
          if (minimum_points && i && (i < minimum_points))
          {
            // the cell should still have the correct parent
            assert(cell->parent == parent);
            // mark the leaf as 'unused' for finalization purposes
            cell->num_points = -cell->num_points;
            // loop over the children again
            for (int s = 1; s < 9; s++)
            {
              // compute index of siblings
              cell_idx = 8*parent->idx + s;
              // is it a real sibling (and not the cell itself)
              if (cell_idx != cell->idx)
              {
                // look for it in the hash
                hash_element = cell_hash->find(cell_idx);
                // does it exist
                if (hash_element != cell_hash->end())
                {
                  sibling = (*hash_element).second;
                  // it should have the correct parent
                  assert(sibling->parent == parent);
                  // mark the leaf as 'unused' for finalization purposes
                  sibling->num_points = -sibling->num_points;
                }
              }
            }
            // the parent should have children 
            assert(parent->num_children > 0);
            // the parent should not yet have points 
            assert(parent->num_points == 0);
            // give the parent the joint point count
            parent->num_points = i;
          }
          // debug assertion
        }
      }
      cell_list =  cell_list_other;
      cell_list_other = 0;
    }
  }

  // make sure the last cell that was created is the root
  assert(parent && parent->idx == 0 || cell && cell->idx == 0);
  
  fprintf(stderr,"created parent grid cells (new total is %d cells)\n",cell_buffer_size);

#ifdef _WIN32
  fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif

  // if we computed an adaptive tree then output the new histogram
  
  if (histogram && minimum_points)
  {
    int min_points = spreader->npoints;
    int max_points = 0;
    int logi;
    int logi_max = 0;
    int logi_histogram[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int logi_numpoints[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // loop over all cells in the hash

    hash_element = cell_hash->begin();
    while (hash_element != cell_hash->end())
    {
      cell = (*hash_element).second;
      if (cell->num_points > 0)
      {
        if (cell->num_points < min_points) min_points = cell->num_points;
        if (cell->num_points > max_points) max_points = cell->num_points;
        i = (cell->num_points)>>1;
        logi = 0;
        while(i)
        {
          i = i >> 1;
          logi++;
        }
        if (logi > logi_max) logi_max = logi;
        logi_histogram[logi]++;
        logi_numpoints[logi] += cell->num_points;
        cell = cell->buffer_next;
      }
      hash_element++;
    }
    int cum_cell = 0;
    int cum_point = 0;
    for (i = 0; i <= logi_max; i++)
    {
      cum_cell += logi_histogram[i];
      cum_point += logi_numpoints[i];
      if (logi_histogram[i]) fprintf(stderr, "%6d cells (%4.1f%% %5.1f%%) have on average %6d points (%4.1f%% %5.1f%% of total points)\n", logi_histogram[i], ((float)logi_histogram[i]*100.0f/cell_buffer_size), ((float)cum_cell*100.0f/cell_buffer_size), logi_histogram[i] ? logi_numpoints[i]/logi_histogram[i] : 0 ,((float)logi_numpoints[i]*100.0f/spreader->npoints),((float)cum_point*100.0f/spreader->npoints));
    }
    assert(cum_point == spreader->npoints);
    fprintf(stderr, "+---> total of %d cells contain %d points with an average (min/max) of %d (%d/%d)\n",cum_cell,cum_point,cum_point/cum_cell, min_points, max_points);
  }

  // open file for second pass

  if (strstr(file_name_in, ".gz"))
  {
#ifdef _WIN32
    if (strstr(file_name_in, ".spa.gz") || strstr(file_name_in, ".sma.gz") || strstr(file_name_in, ".sva.gz") || ispa)
    {
      file_in = fopenGzipped(file_name_in, "r");
    }
    else
    {
      file_in = fopenGzipped(file_name_in, "rb");
    }
#else
    fprintf(stderr,"ERROR: cannot open gzipped file '%s'\n",file_name_in);
    exit(0);
#endif
  }
  else
  {
    if (strstr(file_name_in, ".spa") || strstr(file_name_in, ".sma") || strstr(file_name_in, ".sva") || ispa)
    {
      file_in = fopen(file_name_in, "r");
    }
    else
    {
      file_in = fopen(file_name_in, "rb");
    }
  }

  if (file_in == 0)
  {
    fprintf(stderr,"ERROR: cannot open '%s' for second pass\n", file_name_in);
    exit(0);
  }
  
  if (strstr(file_name_in, ".spa") || ispa)
  {
    SPreader_spa* spreader_spa = (SPreader_spa*)spreader;
    spreader_spa->open(file_in);
  }
  else if (strstr(file_name_in, ".spb") || ispb)
  {
    SPreader_spb* spreader_spb = (SPreader_spb*)spreader;
    spreader_spb->open(file_in);
  }
  else if (strstr(file_name_in, ".node"))
  {
    SPreader_node* spreader_node = (SPreader_node*)spreader;
    spreader_node->open(file_in);
  }
  else if (strstr(file_name_in, ".sma"))
  {
    SMreader_sma* smreader_sma = (SMreader_sma*)smreader;
    smreader_sma->open(file_in);
  }
  else if (strstr(file_name_in, ".smb"))
  {
    SMreader_smb* smreader_smb = (SMreader_smb*)smreader;
    smreader_smb->open(file_in);
  }
  else if (strstr(file_name_in, ".sva"))
  {
    SVreader_sva* svreader_sva = (SVreader_sva*)svreader;
    svreader_sva->open(file_in);
  }
  else if (strstr(file_name_in, ".svb"))
  {
    SVreader_svb* svreader_svb = (SVreader_svb*)svreader;
    svreader_svb->open(file_in);
  }
  else if (strstr(file_name_in, ".raw_d") && tiles)
  {
    SPreader_tiles* spreader_tiles = (SPreader_tiles*)spreader;
    spreader_tiles->open(file_name_in,header_file_name_in, tiles_x, tiles_y);
    spreader = spreader_tiles;
  }
  else if (strstr(file_name_in, ".raw_d"))
  {
    SPreader_raw_d* spreader_raw_d = (SPreader_raw_d*)spreader;
    spreader_raw_d->open(file_in);
  }
  else if (strstr(file_name_in, ".raw"))
  {
    SPreader_raw* spreader_raw = (SPreader_raw*)spreader;
    spreader_raw->open(file_in);
  }
  else
  {
    fprintf(stderr,"ERROR: cannot determine which reader to use for '%s'\n", file_name_in);
    exit(0);
  }

  // open file for output

  SPwriter* spwriter;
  FILE* file_out;

  if (file_name_out || ospa || ospb || ospc)
  {
    if (dry)
    {
      fprintf(stderr,"dry write pass. no output.\n");
      file_out = 0;
      spwriter = new SPwriter_nil();
    }
    else
    {
      if (file_name_out)
      {
        if (strstr(file_name_out, ".spa") || ospa)
        {
          file_out = fopen(file_name_out, "w");
        }
        else if (strstr(file_name_out, ".spb") || strstr(file_name_out, ".spc") || ospb || ospc)
        {
          file_out = fopen(file_name_out, "wb");
        }
        else
        {
          fprintf(stderr,"ERROR: output file name '%s' does not end in .spa or .spb or .spc\n",file_name_out);
          exit(0);
        }
        if (file_out == 0)
        {
          fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_out);
          exit(0);
        }
      }
      else
      {
        file_out = stdout;
      }
    

      if ((file_name_out && strstr(file_name_out, ".spa")) || ospa)
      {
        SPwriter_spa* spwriter_spa = new SPwriter_spa();
        spwriter_spa->open(file_out);
        spwriter = spwriter_spa;
      }
      else if ((file_name_out && strstr(file_name_out, ".spb")) || ospb)
      {
        SPwriter_spb* spwriter_spb = new SPwriter_spb();
        spwriter_spb->open(file_out);
        spwriter = spwriter_spb;
      }
/*
      else if ((file_name_out && strstr(file_name_out, ".spc")) || ospc)
      {
        SPwriter_spc* spwriter_spc = new SPwriter_spc();
        spwriter_spc->open(file_out,bits);
        spwriter = spwriter_spc;
      }
*/
      else
      {
        fprintf(stderr,"ERROR: cannot determine which writer to use for '%s'\n", file_name_out);
        exit(0);
      }
    }
  }
  else
  {
    fprintf(stderr,"no output specified. read pass only.\n");
    spwriter = new SPwriter_nil();
    file_out = 0;
  }

  // do we need the converter again?

  if (smreader)
  {
    SPconverter* spconverter = (SPconverter*)spreader;
    spconverter->open(smreader);
  }

  if (set_npoints != -1)
  {
    fprintf(stderr,"npoints of input file header is corrected from %d to %d\n",spreader->npoints,set_npoints);
    spreader->npoints = set_npoints;
  }

  // start second pass over points
  fprintf(stderr,"third pass: using finalization grid to spatially finalize %d points... \n",spreader->npoints);

  // set the header information
  spwriter->set_npoints(spreader->npoints);
  if (spreader->datatype == SP_FLOAT)
  {
    spwriter->set_datatype(SP_FLOAT);
    spwriter->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
  }
  else
  {
#ifdef OUTPUT_DOUBLE
    fprintf(stderr,"WARNING: output precision is double ... \n");
    spwriter->set_datatype(SP_DOUBLE);
    spwriter->set_boundingbox(spreader->bb_min_d, spreader->bb_max_d);
#else
    fprintf(stderr,"WARNING: forcing output precision to float ... \n");
    spwriter->set_datatype(SP_FLOAT);
    spwriter->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
#endif
  }
  // write the header
  spwriter->write_header();

  // hide whether we do 3d finalization with a -1 finalize call
  if (!terrain) spwriter->write_finalize_cell(-1);

  if (finalize_empty)
  {
    cell_list = parent; // this is the root
    assert(parent->buffer_next == 0);

    fprintf(stderr,"empty finalizations: ");
    while (finalize_empty)
    {
      int count_empty_finalizations = 0;

      parent = cell_list;
      cell_list = 0;
      while (parent)
      {
        if (terrain)
        {
          cell_idx = 4*parent->idx;
          // find and finalize all non-existing children
          for (i = 0; i < 4; i++)
          {
            cell_idx++;
            hash_element = cell_hash->find(cell_idx);
            if (hash_element == cell_hash->end()) // this one gets finalized
            {
              spwriter->write_finalize_cell(cell_idx);
              count_empty_finalizations++;
            }
            else // this one is enlisted to be checked at the next level
            {
              cell = (*hash_element).second;
              cell->buffer_next = cell_list;
              cell_list = cell;
            }
          }
        }
        else
        {
          cell_idx = 8*parent->idx;
          // find and finalize all non-existing children
          for (i = 0; i < 8; i++)
          {
            cell_idx++;
            hash_element = cell_hash->find(cell_idx);
            if (hash_element == cell_hash->end()) // this one gets finalized
            {
              spwriter->write_finalize_cell(cell_idx);
              count_empty_finalizations++;
            }
            else // this one is enlisted to be checked at the next level
            {
              cell = (*hash_element).second;
              cell->buffer_next = cell_list;
              cell_list = cell;
            }
          }
        }
        parent = parent->buffer_next;
      }
      fprintf(stderr,"%d ",count_empty_finalizations);

      finalize_empty--;
    }
    fprintf(stderr,"\n");
  }

  SPFpoint* point;

  if (spreader->datatype == SP_FLOAT)
  {
    while (event = spreader->read_event())
    {
      switch (event)
      {
      case SP_POINT:
        // get index of the grid cell (a leaf) into which the point falls
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_f, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_f, grid_precision);
        // find the grid cell for the grid position of this point
        hash_element = cell_hash->find(cell_idx);
        if (hash_element == cell_hash->end())
        {
          fprintf(stderr,"FATAL ERROR: grid cell of point %d was not in hash\n", spreader->p_count-1);
          exit(0);
        }
        else
        {
          cell = (*hash_element).second;
        }

        // in case we specified a point minimum we may have to move up the tree

        if (minimum_points)
        {
          while (cell->num_points < 0)
          {
            if (cell->num_points == -1)
            {
              // we no longer have children
              assert(cell->num_children == 0);
              // the parent cell must exist
              assert(cell->parent);
              // it should still have children
              assert(cell->parent->num_children > 0);
              // and it loses another child
              cell->parent->num_children--;
              // this cell is not longer needed
              deallocCell(cell);
            }
            else
            {
              cell->num_points++;
            }
            cell = cell->parent;
          }
        }

        if (preserve)
        {
          spwriter->write_point(spreader->p_pos_f);
        }
        else
        {
#ifdef GLOBAL_SCATTER
          if (cell->points == (SPFpoint*)-1) // is this first point for this cell
          {
            // only create it if it was not sprinkleed 
            if (sprinkle)
            {
              cell->points = 0;
            }
            else
            {
              point = allocPoint(spreader->p_pos_f);
              point->buffer_next = 0;
              cell->points = point;
            }
          }
          else
          {
            point = allocPoint(spreader->p_pos_f);
            point->buffer_next = cell->points;
            cell->points = point;
          }
#else // GLOBAL_SCATTER
          point = allocPoint(spreader->p_pos_f);
          point->buffer_next = cell->points;
          cell->points = point;
#endif // GLOBAL_SCATTER
        }

        // finalize if this is the last point of the cell and it has not children
        if (cell->num_points == 1 && cell->num_children == 0)
        {
          assert(cell->parent == 0 || cell->parent->num_points == 0);
#ifdef GLOBAL_SCATTER
          if (sprinkle)
          {
            // write sprinkle points of cell and all ancestors
            if (terrain)
            {
              write_global_sprinkle_points_terrain(cell, spwriter, cell_hash);
            }
            else
            {
              write_global_sprinkle_points(cell, spwriter, cell_hash);
            }
          }
#endif // GLOBAL_SCATTER

//#define REVERSE_POINTS_AGAIN
#define LOCAL_SCATTER_85

#if defined (REVERSE_POINTS_AGAIN)
          // write all the points that were stored with this cell
          SPFpoint* points = 0;
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            point->buffer_next = points;
            points = point;
          }
          while (point = points)
          {
            points = point->buffer_next;
            spwriter->write_point(point->pos);
            deallocPoint(point);
          }
#elif defined (LOCAL_SCATTER_85)
          // organize points into local lists
          SPFpoint* points[3] = {0,0,0};
          int count = 0;
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            if (count % 256 == 0)
            {
              point->buffer_next = points[0];
              points[0] = point;
            }
            else if (count % 32 == 0)
            {
              point->buffer_next = points[1];
              points[1] = point;
            }
            else
            {
              point->buffer_next = points[2];
              points[2] = point;
            }
            count++;
          }
          // write all the points in those lists
          for (int i = 0; i < 3; i++)
          {
            while (point = points[i])
            {
              points[i] = point->buffer_next;
              spwriter->write_point(point->pos);
              deallocPoint(point);
            }
          }
#elif defined (LOCAL_SCATTER_753)
          // organize points into local lists
          SPFpoint* points[4] = {0,0,0,0};
          int count = 0;
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            if (count % 128 == 0)
            {
              point->buffer_next = points[0];
              points[0] = point;
            }
            else if (count % 32 == 0)
            {
              point->buffer_next = points[1];
              points[1] = point;
            }
            else if (count % 8 == 0)
            {
              point->buffer_next = points[2];
              points[2] = point;
            }
            else
            {
              point->buffer_next = points[3];
              points[3] = point;
            }
            count++;
          }
          // write all the points tin those lists
          for (int i = 0; i < 4; i++)
          {
            while (point = points[i])
            {
              points[i] = point->buffer_next;
              spwriter->write_point(point->pos);
              deallocPoint(point);
            }
          }
#else
          // write all the points that were stored with this cell
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            spwriter->write_point(point->pos);
            deallocPoint(point);
          }
#endif

          // finalize cell
          // propagate the finalization event as high as possible
          while (cell->parent && cell->parent->num_children == 1  && cell->parent->num_points == 0)
          {
            parent = cell->parent;
            deallocCell(cell);
            cell = parent;
          }
          if (cell->parent) cell->parent->num_children--;
          spwriter->write_finalize_cell(cell->idx);
          deallocCell(cell);
        }
        else
        {
          cell->num_points--;
        }
        break;
      default:
        break;
      }
    }
  }
  else
  {
    while (event = spreader->read_event())
    {
      switch (event)
      {
      case SP_POINT:
#ifdef OUTPUT_DOUBLE
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_d, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_d, grid_precision);
#else
        VecCopy3fv(spreader->p_pos_f, spreader->p_pos_d);
        if (terrain)
          cell_idx = ss2d->get_idx(spreader->p_pos_f, grid_precision);
        else
          cell_idx = ss3d->get_idx(spreader->p_pos_f, grid_precision);
#endif
        // find the grid cell for the grid position of this point
        hash_element = cell_hash->find(cell_idx);
        if (hash_element == cell_hash->end())
        {
          fprintf(stderr,"FATAL ERROR: grid cell of point %d was not in hash\n", spreader->p_count-1);
          exit(0);
        }
        else
        {
          cell = (*hash_element).second;
        }

        // in case we specified a point minimum we may have to move up the tree

        if (minimum_points)
        {
          while (cell->num_points < 0)
          {
            if (cell->num_points == -1)
            {
              // the parent cell must exist
              assert(cell->parent);
              // it should still have children
              assert(cell->parent->num_children > 0);
              // and it loses another child
              cell->parent->num_children--;
              // this cell is not longer needed
              deallocCell(cell);
            }
            else
            {
              cell->num_points++;
            }
            cell = cell->parent;
          }
        }

        if (preserve)
        {
#ifdef OUTPUT_DOUBLE
          spwriter->write_point(spreader->p_pos_d);
#else
          spwriter->write_point(spreader->p_pos_f);
#endif
        }
        else
        {
#ifdef GLOBAL_SCATTER
          if (cell->points == (SPFpoint*)-1) // is this first point for this cell
          {
            // only create it if it was not sprinkleed 
            if (sprinkle)
            {
              cell->points = 0;
            }
            else
            {
#ifdef OUTPUT_DOUBLE
              point = allocPoint(spreader->p_pos_d);
#else
              point = allocPoint(spreader->p_pos_f);
#endif
              point->buffer_next = 0;
              cell->points = point;
            }
          }
          else
          {
#ifdef OUTPUT_DOUBLE
            point = allocPoint(spreader->p_pos_d);
#else
            point = allocPoint(spreader->p_pos_f);
#endif
            point->buffer_next = cell->points;
            cell->points = point;
          }
#else // GLOBAL_SCATTER
#ifdef OUTPUT_DOUBLE
          point = allocPoint(spreader->p_pos_d);
#else
          point = allocPoint(spreader->p_pos_f);
#endif
          point->buffer_next = cell->points;
          cell->points = point;
#endif // GLOBAL_SCATTER
        }

        // finalize if this is the last point of the cell
        if (cell->num_points == 1)
        {
#ifdef GLOBAL_SCATTER
          if (sprinkle)
          {
            // write sprinkle points of cell and all ancestors
            if (terrain)
            {
              write_global_sprinkle_points_terrain(cell, spwriter, cell_hash);
            }
            else
            {
              write_global_sprinkle_points(cell, spwriter, cell_hash);
            }
          }
#endif // GLOBAL_SCATTER

#if defined (LOCAL_SCATTER_85)
          // organize points into local lists
          SPFpoint* points[3] = {0,0,0};
          int count = 0;
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            if (count % 256 == 0)
            {
              point->buffer_next = points[0];
              points[0] = point;
            }
            else if (count % 32 == 0)
            {
              point->buffer_next = points[1];
              points[1] = point;
            }
            else
            {
              point->buffer_next = points[2];
              points[2] = point;
            }
            count++;
          }
          // write all the points tin those lists
          for (int i = 0; i < 3; i++)
          {
            while (point = points[i])
            {
              points[i] = point->buffer_next;
              spwriter->write_point(point->pos);
              deallocPoint(point);
            }
          }
#elif defined (LOCAL_SCATTER_753)
          // organize points into local lists
          SPFpoint* points[4] = {0,0,0,0};
          int count = 0;
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            if (count % 128 == 0)
            {
              point->buffer_next = points[0];
              points[0] = point;
            }
            else if (count % 32 == 0)
            {
              point->buffer_next = points[1];
              points[1] = point;
            }
            else if (count % 8 == 0)
            {
              point->buffer_next = points[2];
              points[2] = point;
            }
            else
            {
              point->buffer_next = points[3];
              points[3] = point;
            }
            count++;
          }
          // write all the points tin those lists
          for (int i = 0; i < 4; i++)
          {
            while (point = points[i])
            {
              points[i] = point->buffer_next;
              spwriter->write_point(point->pos);
              deallocPoint(point);
            }
          }
#else
          // write all the points that were stored with this cell
          while (point = cell->points)
          {
            cell->points = point->buffer_next;
            spwriter->write_point(point->pos);
            deallocPoint(point);
          }
#endif

          // finalize cell
          // propagate the finalization event as high as possible
          while (cell->parent && cell->parent->num_children == 1 && cell->parent->num_points == 0)
          {
            parent = cell->parent;
            deallocCell(cell);
            cell = parent;
          }
          if (cell->parent) cell->parent->num_children--;
          spwriter->write_finalize_cell(cell->idx);
          deallocCell(cell);
        }
        else
        {
          cell->num_points--;
        }
        break;
      default:
        break;
      }
    }
  }
  
  if (spwriter)
  {
    spwriter->close();
    if (file_out && file_name_out) fclose(file_out);
    delete spwriter;
  }

  if (sprinkle)
  {
    fprintf(stderr,"third pass: done. stored at most %d points (this includes the %d sprinkle points)\n", point_buffer_maxsize, sprinkle_count);
  }
  else
  {
    fprintf(stderr,"third pass: done. stored at most %d points\n",point_buffer_maxsize);
  }

  if (preserve)
  {
    if (point_buffer_maxsize)
    {
      fprintf(stderr,"WARNING: the order was preserved yet some points were stored\n");
    }
    else
    {
      fprintf(stderr,"because the order was preserved no points were stored\n");
    }
  }

#ifdef _WIN32
  fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif

  return 1;
}
