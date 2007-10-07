/*
===============================================================================

  FILE:  SPreader_tiles.cpp
  
  CONTENTS:
  
    see corresponding header file
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2006  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "spreader_tiles.h"

#include <stdlib.h>

#include "vec3dv.h"
#include "vec3fv.h"

bool SPreader_tiles::open(char* file_name_raw, char* file_name_hdr, int x, int y)
{
  if (file_name_raw == 0 || file_name_hdr == 0)
  {
    fprintf(stderr,"ERROR: file_name_raw or file_name_hdr is zero\n");
    return false;
  }

  FILE* file = fopen(file_name_raw, "rb");

  if (file == 0)
  {
    fprintf(stderr,"ERROR: cannot open %s\n", file_name_raw);
    return false;
  }

  this->file_name = file_name_raw;
  this->file = file;

  FILE* file_hdr = fopen(file_name_hdr, "rb");

  if (file_hdr == 0)
  {
    fprintf(stderr,"ERROR: cannot open %s\n", file_name_hdr);
    return false;
  }

  fread(&npoints,sizeof(int),1,file_hdr);
  fread(&datatype,sizeof(SPdatatype),1,file_hdr);
  if (datatype == SP_DOUBLE)
  {
    if (bb_min_d == 0) bb_min_d = new double[3];
    if (bb_max_d == 0) bb_max_d = new double[3];
    fread(bb_min_d,sizeof(double),3,file_hdr);
    fread(bb_max_d,sizeof(double),3,file_hdr);
  }
  else
  {
    fprintf(stderr,"ERROR: reading SP_DOUBLE header with SPreader_tiles\n");
    exit(1);
  }
  fclose(file_hdr);

  fprintf(stderr,"original points %d times %dx%d tiling makes %d points\n", npoints, x, y, npoints*x*y);
  npoints = npoints*x*y;

  fprintf(stderr,"original bounding box %f %f %f %f\n", bb_min_d[0], bb_min_d[1], bb_max_d[0], bb_max_d[1]);

  bb_range_d[0] = bb_max_d[0] - bb_min_d[0];
  bb_range_d[1] = bb_max_d[1] - bb_min_d[1];

  bb_max_d[0] += (x-1)*bb_range_d[0];
  bb_max_d[1] += (y-1)*bb_range_d[1];

  fprintf(stderr,"enlarged to %f %f %f %f\n", bb_min_d[0], bb_min_d[1], bb_max_d[0], bb_max_d[1]);

  x_max = x;
  y_max = y;

  x_cur = 0;
  y_cur = 0;

  p_count = 0;

  return true;
}

void SPreader_tiles::close()
{
  // close of SPreader interface
  p_count = -1;

  // close of SPreader_tiles
  fclose(file);
  file = 0;
}

SPevent SPreader_tiles::read_event()
{
  if (fread(p_pos_d, sizeof(double), 3, file) == 3)
  {
    p_pos_d[0] += x_cur*bb_range_d[0];
    p_pos_d[1] += y_cur*bb_range_d[1];
    p_count++;
    return SP_POINT;
  }
  else if (x_cur < x_max-1)
  {
    x_cur++;
    // move file pointer back to the beginning
//    fseek(file, 0, SEEK_SET);
    fclose(file);
    file = fopen(file_name, "rb");
    fprintf(stderr,"advancing to tile %d x %d after point %d\n",x_cur, y_cur, p_count);
    return read_event();
  }
  else if (y_cur < y_max-1)
  {
    x_cur = 0;
    y_cur++;
    // move file pointer back to the beginning
//    fseek(file, 0, SEEK_SET);
    fclose(file);
    file = fopen(file_name, "rb");
    fprintf(stderr,"advancing to tile %d x %d after point %d\n",x_cur, y_cur, p_count);
    return read_event();
  }
  return SP_EOF;
}

SPreader_tiles::SPreader_tiles()
{
  // init of SPreader interface
  ncomments = 0;
  comments = 0;

  npoints = -1;
  p_count = -1;

  datatype = SP_DOUBLE;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPreader_tiles
  file = 0;
}

SPreader_tiles::~SPreader_tiles()
{
  // clean-up for SPreader interface
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
  }

  if (bb_min_d) delete [] bb_min_d;
  if (bb_max_d) delete [] bb_max_d;
}
