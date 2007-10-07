/*
===============================================================================

  FILE:  SPreader_raw_d.cpp
  
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
#include "spreader_raw_d.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3dv.h"
#include "vec3fv.h"

bool SPreader_raw_d::open(FILE* file, bool precompute_bounding_box)
{
  if (file == 0)
  {
    fprintf(stderr,"ERROR: file pointer is zero\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to binary (untranslated) mode\n");
    }
  }
#endif

  this->file = file;

  if (precompute_bounding_box)
  {
    // compute the bounding box
    compute_bounding_box();
    // move file pointer back to the beginning
    fseek(file, 0, SEEK_SET);
  }

  p_count = 0;

  return true;
}

void SPreader_raw_d::close()
{
  // close of SPreader interface
  p_count = -1;

  // close of SPreader_raw_d
  file = 0;
}

SPevent SPreader_raw_d::read_event()
{
  if (fread(p_pos_d, sizeof(double), 3, file) == 3)
  {
    p_count++;
    return SP_POINT;
  }
  return SP_EOF;
}

void SPreader_raw_d::compute_bounding_box()
{
  if (bb_min_d == 0) bb_min_d = new double[3]; 
  if (bb_max_d == 0) bb_max_d = new double[3];
  
  if (fread(bb_min_d, sizeof(double), 3, file) == 3)
  {
    p_count = 1;
    VecCopy3dv(bb_max_d, bb_min_d);
    while (fread(p_pos_d, sizeof(double), 3, file) == 3)
    {
      VecUpdateMinMax3dv(bb_min_d, bb_max_d, p_pos_d);
      p_count++;
    }
    npoints = p_count;
  }
  else
  {
    fprintf(stderr,"ERROR: cannot read first point when computing boundingbox\n");
  }
}

bool SPreader_raw_d::load_header(FILE* file)
{
  fread(&npoints,sizeof(int),1,file);
  fread(&datatype,sizeof(SPdatatype),1,file);
  if (datatype == SP_DOUBLE)
  {
    if (bb_min_d == 0) bb_min_d = new double[3];
    if (bb_max_d == 0) bb_max_d = new double[3];
    fread(bb_min_d,sizeof(double),3,file);
    fread(bb_max_d,sizeof(double),3,file);
    if (bb_min_f == 0) bb_min_f = new float[3];
    if (bb_max_f == 0) bb_max_f = new float[3];
    VecCopy3fv(bb_min_f, bb_min_d);
    VecCopy3fv(bb_max_f, bb_max_d);
  }
  else if (datatype == SP_FLOAT)
  {
    if (bb_min_f == 0) bb_min_f = new float[3];
    if (bb_max_f == 0) bb_max_f = new float[3];
    fread(bb_min_f,sizeof(float),3,file);
    fread(bb_max_f,sizeof(float),3,file);
  }
  else if (datatype == SP_INT)
  {
    if (bb_min_i == 0) bb_min_i = new int[3];
    if (bb_max_i == 0) bb_max_i = new int[3];
    fread(bb_min_i,sizeof(int),3,file);
    fread(bb_max_i,sizeof(int),3,file);
  }
  else
  {
    return false;
  }
  return true;
}

SPreader_raw_d::SPreader_raw_d()
{
  // init of SPreader interface
  ncomments = 0;
  comments = 0;

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

  // init of SPreader_raw_d
  file = 0;
}

SPreader_raw_d::~SPreader_raw_d()
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
