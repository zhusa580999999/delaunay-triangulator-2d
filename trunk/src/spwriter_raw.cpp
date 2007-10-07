/*
===============================================================================

  FILE:  SPwriter_raw.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2005  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "spwriter_raw.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vec3fv.h"

bool SPwriter_raw::open(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not allowed for SPwriter_raw\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdout to binary (untranslated) mode\n");
    }
  }
#endif

  this->file = file;

  ncomments = 0;
  p_count = 0;
  datatype = SP_FLOAT;

  return true;
}

void SPwriter_raw::add_comment(const char* comment)
{
  fprintf(stderr, "WARNING: add_comment not supported by SPwriter_raw\n");
}

void SPwriter_raw::set_npoints(int npoints)
{
  fprintf(stderr, "WARNING: set_npoints not supported by SPwriter_raw\n");
  this->npoints = npoints;
}

void SPwriter_raw::set_datatype(SPdatatype datatype)
{
  if (datatype != SP_FLOAT)
  {
    fprintf(stderr, "WARNING: only SP_FLOAT is supported by SPwriter_raw\n");
  }
}

void SPwriter_raw::set_boundingbox(const double* bb_min_d, const double* bb_max_d)
{
  fprintf(stderr, "WARNING: set_boundingbox not supported by SPwriter_raw\n");
}

void SPwriter_raw::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  fprintf(stderr, "WARNING: set_boundingbox not supported by SPwriter_raw\n");
}

void SPwriter_raw::set_boundingbox(const int* bb_min_i, const int* bb_max_i)
{
  fprintf(stderr, "WARNING: set_boundingbox not supported by SPwriter_raw\n");
}

void SPwriter_raw::set_finalizemethod(SPfinalizemethod finalizemethod)
{
  fprintf(stderr, "WARNING: set_finalizemethod not supported by SPwriter_raw\n");
}

void SPwriter_raw::write_header()
{
  fprintf(stderr, "WARNING: write_header is meaningless for SPwriter_raw\n");
}

void SPwriter_raw::write_point(const double* p_pos_d)
{
  float p_pos_f[3];
  p_pos_f[0] = (float)p_pos_d[0];
  p_pos_f[1] = (float)p_pos_d[1];
  p_pos_f[2] = (float)p_pos_d[2];
  fwrite(p_pos_f, sizeof(float), 3, file);
  p_count++;
}

void SPwriter_raw::write_point(const float* p_pos_f)
{
  fwrite(p_pos_f, sizeof(float), 3, file);
  p_count++;
}

void SPwriter_raw::write_point(const int* p_pos_i)
{
  float p_pos_f[3];
  p_pos_f[0] = (float)p_pos_i[0];
  p_pos_f[1] = (float)p_pos_i[1];
  p_pos_f[2] = (float)p_pos_i[2];
  fwrite(p_pos_f, sizeof(float), 3, file);
  p_count++;
}

void SPwriter_raw::write_finalize_cell(int idx)
{
  if (!final_warn)
  {
    final_warn = true;
    fprintf(stderr, "WARNING: write_finalize_cell is not supported by SPwriter_raw\n");
  }
}

void SPwriter_raw::close()
{
  file = 0;

  if (npoints != -1)
  {
    if (npoints != p_count)  fprintf(stderr,"WARNING: set npoints to %d but p_count was %d\n",npoints,p_count);
  }
  npoints = p_count;
  p_count = -1;
}

SPwriter_raw::SPwriter_raw()
{
  // init of SPwriter interface
  ncomments = 0;
  comments = 0;

  datatype = SP_FLOAT;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPwriter_raw interface
  file = 0;
  final_warn = false;
}

SPwriter_raw::~SPwriter_raw()
{
  // clean-up for SPwriter interface
  if (p_count != -1)
  {
    close(); // user must have forgotten to close the mesh
  }
}
