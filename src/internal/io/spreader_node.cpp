/*
===============================================================================

  FILE:  SPreader_node.cpp
  
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
#include "spreader_node.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

bool SPreader_node::open(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPreader_node\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
    }
  }
#endif

  this->file = file;

  line = (char*)malloc(sizeof(char)*256);
  if (fgets(line, sizeof(char) * 256, file) == 0)
  {
    free(line);
    return false;
  }

  // read header information
  int nattrib;
  
  if (sscanf(&(line[0]), "%d %d %d", &npoints, &ncoords, &nattrib) == 3)
  {
    ncoords += nattrib;
  }

  if (fgets(line, sizeof(char) * 256, file) == 0)
  {
    free(line);
    return false;
  }

  p_count = 0;

  return true;
}

void SPreader_node::force_double()
{
  this->datatype = SP_DOUBLE;
}

void SPreader_node::close()
{
  p_count = -1;
}

SPevent SPreader_node::read_event()
{
  while (line)
  {
    int idx;

    if (datatype == SP_FLOAT)
    {
      if (sscanf(&(line[0]), "%d %f %f %f", &idx, &(p_pos_f[0]), &(p_pos_f[1]), &(p_pos_f[2])) == 3)
      {
        p_pos_f[2] = 0;
      }
    }
    else if (datatype == SP_DOUBLE)
    {
      if (sscanf(&(line[0]), "%d %lf %lf %lf", &idx, &(p_pos_d[0]), &(p_pos_d[1]), &(p_pos_d[2])) == 3)
      {
        p_pos_d[2] = 0;
      }
    }
    else
    {
      if (sscanf(&(line[0]), "%d %d %d %d", &idx, &(p_pos_i[0]), &(p_pos_i[1]), &(p_pos_i[2])) == 3)
      {
        p_pos_i[2] = 0;
      }
    }

    p_count++;
    if (fgets(line, sizeof(char) * 256, file) == 0)
    {
      free(line);
      line = 0;
    }
    return SP_POINT;
  }

  if (npoints == -1)
  {
    npoints = p_count;
  }
  else
  {
    if (p_count != npoints)
    {
      fprintf(stderr,"WARNING: wrong point count: p_count (%d) != npoints (%d)\n", p_count, npoints);
    }
  }
  return SP_EOF;
}

SPreader_node::SPreader_node()
{
  // init of SPreader interface
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

  // init of SPreader_node
  file = 0;
  line = 0;
}

SPreader_node::~SPreader_node()
{
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
  if (bb_min_f) delete [] bb_min_f;
  if (bb_max_f) delete [] bb_max_f;
  if (bb_min_i) delete [] bb_min_i;
  if (bb_max_i) delete [] bb_max_i;
}
