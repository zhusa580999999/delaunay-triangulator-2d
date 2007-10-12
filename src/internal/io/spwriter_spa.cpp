/*
===============================================================================

  FILE:  SPwriter_spa.cpp
  
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
#include "spwriter_spa.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

bool SPwriter_spa::open(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPwriter_spa\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
    }
  }
#endif

  this->file = file;
  ncomments = 0;
  p_count = 0;
  datatype = SP_FLOAT; // default data type
  finalizemethod = SP_NONE; // default finalize method
  return true;
}

void SPwriter_spa::add_comment(const char* comment)
{
  if (comments == 0)
  {
    comments = (char**)malloc(sizeof(char*)*10);
    comments[9] = (char*)-1;
  }
  else if (comments[ncomments] == (char*)-1)
  {
    comments = (char**)realloc(comments,sizeof(char*)*ncomments*2);
    comments[ncomments*2-1] = (char*)-1;
  }
  comments[ncomments] = strdup(comment);
  ncomments++;
}

void SPwriter_spa::set_datatype(SPdatatype datatype)
{
  this->datatype = datatype;
}

void SPwriter_spa::set_finalizemethod(SPfinalizemethod finalizemethod)
{
  this->finalizemethod = finalizemethod;
}

void SPwriter_spa::set_npoints(int npoints)
{
  this->npoints = npoints;
}

void SPwriter_spa::set_boundingbox(const double* bb_min_d, const double* bb_max_d)
{
  if (this->bb_min_d == 0) this->bb_min_d = new double[3];
  if (this->bb_max_d == 0) this->bb_max_d = new double[3];
  VecCopy3dv(this->bb_min_d, bb_min_d);
  VecCopy3dv(this->bb_max_d, bb_max_d);
}

void SPwriter_spa::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];
  VecCopy3fv(this->bb_min_f, bb_min_f);
  VecCopy3fv(this->bb_max_f, bb_max_f);
}

void SPwriter_spa::set_boundingbox(const int* bb_min_i, const int* bb_max_i)
{
  if (this->bb_min_i == 0) this->bb_min_i = new int[3];
  if (this->bb_max_i == 0) this->bb_max_i = new int[3];
  VecCopy3iv(this->bb_min_i, bb_min_i);
  VecCopy3iv(this->bb_max_i, bb_max_i);
}

void SPwriter_spa::write_header()
{
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fprintf(file, "# %s\012",comments[i]);
    }
  }
  if (npoints != -1) fprintf(file, "# npoints %d\012",npoints);
  if (bb_min_d) fprintf(file, "# bb_min_d %lg %lg %lg\012",bb_min_d[0],bb_min_d[1],bb_min_d[2]);
  if (bb_max_d) fprintf(file, "# bb_max_d %lg %lg %lg\012",bb_max_d[0],bb_max_d[1],bb_max_d[2]);
  if (bb_min_f) fprintf(file, "# bb_min_f %g %g %g\012",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
  if (bb_max_f) fprintf(file, "# bb_max_f %g %g %g\012",bb_max_f[0],bb_max_f[1],bb_max_f[2]);
  if (bb_min_i) fprintf(file, "# bb_min_i %d %d %d\012",bb_min_i[0],bb_min_i[1],bb_min_i[2]);
  if (bb_max_i) fprintf(file, "# bb_max_i %d %d %d\012",bb_max_i[0],bb_max_i[1],bb_max_i[2]);
  if (bb_max_i) fprintf(file, "# bb_max_i %d %d %d\012",bb_max_i[0],bb_max_i[1],bb_max_i[2]);
  switch (datatype)
  {
  case SP_FLOAT:
    fprintf(file, "# datatype SP_FLOAT\012");
    break;
  case SP_DOUBLE:
    fprintf(file, "# datatype SP_DOUBLE\012");
    break;
  case SP_INT:
    fprintf(file, "# datatype SP_INT\012");
    break;
  }
  switch (finalizemethod)
  {
  case SP_QUAD_TREE:
    fprintf(file, "# finalizemethod SP_QUAD_TREE\012");
    break;
  case SP_OCT_TREE:
    fprintf(file, "# finalizemethod SP_OCT_TREE\012");
    break;
  case SP_CLARKSON_2D:
    fprintf(file, "# finalizemethod SP_CLARKSON_2D\012");
    break;
  case SP_CLARKSON_3D:
    fprintf(file, "# finalizemethod SP_CLARKSON_3D\012");
    break;
  }
}

void SPwriter_spa::write_point(const double* p_pos_d)
{
  fprintf(file, "v %lg %lg %lg\012",p_pos_d[0],p_pos_d[1],p_pos_d[2]);
  p_count++;
}

void SPwriter_spa::write_point(const float* p_pos_f)
{
  fprintf(file, "v %g %g %g\012",p_pos_f[0],p_pos_f[1],p_pos_f[2]);
  p_count++;
}

void SPwriter_spa::write_point(const int* p_pos_i)
{
  fprintf(file, "v %d %d %d\012",p_pos_i[0],p_pos_i[1],p_pos_i[2]);
  p_count++;
}

void SPwriter_spa::write_finalize_cell(int idx)
{
  fprintf(file, "x cell %d\012",idx);
}

void SPwriter_spa::close()
{
  file = 0;

  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
    comments = 0;
    ncomments = -1;
  }

  p_count = -1;

  datatype = SP_VOID;

  if (bb_min_d) {delete [] bb_min_d; bb_min_d = 0;}
  if (bb_max_d) {delete [] bb_max_d; bb_max_d = 0;}
  if (bb_min_f) {delete [] bb_min_f; bb_min_f = 0;}
  if (bb_max_f) {delete [] bb_max_f; bb_max_f = 0;}
  if (bb_min_i) {delete [] bb_min_i; bb_min_i = 0;}
  if (bb_max_i) {delete [] bb_max_i; bb_max_i = 0;}
}

SPwriter_spa::SPwriter_spa()
{
  // init of SPwriter interface
  ncomments = -1;
  comments = 0;

  datatype = SP_VOID;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPwriter_spa interface
  file = 0;
}

SPwriter_spa::~SPwriter_spa()
{
  // clean-up for SPwriter interface
  if (p_count != -1)
  {
    close(); // user must have forgotten to close the mesh
  }
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
