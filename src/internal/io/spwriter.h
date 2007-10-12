/*
===============================================================================

  FILE:  SPwriter.h
  
  CONTENTS:
  
    Streaming Point Writer
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    09 February 2005 -- created on my ultimate muskelkater day at MPI
  
===============================================================================
*/
#ifndef SPWRITER_H
#define SPWRITER_H

#ifndef SPREADER_H

#include <stdio.h>


// data type of coordinates 
typedef enum {
  SP_VOID = -1,
  SP_FLOAT = 0,
  SP_INT = 1,
  SP_DOUBLE = 2
} SPdatatype;

// method of spatial finalization  
typedef enum {
  SP_NONE = 0,
  SP_QUAD_TREE = 1,
  SP_OCT_TREE = 2,
  SP_CLARKSON_2D = 3,
  SP_CLARKSON_3D = 4
} SPfinalizemethod;

#endif

class SPwriter
{
public:

  // mesh variables

  int ncomments;
  char** comments;

  int npoints;
  int p_count;

  SPdatatype datatype;
  SPfinalizemethod finalizemethod;

  double* bb_min_d;
  double* bb_max_d;
  float* bb_min_f;
  float* bb_max_f;
  int* bb_min_i;
  int* bb_max_i;

  // virtual functions

  virtual void add_comment(const char* comment)=0;

  virtual void set_datatype(SPdatatype datatype)=0;
  virtual void set_finalizemethod(SPfinalizemethod finalizemethod)=0;

  virtual void set_npoints(int npoints)=0;
  virtual void set_boundingbox(const double * bb_min_d, const double* bb_max_d)=0;
  virtual void set_boundingbox(const float* bb_min_f, const float* bb_max_f)=0;
  virtual void set_boundingbox(const int* bb_min_i, const int* bb_max_i)=0;

  virtual void write_header()=0;

  virtual void write_point(const double* p_pos_d)=0;
  virtual void write_point(const float* p_pos_f)=0;
  virtual void write_point(const int* p_pos_i)=0;

  virtual void write_finalize_cell(int idx)=0;

  virtual bool open(FILE* file) {
    fprintf(stderr, "Instance of SPWriter does not support open(FILE *)\n");
    return false;
  };
  virtual void close()=0;

  virtual ~SPwriter(){};
};

#endif
