/*
===============================================================================

  FILE:  SPreader.h
  
  CONTENTS:
  
    Streaming Point Reader
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    12 January 2006 -- added double precision fields and SPdatatype enum
    09 February 2005 -- created on my ultimate muskelkater day at MPI
  
===============================================================================
*/
#ifndef SPREADER_H
#define SPREADER_H

#include <stdio.h>

// events 
typedef enum {
  SP_ERROR = -1,
  SP_EOF = 0,
  SP_POINT = 1,
  SP_FINALIZED_CELL = 2,
} SPevent;

#ifndef SPWRITER_H

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

class SPreader
{
public:
  // per point variables
  double p_pos_d[3];
  float p_pos_f[3];
  float p_pos_i[3];

  // space finalization
  int final_idx;

  // per data set variables

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

  virtual SPevent read_event()=0;

  virtual bool open(FILE* file, bool skip_finalize_header = true) {
    fprintf(stderr, "Instance of SPreader does not support open(FILE *)\n");
    return false;
  };

  virtual void close()=0;

  virtual ~SPreader(){};
};

#endif
