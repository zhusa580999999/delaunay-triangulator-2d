/*
===============================================================================

  FILE:  SVreader.h
  
  CONTENTS:
  
    Streaming Volume Mesh Reader
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    20 February 2005 -- created after streaming points turned out too tricky
  
===============================================================================
*/
#ifndef SVREADER_H
#define SVREADER_H

// types 
typedef enum {
  SV_UNDEF = -1,
  SV_TET = 0,
  SV_HEX = 1,
} SVtype;

// events 
typedef enum {
  SV_ERROR = -1,
  SV_EOF = 0,
  SV_VERTEX = 1,
  SV_TETRAHEDRON = 2,
  SV_FINALIZED = 3
} SVevent;

class SVreader
{
public:

  // vertex variables

  float* v_prop_f;
  int*   v_prop_i;

  // cell variables

  int    c_type;
  int*   c_idx;
  bool*  c_final;
  float* c_prop_f;
  int*   c_prop_i;

  // explicit finalization variables

  int final_idx;

  // mesh variables

  SVtype type;

  int ncomments;
  char** comments;

  int v_pnum_f;
  float* v_pmin_f;
  float* v_pmax_f;

  int v_pnum_i;
  int* v_pmin_i;
  int* v_pmax_i;

  int c_pnum_f;
  float* c_pmin_f;
  float* c_pmax_f;

  int c_pnum_i;
  int* c_pmin_i;
  int* c_pmax_i;

  int nverts;
  int ncells;

  int v_count;
  int c_count;

  bool post_order;

  // functions

  virtual SVevent read_element()=0;
  virtual SVevent read_event()=0;

  virtual void close()=0;

  virtual ~SVreader(){};
};

#endif
