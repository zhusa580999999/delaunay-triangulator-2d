/*
===============================================================================

  FILE:  SPwriter_nil.h
  
  CONTENTS:
  
    A Streaming Point writer that is a dummy. It does not really do anything.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    14 January 2006 -- created with fire in living room and xena on table
  
===============================================================================
*/
#ifndef SPWRITER_NIL_H
#define SPWRITER_NIL_H

#include "spwriter.h"

class SPwriter_nil : public SPwriter
{
public:

  // spwriter interface function implementations

  void add_comment(const char* comment) {};

  void set_npoints(int npoints);
  void set_datatype(SPdatatype datatype) {};
  void set_boundingbox(const double* bb_min_d, const double* bb_max_d) {};
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f) {};
  void set_boundingbox(const int* bb_min_i, const int* bb_max_i) {};
  void set_finalizemethod(SPfinalizemethod finalizemethod) {};

  void write_header() {};

  void write_point(const double* p_pos_d);
  void write_point(const float* p_pos_f);
  void write_point(const int* p_pos_i);

  void write_finalize_cell(int idx) {};

  void close() {};

  // spwriter_nil functions

  SPwriter_nil() {npoints = -1; p_count = 0;};
  ~SPwriter_nil() {};
};

void SPwriter_nil::set_npoints(int npoints)
{
  this->npoints = npoints;
}

void SPwriter_nil::write_point(const double* p_pos_d)
{
  p_count++;
}

void SPwriter_nil::write_point(const float* p_pos_f)
{
  p_count++;
}

void SPwriter_nil::write_point(const int* p_pos_i)
{
  p_count++;
}

#endif
