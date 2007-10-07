/*
===============================================================================

  FILE:  SPwriter_raw.h
  
  CONTENTS:
  
    Writes raw single-precision point data in an efficient binary format.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2004  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    19 January 2006 -- created for fast read/write of the other data sets
  
===============================================================================
*/
#ifndef SPWRITER_RAW_H
#define SPWRITER_RAW_H

#include "spwriter.h"

#include <stdio.h>

class SPwriter_raw : public SPwriter
{
public:

  // spwriter interface function implementations

  void add_comment(const char* comment);

  void set_npoints(int npoints);
  void set_datatype(SPdatatype datatype);
  void set_boundingbox(const double* bb_min_d, const double* bb_max_d);
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f);
  void set_boundingbox(const int* bb_min_i, const int* bb_max_i);
  void set_finalizemethod(SPfinalizemethod finalizemethod);

  void write_header();

  void write_point(const double* p_pos_d);
  void write_point(const float* v_pos_f);
  void write_point(const int* v_pos_i);

  void write_finalize_cell(int idx);

  void close();

  // spwriter_raw functions

  bool open(FILE* file);

  SPwriter_raw();
  ~SPwriter_raw();

private:
  FILE* file;
  bool final_warn;
};

#endif
