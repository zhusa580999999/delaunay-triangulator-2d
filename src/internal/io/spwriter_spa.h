/*
===============================================================================

  FILE:  SPwriter_spa.h
  
  CONTENTS:
  
    Writes Streaming Points in a simple ASCII format. 
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 March 2006 -- added set_finalizemethod which is included in header 
    09 February 2005 -- created on my ultimate muskelkater day at MPI
  
===============================================================================
*/
#ifndef SPWRITER_SPA_H
#define SPWRITER_SPA_H

#include "spwriter.h"

#include <stdio.h>

class SPwriter_spa : public SPwriter
{
public:

  // spwriter interface function implementations

  void add_comment(const char* comment);

  void set_datatype(SPdatatype datatype);
  void set_finalizemethod(SPfinalizemethod finalizemethod);

  void set_npoints(int npoints);
  void set_boundingbox(const double* bb_min_d, const double* bb_max_d);
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f);
  void set_boundingbox(const int* bb_min_i, const int* bb_max_i);

  void write_header();

  void write_point(const double* p_pos_d);
  void write_point(const float* p_pos_f);
  void write_point(const int* p_pos_i);

  void write_finalize_cell(int idx);

  void close();

  // spwriter_spa functions

  bool open(FILE* file);

  SPwriter_spa();
  ~SPwriter_spa();

private:
  FILE* file;
};

#endif
