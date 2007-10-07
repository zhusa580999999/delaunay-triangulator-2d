/*
===============================================================================

  FILE:  SPwriter_spb.h
  
  CONTENTS:
  
    Writes Streaming Points in an efficient binary format.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2004  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 March 2006 -- added set_finalizemethod which is included in header 
    14 June 2005 -- initial version created for Pacific Graphics review purpose
  
===============================================================================
*/
#ifndef SPWRITER_SPB_H
#define SPWRITER_SPB_H

#include "spwriter.h"

#include <stdio.h>

class SPwriter_spb : public SPwriter
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
  void write_point(const float* v_pos_f);
  void write_point(const int* v_pos_i);

  void write_finalize_cell(int idx);

  void close();

  // spwriter_spb functions

  void set_endianness(bool big_endian);

  bool open(FILE* file);

  SPwriter_spb();
  ~SPwriter_spb();

private:
  FILE* file;

  void write_buffer();
  void write_buffer_remaining();

  bool endian_swap;

  int element_size;
  int element_number;
  unsigned int element_descriptor;
  int* element_buffer;
};

#endif
