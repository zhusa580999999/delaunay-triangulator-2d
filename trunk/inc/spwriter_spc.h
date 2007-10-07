/*
===============================================================================

  FILE:  SPwriter_spc.h
  
  CONTENTS:
  
    Writes a Streaming Mesh in a losslessly compressed binary format using a
    fairly "light-weight" compression scheme for the connectivity.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    14 February 2005 -- created after receiving Henna's sweet Valentine card
  
===============================================================================
*/
#ifndef SPWRITER_SPC_H
#define SPWRITER_SPC_H

#include "spwriter.h"

#include <stdio.h>

class SPwriter_spc : public SPwriter
{
public:

  // spwriter interface function implementations

  void add_comment(const char* comment);

  void set_npoints(int npoints);
  void set_datatype(SPdatatype datatype){};
  void set_boundingbox(const double* bb_min_d, const double* bb_max_d){};
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f);
  void set_boundingbox(const int* bb_min_i, const int* bb_max_i);

  void write_header();

  void write_point(const double* p_pos_d){};
  void write_point(const float* p_pos_f);
  void write_point(const int* p_pos_i);

  void write_finalize(); // finalizes all points written so far
  void write_finalize_cell(int idx){};

  void close();

  // spwriter_spc functions

  bool open(FILE* fd, int nbits=16);

  SPwriter_spc();
  ~SPwriter_spc();

private:
  int nbits;
  int method;
};

#endif
