/*
===============================================================================

  FILE:  SVwriter_off.h
  
  CONTENTS:
  
    Writes a Streaming Volume Mesh in a simple ASCII format.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 June 2005 -- created after a yummy "Kartoffelauflauf" made by Muran

===============================================================================
*/
#ifndef SVWRITER_OFF_H
#define SVWRITER_OFF_H

#include "svwriter.h"

#include <stdio.h>

class SVwriter_off : public SVwriter
{
public:

  // svwriter interface function implementations

  void add_comment(const char* comment);

  void set_nverts(int nverts);
  void set_ncells(int ncells);

  void set_v_pnum_f(int v_pnum_f);
  void set_v_pbox_f(const float* v_pmin_f, const float* v_pmax_f);
  void set_v_pnum_i(int v_pnum_i);
  void set_v_pbox_i(const int* v_pmin_i, const int* v_pmax_i);
  void set_c_pnum_f(int c_pnum_f);
  void set_c_pbox_f(const float* c_pmin_f, const float* c_pmax_f);
  void set_c_pnum_i(int c_pnum_i);
  void set_c_pbox_i(const int* c_pmin_i, const int* c_pmax_i);

  void write_vertex(const float* v_prop_f);
  void write_vertex(const int* v_prop_i);
  void write_vertex(const float* v_prop_f, const int* v_prop_i);

  void write_tetrahedron(const int* c_idx, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const float* c_prop_f, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const int* c_prop_i, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const float* c_prop_f, const int* c_prop_i, const bool* c_final=0);

  void write_finalized(int final_idx);

  void close();

  // svwriter_off functions

  bool open(FILE* file);

  SVwriter_off();
  ~SVwriter_off();

private:
  FILE* file;
  int vertex_buffer_alloc;
  float* vertex_buffer;
  int cell_buffer_alloc;
  int* cell_buffer;
};

#endif
