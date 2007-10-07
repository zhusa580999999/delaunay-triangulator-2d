/*
===============================================================================

  FILE:  SVwriter_nil.h
  
  CONTENTS:
  
    A Streaming Volume Mesh writer that is a dummy. It does not do anything.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2006  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    03 January 2006 -- created after eating the gummybears henna got in eugene  

===============================================================================
*/
#ifndef SVWRITER_NIL_H
#define SVWRITER_NIL_H

#include "svwriter.h"

class SVwriter_nil : public SVwriter
{
public:

  // svwriter interface function implementations

  void add_comment(const char* comment) {};

  void set_nverts(int nverts);
  void set_ncells(int ncells);

  void set_v_pnum_f(int v_pnum_f) {};
  void set_v_pbox_f(const float* v_pmin_f, const float* v_pmax_f) {};
  void set_v_pnum_i(int v_pnum_i) {};
  void set_v_pbox_i(const int* v_pmin_i, const int* v_pmax_i) {};
  void set_c_pnum_f(int c_pnum_f) {};
  void set_c_pbox_f(const float* c_pmin_f, const float* c_pmax_f) {};
  void set_c_pnum_i(int c_pnum_i) {};
  void set_c_pbox_i(const int* c_pmin_i, const int* c_pmax_i) {};

  void write_vertex(const float* v_prop_f);
  void write_vertex(const int* v_prop_i);
  void write_vertex(const float* v_prop_f, const int* v_prop_i);

  void write_tetrahedron(const int* c_idx, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const float* c_prop_f, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const int* c_prop_i, const bool* c_final=0);
  void write_tetrahedron(const int* c_idx, const float* c_prop_f, const int* c_prop_i, const bool* c_final=0);

  void write_finalized(int final_idx) {};

  void close() {};

  // svwriter_sva functions

  SVwriter_nil() {nverts = -1; ncells = -1; v_count = 0; c_count = 0;};
  ~SVwriter_nil() {};
};

void SVwriter_nil::set_nverts(int nverts)
{
  this->nverts = nverts;
}
void SVwriter_nil::set_ncells(int ncells)
{
  this->ncells = ncells;
}
void SVwriter_nil::write_vertex(const float* v_prop_f)
{
  v_count++;
}
void SVwriter_nil::write_vertex(const int* v_prop_i)
{
  v_count++;
}
void SVwriter_nil::write_vertex(const float* v_prop_f, const int* v_prop_i)
{
  v_count++;
}
void SVwriter_nil::write_tetrahedron(const int* c_idx, const bool* c_final)
{
  c_count++;
}
void SVwriter_nil::write_tetrahedron(const int* c_idx, const float* c_prop_f, const bool* c_final)
{
  c_count++;
}
void SVwriter_nil::write_tetrahedron(const int* c_idx, const int* c_prop_i, const bool* c_final)
{
  c_count++;
}
void SVwriter_nil::write_tetrahedron(const int* c_idx, const float* c_prop_f, const int* c_prop_i, const bool* c_final)
{
  c_count++;
}
#endif
