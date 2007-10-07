/*
===============================================================================

  FILE:  SVwriter.h
  
  CONTENTS:
  
    Streaming Mesh Writer
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    21 February 2005 -- created after booking Henna's ticket to APA conference
  
===============================================================================
*/
#ifndef SVWRITER_H
#define SVWRITER_H

class SVwriter
{
public:

  // mesh variables

  int ncomments;
  char** comments;

  int nverts;
  int ncells;

  int v_count;
  int c_count;

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

  // functions

  virtual void add_comment(const char* comment)=0;

  virtual void set_nverts(int nverts)=0;
  virtual void set_ncells(int ncells)=0;
  virtual void set_v_pnum_f(int v_pnum_f)=0;
  virtual void set_v_pbox_f(const float* v_pmin_f, const float* v_pmax_f)=0;
  virtual void set_v_pnum_i(int v_pnum_i)=0;
  virtual void set_v_pbox_i(const int* v_pmin_i, const int* v_pmax_i)=0;
  virtual void set_c_pnum_f(int c_pnum_f)=0;
  virtual void set_c_pbox_f(const float* c_pmin_f, const float* c_pmax_f)=0;
  virtual void set_c_pnum_i(int c_pnum_i)=0;
  virtual void set_c_pbox_i(const int* c_pmin_i, const int* c_pmax_i)=0;

  virtual void write_vertex(const float* v_data_f)=0;
  virtual void write_vertex(const int* v_data_i)=0;
  virtual void write_vertex(const float* v_data_f, const int* v_data_i)=0;

  virtual void write_tetrahedron(const int* t_idx, const bool* t_final=0)=0;
  virtual void write_tetrahedron(const int* t_idx, const float* c_data_f, const bool* t_final=0)=0;
  virtual void write_tetrahedron(const int* t_idx, const int* c_data_i, const bool* t_final=0)=0;
  virtual void write_tetrahedron(const int* t_idx, const float* c_data_f, const int* c_data_i, const bool* t_final=0)=0;

  virtual void write_finalized(int final_idx)=0;

  virtual void close()=0;

  virtual ~SVwriter(){};
};

#endif
