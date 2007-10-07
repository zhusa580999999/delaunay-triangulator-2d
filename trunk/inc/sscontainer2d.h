/*
===============================================================================

  FILE:  SSmanager2D.h
  
  CONTENTS:
  
    Streaming Space Manager for two-dimensional spatial streams.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    16 January 2006 -- double-precision points find the correct cell
    09 January 2006 -- class SScircle created to make interface cleaner
    08 January 2006 -- integrated correct floating-point error treatment 
    29 August 2005 -- early intersect by storing cell idxs with the triangles
    17 August 2005 -- no real speed-up from srqt()-less squared-radii circles
    15 August 2005 -- more efficient intersection tests using tree structure
    27 July 2005 -- created before cooking my first 2611 etna street dinner
  
===============================================================================
*/
#ifndef SSCONTAINER2D_H
#define SSCONTAINER2D_H

class SScircle
{
public:
  float cen[2];
  float cen_err;
  float rad;  // includes the rad error
  float rad2; // includes the rad2 error
  int cell_idx;
  void initialize(const double* v0, const double* v1, const double* v2);
};

class SScontainer2D
{
public:

  // functions for setup

  void open(const double* bb_min_f, const double* bb_max_f);
  void open(const float* bb_min_d, const float* bb_max_d);
  void close();

  // functions for finalizing space

  void finalize_cell(int idx);

  // functions for getting the index of the cell containing a position at a level
  
  int get_idx(const float* pos_f, int level);
  int get_idx(const double* pos_d, int level);

  // get the middle of the cell containing idx

  void get_mid(float* mid_f, const float* pos_f, int idx);
  void get_mid(float* mid_f, const double* pos_d, int idx);

  // get the min and max of the cell containing idx
  
  void get_min_max(float* min_f, float* max_f, const float* pos_f, int idx);
  void get_min_max(float* min_f, float* max_f, int idx);

  // functions for testing whether spatial elements are finalized

  bool is_finalized(const float* p_pos);
  bool is_finalized(const double* p_pos);
  bool is_finalized(SScircle* circle);
  bool is_still_finalized(SScircle* circle);

  // iterate over unfinalized space

  const float* r_min_f;
  const float* r_max_f;

  void iterateInit();
  bool iterateNext();

  // functions for testing whether spatial elements are finalized

  SScontainer2D();
  ~SScontainer2D();

  // manager variables

  void* root;
  float bb_min_f[2];
  float bb_max_f[2];
  int level_offset[20];
};

#endif
