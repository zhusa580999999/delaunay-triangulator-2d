/*
===============================================================================

  FILE:  SSmanager3D.h
  
  CONTENTS:
  
    Streaming Space Manager for three-dimensional spatial streams.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    16 January 2006 -- linked list managing of spheres moved here
    14 January 2006 -- Leo's heroic hacking on correcting the errors offsets  
    10 January 2006 -- SSsphere and floating-point error treatment now included
    05 January 2006 -- use of grid in case brute force search fails
    04 January 2006 -- parental cell_overlaps_sphere test if (radius > cell size)
    03 January 2006 -- finalization check for finalized leaf cells and parents
    01 January 2006 -- faster sphere intersection test with num_children bitmask
    05 October 2005 -- store box of most recent finalized cell for quick accept
    04 October 2005 -- return index on intersecting cell for early reject
    12 September 2005 -- created after brit made dinner but left lots of dishes
  
===============================================================================
*/
#ifndef SSCONTAINER3D_H
#define SSCONTAINER3D_H

class SScell3D
{
public:
  SScell3D* buffer_next;     // used for efficient memory management
  SScell3D* parent;
  int idx;
  int level;
  SScell3D* child[8];
  int num_children;
  float r_min[3];
  float r_mid[3];
  float r_max[3];
  int data;
};

class SSsphere
{  
public:
	int next;
	int prev;
  SScell3D* cell;
  float cen[3];
  float cen_err;
  float rad;  // includes the rad error
  float rad2; // includes the rad2 error
  void initialize(const double* v0, const double* v1, const double* v2, const double* v3);
};

class SScontainer3D
{
public:

  float bb_min_f[3];
  float bb_max_f[3];

  // functions for setup

  void open(const double* bb_min_f, const double* bb_max_f);
  void open(const float* bb_min_d, const float* bb_max_d);
  void close();

  // functions for finalizing space

  int finalize_cell(int cell_idx);

  // function for getting the cell index of a given position at a certain level
  
  int get_idx(const float* pos_f, int level);
  int get_idx(const double* pos_d, int level);

  // get the middle point of the cell containing idx

  void get_mid(float* mid_f, const float* pos_f, int idx);
  void get_mid(float* mid_f, const double* pos_f, int idx);

  // functions for getting the data of an *existing* octree leaf for a given position

  int get_leaf_data(const float* pos_f);

  // functions for testing whether spatial elements are finalized

  bool is_point_finalized(const double * p_pos) const;
  bool is_point_finalized(const float* p_pos) const;

  bool was_sphere_just_finalized(const SSsphere* sphere) const;
  bool is_sphere_finalized(SSsphere* sphere, int s_idx) const; 
  bool is_sphere_finalized_parent(SSsphere* sphere, int s_idx) const; 
  bool is_sphere_finalized_leaf(SSsphere* sphere, int s_idx) const; 

  bool is_halfplane_finalized(const float* p0_pos, const float* p1_pos, const float* p2_pos) const;
  bool is_plane_finalized(const float* p_nor, const float* p_pnt);

  // allow user to store and retrieve data (-> tetrahedra) with cells

  int data;
  int idx;

  bool write_data_to_cell(int idx, int data);
  bool get_data_from_cell(int idx);

  int prepareParent();
  int prepareLeaf();

  // iterate over unfinalized space (for visualization)

  const float* r_min_f;
  const float* r_max_f;

  void iterateInit();
  bool iterateNext();
  bool noiterate(int cell_idx);

  // constructor and destructor

  SScontainer3D();
  ~SScontainer3D();

private:

  SScell3D* root;
  int level_offset[20];

  SScell3D* create_parent_and_siblings(int idx, bool finalize);
  SScell3D* finalized;
  SScell3D* finalized_parent;
  SScell3D* finalized_leaf;
};

#endif
