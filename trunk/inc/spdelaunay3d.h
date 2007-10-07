/*
===============================================================================

  FILE:  spdelaunay3d.h
  
  CONTENTS:
  
    Expects a streaming 3D point cloud as input and produces a Delaunay
    tetrahedralized volume mesh as output.
  
	  Invariant maintained during Delaunay tetrahedralization: A tetrahedron
    is output whenever its circumsphere is spatially finalized. A vertex
    is output whenever the first tetrahedron that references it is output.

	  The use_count of a vertex is the number of un-finalized tetrahedra 
	  referencing it. A vertex is finalized as soon as its use_count goes
    down to zero

    A vertex is deallocated whenever it is topologically finalized. 
    A tetrahedron is deallocated whenever all its vertices have been
    topologically finalized.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
    yuanxin liuy@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:

    21 January 2006 -- finalize infinite tets with base triangle on bounding box
    03 January 2006 -- finalization distinguishes finalized leaf cells and parents
    31 December 2005 -- finalization of new tetrahedra happens now in selective chunks
    28 December 2005 -- finalized tets are now stored (and referenced) by the grid
    25 September 2005 -- adapted from the 2D code after two chicken were killed
  
===============================================================================
*/

#ifndef SPDELAUNAY3D_H
#define SPDELAUNAY3D_H

#include "spwriter.h"
#include "svwriter.h"
#include "sscontainer3d.h"
#include "delaunay3.h"

class SPdelaunay3D : public SPwriter
{
public:

  // SPwriter interface
  void add_comment(const char* comment){};

  void set_npoints(int npoints){};
  void set_datatype(SPdatatype datatype){};
  void set_boundingbox(const double * bb_min_d, const double* bb_max_d);
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f);
  void set_boundingbox(const int* bb_min_i, const int* bb_max_i){};
  void set_finalizemethod(SPfinalizemethod finalizemethod);

  void write_header(){};

  void write_point(const double* p_pos_d){};
  void write_point(const float* p_pos_f);
  void write_point(const int* p_pos_i){};

  void write_finalize_cell(int idx);

  void close();

  ~SPdelaunay3D(){};

  // SPdelaunay3D interface

  bool open(SVwriter* svwriter = 0);

  int active_tetrahedra_next;
  void getActiveTetrahedraInit();
  int getActiveTetrahedraNext(float* v0, float* v1, float* v2, float* v3);
  bool getActiveTetrahedraCoords(int idx, float* v0, float* v1, float* v2, float* v3);
  bool getActiveTetrahedraSphere(int idx, float* rad_cen);
  bool getActiveTetrahedraCellIdx(int idx, int* cell_idx);

  int infinite_tetrahedra_next;
  void getInfiniteTetrahedraBaseTriangleInit();
  int getInfiniteTetrahedraBaseTriangleNext(float* v0, float* v1, float* v2);

  SPdelaunay3D();

  // SPdelaunay3D variables

  SVwriter* svwriter;
  SScontainer3D* ss3d;
  Delaunay3* dt;
};

#endif
