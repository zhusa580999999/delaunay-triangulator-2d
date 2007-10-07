/*
===============================================================================

  FILE:  spdelaunay2d.h
  
  CONTENTS:
  
    Expects a streaming 2D point cloud as input and produces a 
    streaming Delaunay triangle mesh as output. The Delaunay
    triangulation is computed incrementally as points stream
    in using the standard "locate" and "update: approach.
    
    We locate the triangle containing the next point by walking
    there. We start the walk at the triangle that was created
    last. This walk may fail if the walk leads us through an
    already finalized part of the triangulation. In this case
    we search brute force over all triangles.

    We update the triangulation to reestablish the Delaunay
    property using the Bowyer-Watson method that deletes all
    violated triangles and reconnects the resulting horizon
    of edges to the new point.

    The main new part is the concept of finalizing parts of
    the alreay computed Delaunay triangulation and outputting
    all finalized Delaunay triangles and their vertices in
    form of a streaming mesh.

    A Delaunay triangle can be finalized when its circumcircle
    is completely inside the spatially finalized region. It is
    then guaranteed to persist and may already be output and
    deallocated. A vertex is output just before the first triangle
    that references it and is topologically finalized in the moment
    is has no more triangles. In this moment it can be deallocated.
	  For this we keep a use_count with each such active vertex that
    at any time reflects the number of active triangles referencing
    it. As soon as this use_count goes down to zero the vertex may
    be finalized.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
    yuanxin liuy@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:

    21 January 2006 - delete infinite tris adjacent to finite & finalized tris
    30 August 2005 - added leo's fix for finding all duplicate points
    29 August 2005 - early intersect by storing cell idxs with the triangles
    21 August 2005 - improved handling of degenerate case for infinite circles.
    19 August 2005 - gotten rid of 'zombie' triangles.
    16 August 2005 - removed unnecessary incircle checks.
    15 August 2005 - significancly improved point location using walk when possible.
    14 August 2005 - more efficient intersection tests using tree structure
    12 August 2005 - initial prototype completed.
    07 August 2005 - adapted from Leo's code to support *real* streaming
    22 July 2005 -- last update by Yuanxin Liu
  
===============================================================================
*/

#ifndef SPDELAUNAY2D_H
#define SPDELAUNAY2D_H

#include "spwriter.h"
#include "smwriter.h"
#include "sscontainer2d.h"
#include "delaunay2.h"

class SPdelaunay2D : public SPwriter
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

  void write_point(const double* p_pos_d);
  void write_point(const float* p_pos_f);
  void write_point(const int* p_pos_i){};

  void write_finalize_cell(int idx);

  void close();

  ~SPdelaunay2D(){};

  // SPdelaunay2D interface
  bool open(SMwriter* smwriter = 0);

  int active_triangle_next;
  void getActiveTrianglesInit();
  int getActiveTrianglesNext(float* v0, float* v1, float* v2);
  bool getActiveTrianglesIdx(int idx, float* v0, float* v1, float* v2);
  bool getActiveTrianglesIdx(int idx, float* rad_cen);

  int infinite_triangle_next;
  void getInfiniteTrianglesInit();
  bool getInfiniteTrianglesNext(float* v0, float* v1, float* v2);

  SPdelaunay2D();

  // SPdelaunay2D variables

  SMwriter* smwriter;
  SScontainer2D* ss2d;
  Delaunay2* dt;
  int p_count_init;
};

#endif
