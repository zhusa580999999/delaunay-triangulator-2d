/*
===============================================================================

  FILE:  delaunay3.h
  
  CONTENTS:
  
    This is originally Leo's code and it is getting slowly modified for speed,
    memory efficiency, correctness, and readability.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
    yuanxin liuy@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:

     9 January 2006 -- adaptive offset for circumsphere counters floating-point error 
     5 January 2006 -- before going brute force we allow a second try with a grid tet
    31 December 2005 -- changed allocation of Delaunay Tetrahedron to happen in chunks
    24 December 2005 -- reworked the code and optimized bits and pieces everywhere
    14 December 2005 -- integrating Leo's new code that has everything as one struct
    20 October 2005 -- delay computation of sphere until sphere is needed
    28 September 2005 -- fewer brute_force searches by replacing "<=" with "<"
    27 September 2005 -- inserted compile switch counters for statistics
    26 September 2005 -- memory management for vertices and initial profiling
    25 September 2005 -- adapted from Leo's code after two chicken were killed
  
===============================================================================
*/

#ifndef DELAUNAY3_H
#define DELAUNAY3_H

#include <stdio.h>
#include "sscontainer3d.h"

#define D3_NULL_INDEX -19999999
#define D3_TETRA(corner) ((corner) >> 2)
#define D3_INDEX(corner) ((corner) & 3)

typedef struct Delaunay3Vertex
{
  // used for memory management
  Delaunay3Vertex* buffer_next;

  // input coordinates
	double x[3];

  // lifted coordinate (stored for efficiency)
  double sq;

	// for topological finalization
	int ref_count;

  // for streaming output
  int index;
} Delaunay3Vertex; 

#define RAD2_DEAD rad2
#define RAD2_DEAD_TRUE -2.0f
#define RAD2_DEAD_FALSE -1.0f

class Delaunay3Tetrahedron : public SSsphere
{
 public:
	Delaunay3Vertex* V[4];
	int N[4]; //opposite corner
};

class Delaunay3Chunk
{
 public:
  int available;
  int next;
  int full;
};

#define DELAUNAY3TETRAHEDRON_CHUNK_SIZE  1024
//#define DELAUNAY3TETRAHEDRON_CHUNK_EMPTY 512
#define DELAUNAY3TETRAHEDRON_CHUNK_EMPTY 756
#define DELAUNAY3TETRAHEDRON_CHUNK_ALLOC 16

class Delaunay3
{
public:
	// efficient memory management of tetrahedra
  int* tetrahedron_newchunks;
  int tetrahedron_newchunk_number;
  int tetrahedron_newchunk_alloc;

  Delaunay3Chunk* tetrahedron_chunks;
  int tetrahedron_chunk_current;
  int tetrahedron_chunk_alloc;
  int tetrahedron_chunk_available;

  Delaunay3Tetrahedron* tpool;
  int tetrahedron_buffer_used;
  int tetrahedron_buffer_maxsize;
 
  int allocTetrahedron();
	void deallocTetrahedron(int idx);

  void resetDelaunay3TetrahedronNewChunks();

  // efficient memory management of vertices

  Delaunay3Vertex* delaunay_vertex_buffer;
  int delaunay_vertex_buffer_size;
  int delaunay_vertex_buffer_maxsize;
  int delaunay_vertex_buffer_alloc;
	
  Delaunay3Vertex* allocVertex();
  Delaunay3Vertex* allocVertex(const float* xyz);
  Delaunay3Vertex* allocVertex(float x, float y, float z);
  void deallocVertex(Delaunay3Vertex* v);

  // the point at infinity
  Delaunay3Vertex* pinf;

  // does the tetrahedron have an infinite vertex
	bool isInf(int t);

  // point location variables
	int t_start; // the seed tetra for the point location walk
	
  // initialize incremental contruction with four (non-cocircular) vertices
  bool initialize(Delaunay3Vertex* V[]); 

  // insert point after point
  void insert(Delaunay3Vertex* p);

  // find a tetrahedra whose circumsphere contains the point
	int locate(Delaunay3Vertex* v);

  // helper routines for locate
  bool DuplicateVertex(Delaunay3Tetrahedron& t, Delaunay3Vertex* p);
	bool intersectTriangle(int c, Delaunay3Vertex* q, Delaunay3Vertex* p);
  int otherStartForLocate(Delaunay3Vertex* v);
  bool secondTry;
  int bruteLocate(Delaunay3Vertex* v); 

  //the insphere predicates with and without perturbation
	//the return value is the insphere determinant
	//it is negative, zero, or positive, when v is inside, on or outside the sphere.
	double inSphere(int t, Delaunay3Vertex* v);      // with perturbation, never 0
	double inSphereExact(int t, Delaunay3Vertex* v); // without perturbation

  // debug code (checks geometry and connectivity of the mesh)
 	void audit();
  void cornerPrint(int c);
	void cornerPrint4(int c);

	// constructor and destructor
	Delaunay3(); 
  ~Delaunay3();
};

#endif


