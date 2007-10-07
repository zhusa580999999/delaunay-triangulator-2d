/*
===============================================================================

  FILE:  delaunay3.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2005  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "delaunay3.h"
 
extern double insphere(double* pa, double* pb, double* pc, double* pd, double* pe);
extern double orient3d(double* pa, double* pb, double* pc, double* pd);
extern void exactinit();

//#define AUDIT_INSERT

#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

#ifdef COLLECT_STATISTICS
static int tetrahedron_buffer_maxused = 0;

static int stat_duplicate_points = 0;

static int stat_new_tetrahedra = 0;
static int stat_deleted_tetrahedra = 0;
extern int stat_delayed_deleted_tetrahedra;
static int stat_delayed_deleted_tetrahedra_max = 0;

static int stat_locate_calls = 0;
static int stat_t_start_was_dead = 0;
static int stat_t_start_no_other = 0;
static int stat_locate_after_0walk = 0;
static int stat_locate_after_1walk = 0;
static int stat_locate_active_run = 0;
static int stat_locate_active_maxrun = 0;
static __int64 stat_locate_active_checked = 0;
static int stat_locate_failed_steps = 0;
static int stat_locate_failed_deleted = 0;
static int stat_locate_failed_infinite = 0;
static int stat_locate_failed_noexit = 0;
static int stat_locate_failed_notsurewhy = 0;

static int stat_locateold_calls = 0;
static __int64 stat_locateold_active_total = 0;
static __int64 stat_locateold_active_checked = 0;
#endif

//----------- macro functions -------------------------------------------------------

//  assert macro
#ifndef NDEBUG
#define ASSERT(bool, string) if (!(bool)) {\
  fprintf(stderr,"\nASSERT FAILED (line %d of %s ): %s\n", \
__LINE__, __FILE__, string); }
#else
#define ASSERT(bool, string) 
#endif

// Stack data structure operations
#define STACKMAX 50000
#define POP(stack) (stack##st[stack##sp--])
#define isEMPTY(stack) (stack##sp < 0)
#define stkINIT(stack) {stack##sp = -1; }

#ifdef NDEBUG
#define PUSH(value, stack) { stack##st[++stack##sp] = value; }
#define stkDECLARE(stack,stn) static int stack##sp, stack##st[STACKMAX];    
#else
#define PUSH(value, stack) { \
	stack##st[++stack##sp] = value; \
	if (stack##max < stack##sp) { stack##max = stack##sp; \
	if (stack##max >= STACKMAX) { \
	fprintf(stderr,"ERROR: overflow stack ## %d pushing %d",  stack##max, value); exit(EXIT_FAILURE); } } /**/\
}
#define stkDECLARE(stack,stn) int stack##sp, stack##st[STACKMAX]; int stack##max; //AUDIT /**/ 
#endif

#define D3_CORNER(tetra,index) (((tetra)<<2)+(index))
#define D3_LIFT(x,y,z) (x*x+y*y+z*z)
#define COMPLEMENT(i) (-((i)+1)) 
#define DET2(p,q,i,j)((i##p)*(j##q) - (j##p)*(i##q))

// InSphere dot products
#define spdot(sp,pv,sv) ((sp)->x*((pv)->x[0]-(sv)->x[0])+(sp)->y*((pv)->x[1]-(sv)->x[1])\
         +(sp)->z*((pv)->x[2]-(sv)->x[2])+(sp)->sq*((pv)->sq-(sv)->sq))
 
#define SWAP(a,b,t) {t=a; a=b; b=t;}  
 
//----------- global variables and constants  -----------------------------------

#define DUPLICATE_POINT	-1
#define LOCATE_FAIL			-2

//the following two tables initialv and initialopp are used 
//to set up the initial corner table, which represents 
//five tetra, the first one finite, the other four infinite
static const int initialv[2][5][4] = {{{1,2,3,4}, {2,0,3,4}, {0,1,3,4}, {1,0,2,4}, {0,1,2,3}},
                                      {{0,2,3,4}, {2,1,3,4}, {1,0,3,4}, {0,1,2,4}, {1,0,2,3}}};

static const int initialopp[5][4] = {
  {D3_CORNER(1,1), D3_CORNER(2,0), D3_CORNER(3,1), D3_CORNER(4,0)}, 
  {D3_CORNER(2,1), D3_CORNER(0,0), D3_CORNER(3,0), D3_CORNER(4,1)},  
  {D3_CORNER(0,1), D3_CORNER(1,0), D3_CORNER(3,2), D3_CORNER(4,2)}, 
  {D3_CORNER(1,2), D3_CORNER(0,2), D3_CORNER(2,2), D3_CORNER(4,3)}, 
  {D3_CORNER(0,3), D3_CORNER(1,3), D3_CORNER(2,3), D3_CORNER(3,3)}};

// Tetrahedron manipulation tables

// c+offset[i][D3_INDEX(c)] advances c to (c+i)mod5
static const short offset[4][4] = {/*0*/{0,0,0,0}, /*1*/{1,1,1,-3}, /*2*/{2,2,-2,-2}, /*3*/{3,-1,-1,-1}};
// drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side.
static const short drop[4][3] = {{2,1,3}, {0,2,3}, {1,0,3}, {0,1,2}};
// offdr[i] contains drop(i)-index(i)
static const short offdr[4][3] = {{2,1,3}, {-1,1,2}, {-1,-2,1}, {-3,-2,-1}};
// invdrop[i][k] = j whenever drop[i][j] = k.  4s signal i=k; bad because i is dropped. 
static const short invdrop[4][4] = {{4,1,0,2}, {0,4,1,2}, {1,0,4,2}, {0,1,2,4}};

static short indoff[4][4][4] 
= {/*        0ABCD       1ABCD         2ACBD          3ABCD     */
	/* 0: */ {{5,5,5,5},  { 0,-1, 1, 2},  { 0,-1,-2, 1},  { 0,-3,-2,-1}},
	/* 1: */ {{1,0,2,3},  { 5, 5, 5, 5},  {-2, 0,-1, 1},  {-2, 0,-3,-1}}, 
	/* 2: */ {{2,1,0,3},  {-1, 1, 0, 2},  {-1,-2, 0, 1},  {-3,-2, 0,-1}}, 
	/* 3: */ {{1,2,3,0},  { 1,-1, 2, 0},  {-2,-1, 1, 0},  {-2,-3,-1, 0}}}; 
  
stkDECLARE(dfs, "dfs")  // DFS stack to find dead tetras
stkDECLARE(nhbr, "nhbr") // stack for dead corners with live neighbors
stkDECLARE(kill, "kill") // stack for base corners of tetras to recycle		

// efficient memory allocation

int Delaunay3::allocTetrahedron()
{
  // do we need to switch chunks or allocate more chunks ?
	if (tetrahedron_chunks[tetrahedron_chunk_current].next == D3_NULL_INDEX)
  {
    // set the full indicator of the current chunk to full
    tetrahedron_chunks[tetrahedron_chunk_current].full = DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
    // are there more chunks available
    if (tetrahedron_chunk_available == D3_NULL_INDEX)
    {
      // we need to allocate more chunks ... this will be the index of the next current chunk
      tetrahedron_chunk_current = tetrahedron_chunk_alloc;
      // and this will be the index of the next available chunk
      tetrahedron_chunk_available = tetrahedron_chunk_current+1;
      // double the number of chunks
		  tetrahedron_chunk_alloc = (3*tetrahedron_chunk_alloc)/2;
      // reallocate the array of tetrahedra
		  tpool = (Delaunay3Tetrahedron*) realloc(tpool, tetrahedron_chunk_alloc*DELAUNAY3TETRAHEDRON_CHUNK_SIZE*sizeof(Delaunay3Tetrahedron));
		  if (tpool == NULL)
      {
			  fprintf(stderr,"ERROR: cannot reallocate memory for %d chunks of %d tetrahedra", tetrahedron_chunk_alloc, DELAUNAY3TETRAHEDRON_CHUNK_SIZE); 
			  exit(EXIT_FAILURE);
		  }
      // reallocate the array of chunks
		  tetrahedron_chunks = (Delaunay3Chunk*) realloc(tetrahedron_chunks, tetrahedron_chunk_alloc*sizeof(Delaunay3Chunk));
		  if (tetrahedron_chunks == NULL)
      {
			  fprintf(stderr,"ERROR: cannot reallocate memory for %d tetrahedron chunks ", tetrahedron_chunk_alloc); 
			  exit(EXIT_FAILURE);
		  }
      // link the chunks into a single-linked list
      int i,k,j;
		  for (i = tetrahedron_chunk_current; i < tetrahedron_chunk_alloc; i++)
      {
        // make a link from the chunk to the first tetrahedron
        k = i*DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
        tetrahedron_chunks[i].next = k; 
        // link the tetrahedra of each chunk into a single-linked list
  		  for (j = 0; j < DELAUNAY3TETRAHEDRON_CHUNK_SIZE; j++)
        {
  			  tpool[k].RAD2_DEAD = RAD2_DEAD_TRUE;
	  		  tpool[k].next = k+1;
          k++;
        }
        tpool[k-1].next = D3_NULL_INDEX;
        // link the chunks
        tetrahedron_chunks[i].available = i+1;
      }
      tetrahedron_chunks[i-1].available = D3_NULL_INDEX; 
    }
    else
    {
      // we simply switch to the next available chunk
      tetrahedron_chunk_current = tetrahedron_chunk_available;
      tetrahedron_chunk_available = tetrahedron_chunks[tetrahedron_chunk_current].available;
    }
    // set the full indicator of the new chunk to empty
    tetrahedron_chunks[tetrahedron_chunk_current].full = -1;
    // assert that we do have a tetrahedron now
    assert(tetrahedron_chunks[tetrahedron_chunk_current].next != D3_NULL_INDEX);
    // add the new chunk of tetrahedra to the list of new chunks
    if (tetrahedron_newchunk_number == tetrahedron_newchunk_alloc)
    {
      // need to alloc more memory first
      tetrahedron_newchunk_alloc = tetrahedron_newchunk_alloc * 2;
  	  tetrahedron_newchunks = (int*)realloc(tetrahedron_newchunks, tetrahedron_newchunk_alloc*sizeof(int));
    }
    tetrahedron_newchunks[tetrahedron_newchunk_number] = tetrahedron_chunk_current;
    tetrahedron_newchunk_number++;
  }

  // get the next free tetrahedron

  int tetrahedron = tetrahedron_chunks[tetrahedron_chunk_current].next;  
  tetrahedron_chunks[tetrahedron_chunk_current].next = tpool[tetrahedron].next;

  // make sure it is dead

  assert(tpool[tetrahedron].RAD2_DEAD == RAD2_DEAD_TRUE);
  
  // initialize it
  tpool[tetrahedron].next = -1;
	tpool[tetrahedron].RAD2_DEAD = RAD2_DEAD_FALSE;
  tpool[tetrahedron].cell = 0;

#ifdef COLLECT_STATISTICS
  tetrahedron_buffer_used++;
  if (tetrahedron_buffer_used > tetrahedron_buffer_maxused) tetrahedron_buffer_maxused = tetrahedron_buffer_used;
#endif

  if (tetrahedron > tetrahedron_buffer_maxsize) tetrahedron_buffer_maxsize = tetrahedron;

  return tetrahedron;
}

void Delaunay3::deallocTetrahedron(int idx)
{
  // this sucker is dead
  tpool[idx].RAD2_DEAD = RAD2_DEAD_TRUE;
  // get the chunk that this tetrahedron is part of 
  int chunk = idx / DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
  // link it back into this chunk
  tpool[idx].next = tetrahedron_chunks[chunk].next;
  tetrahedron_chunks[chunk].next = idx;
  // decrement this chunks full counter
  tetrahedron_chunks[chunk].full--;
  // is this chunk empty enough to be used for allocating
  if (tetrahedron_chunks[chunk].full == DELAUNAY3TETRAHEDRON_CHUNK_EMPTY)
  {
    // then put it into the list of available chunks for allocating
    tetrahedron_chunks[chunk].available = tetrahedron_chunk_available;
    tetrahedron_chunk_available = chunk;
  }
#ifdef COLLECT_STATISTICS
  tetrahedron_buffer_used--;
#endif
}

void Delaunay3::resetDelaunay3TetrahedronNewChunks()
{
  tetrahedron_newchunks[0] = tetrahedron_chunk_current;
  tetrahedron_newchunk_number = 1;
}

Delaunay3Vertex* Delaunay3::allocVertex()
{
  if (delaunay_vertex_buffer == 0)
  {
    delaunay_vertex_buffer = (Delaunay3Vertex*)malloc(sizeof(Delaunay3Vertex)*delaunay_vertex_buffer_alloc);
    if (delaunay_vertex_buffer == 0)
    {
      fprintf(stderr,"malloc for delaunay vertex buffer failed\n");
      return 0;
    }
    for (int i = 0; i < delaunay_vertex_buffer_alloc; i++)
    {
      delaunay_vertex_buffer[i].buffer_next = &(delaunay_vertex_buffer[i+1]);
    }
    delaunay_vertex_buffer[delaunay_vertex_buffer_alloc-1].buffer_next = 0;
    delaunay_vertex_buffer_alloc = 2*delaunay_vertex_buffer_alloc;
  }
  // get pointer to next available delaunay vertex
  Delaunay3Vertex* vertex = delaunay_vertex_buffer;
  delaunay_vertex_buffer = vertex->buffer_next;

  // init delaunay vertex
  vertex->index = -1;
  vertex->ref_count = 0;

#ifdef COLLECT_STATISTICS
  delaunay_vertex_buffer_size++; if (delaunay_vertex_buffer_size > delaunay_vertex_buffer_maxsize) delaunay_vertex_buffer_maxsize = delaunay_vertex_buffer_size;
#endif

  return vertex;
}

Delaunay3Vertex* Delaunay3::allocVertex(const float* xyz)
{
  Delaunay3Vertex* vertex = allocVertex();

  vertex->x[0] = xyz[0];
  vertex->x[1] = xyz[1];
  vertex->x[2] = xyz[2];	
  vertex->sq = D3_LIFT(vertex->x[0],vertex->x[1],vertex->x[2]);

  return vertex;
}

Delaunay3Vertex* Delaunay3::allocVertex(float x, float y, float z)
{
  Delaunay3Vertex* vertex = allocVertex();

  vertex->x[0] = x;
  vertex->x[1] = y;
  vertex->x[2] = z;	
  vertex->sq = D3_LIFT(vertex->x[0],vertex->x[1],vertex->x[2]);

  return vertex;
}

void Delaunay3::deallocVertex(Delaunay3Vertex* vertex)
{
  vertex->buffer_next = delaunay_vertex_buffer;
  delaunay_vertex_buffer = vertex;
#ifdef COLLECT_STATISTICS
  delaunay_vertex_buffer_size--;
#endif
}

// a = (b1+b2+b3)/3
static inline void Average(Delaunay3Vertex* a, Delaunay3Vertex* b1, Delaunay3Vertex* b2, Delaunay3Vertex* b3)
{
  a->x[0] = (b1->x[0] + b2->x[0]+ b3->x[0])/3.0;
  a->x[1] = (b1->x[1] + b2->x[1]+ b3->x[1])/3.0;
  a->x[2] = (b1->x[2] + b2->x[2]+ b3->x[2])/3.0;
	a->sq = D3_LIFT(a->x[0], a->x[1],a->x[2]);
}
 
// signed volume of the parallel piped given by v0..v3
double volume(Delaunay3Vertex* v0, Delaunay3Vertex* v1, Delaunay3Vertex* v2, Delaunay3Vertex* v3)
{
	return orient3d(v0->x, v1->x, v2->x, v3->x); 
}

// inSphere predicate without perturbation
// the insphere determinant is always returned
double Delaunay3::inSphereExact(int t, Delaunay3Vertex* v)
{
	double d;

  assert(tpool[t].V[2] != pinf && tpool[t].V[3] != pinf);

  if (tpool[t].V[0] == pinf) // if v0 is infinite
  {
    assert(tpool[t].V[1] != pinf);
		d = -orient3d(v->x, tpool[t].V[1]->x, tpool[t].V[2]->x, tpool[t].V[3]->x);
  }
  else if (tpool[t].V[1] == pinf) // if v1 is infinite
  {
		d = -orient3d(tpool[t].V[0]->x, v->x, tpool[t].V[2]->x, tpool[t].V[3]->x);
  }
  else // finite sphere
  { 
		d = -insphere(tpool[t].V[0]->x, tpool[t].V[1]->x, tpool[t].V[2]->x, tpool[t].V[3]->x, v->x);									
	}

  return d;
}

// given four vertices one of which (V0 here) is infinite, return the normal of the plane
static void Normal_V0(Delaunay3Vertex * v0, Delaunay3Vertex * v1, Delaunay3Vertex * v2, Delaunay3Vertex * pv, float N[])
{
  double x0, y0, z0, sq0, x1, y1, z1, sq1, x2, y2, z2, sq2;
  // v0 is infinite
  x0 = v0->x[0]; y0 = v0->x[1]; z0 = v0->x[2]; sq0 = v0->sq;
  x1 = v1->x[0] - pv->x[0]; y1 = v1->x[1] - pv->x[1]; z1 = v1->x[2] - pv->x[2]; sq1 = v1->sq - pv->sq;
  x2 = v2->x[0] - pv->x[0]; y2 = v2->x[1] - pv->x[1]; z2 = v2->x[2] - pv->x[2]; sq2 = v2->sq - pv->sq;
  // 2x2 minors
  double xy = DET2(0,1,x,y); 
  double xz = DET2(0,1,x,z); 
  double yz = DET2(0,1,y,z); 
  double xs = DET2(0,1,x,sq); 
  double ys = DET2(0,1,y,sq);
  double zs = DET2(0,1,z,sq);
  N[0] = (float)(-y2*zs + z2*ys - sq2*yz);
  N[1] = (float)( x2*zs - z2*xs + sq2*xz);
  N[2] = (float)(-x2*ys + y2*xs - sq2*xy);
}

// given four vertices one of which (V1 here) is infinite, return the normal of the plane
static void Normal_V1(Delaunay3Vertex * v0, Delaunay3Vertex * v1, Delaunay3Vertex * v2, Delaunay3Vertex * pv, float N[])
{
  double x0, y0, z0, sq0, x1, y1, z1, sq1, x2, y2, z2, sq2;
  // v1 is infinite
  x0 = v0->x[0] - pv->x[0]; y0 = v0->x[1] - pv->x[1]; z0 = v0->x[2] - pv->x[2]; sq0 = v0->sq - pv->sq;
  x1 = v1->x[0]; y1 = v1->x[1]; z1 = v1->x[2]; sq1 = v1->sq;
  x2 = v2->x[0] - pv->x[0]; y2 = v2->x[1] - pv->x[1]; z2 = v2->x[2] - pv->x[2]; sq2 = v2->sq - pv->sq;
  // 2x2 minors
  double xy = DET2(0,1,x,y); 
  double xz = DET2(0,1,x,z); 
  double yz = DET2(0,1,y,z); 
  double xs = DET2(0,1,x,sq); 
  double ys = DET2(0,1,y,sq);
  double zs = DET2(0,1,z,sq);
  N[0] = (float)(-y2*zs + z2*ys - sq2*yz);
  N[1] = (float)( x2*zs - z2*xs + sq2*xz);
  N[2] = (float)(-x2*ys + y2*xs - sq2*xy);
}

// inSphere predicate with perturbation, if no perturbation
// is applied, it returns the inSphere determinant with respect to t
// it might return the insphere determinant with respect to a constructed sphere
double Delaunay3::inSphere(int t, Delaunay3Vertex* v)
{
  double d;	 
	Delaunay3Tetrahedron& tet = tpool[t];

  if (tet.V[0] == pinf) // infinite circle because of v0
  {
		d = -orient3d(v->x, tet.V[1]->x, tet.V[2]->x, tet.V[3]->x);
  	if (d == 0)
    {
  		float N[3];
  		Normal_V0(tet.V[0], tet.V[1], tet.V[2], tet.V[3],N);
  		double x[3];
      x[0] = v->x[0] - N[0];
      x[1] = v->x[1] - N[1];
      x[2] = v->x[2] - N[2];
  		d = -insphere(x, tet.V[1]->x, tet.V[2]->x, tet.V[3]->x, v->x);   
    }
  }
  else if (tet.V[1] == pinf) // infinite circle because of v1
  {
		d = -orient3d(tet.V[0]->x, v->x, tet.V[2]->x, tet.V[3]->x);
    if (d == 0)
    {
  		float N[3];
  		Normal_V1(tet.V[0], tet.V[1], tet.V[2], tet.V[3],N);
  		double x[3];
      x[0] = v->x[0] - N[0];
      x[1] = v->x[1] - N[1];
      x[2] = v->x[2] - N[2];
  		d = -insphere(tet.V[0]->x, x, tet.V[2]->x, tet.V[3]->x, v->x);   
    }
  }
  else // finite circle
  {
    d = -insphere(tet.V[0]->x, tet.V[1]->x, tet.V[2]->x, tet.V[3]->x,v->x);
  }
  return d;
}

Delaunay3::Delaunay3()
{
  // initialize the vertex pool
	delaunay_vertex_buffer = 0;
	delaunay_vertex_buffer_size = 0;
	delaunay_vertex_buffer_maxsize = 0;
	delaunay_vertex_buffer_alloc = 512;

  // initialize the tetrahedra pool
#ifdef COLLECT_STATISTICS
  tetrahedron_buffer_used = 0;
  tetrahedron_buffer_maxused = 0;
#endif

  tetrahedron_buffer_maxsize = 0;

  // the initial number of chunks of tetrahedra
	tetrahedron_chunk_alloc = DELAUNAY3TETRAHEDRON_CHUNK_ALLOC;
  // allocate the array of tetrahedra
	tpool = (Delaunay3Tetrahedron*)malloc(tetrahedron_chunk_alloc*DELAUNAY3TETRAHEDRON_CHUNK_SIZE*sizeof(Delaunay3Tetrahedron));
	if (tpool == NULL)
  {
		fprintf(stderr,"ERROR: cannot allocate memory for %d chunks of %d tetrahedra", tetrahedron_chunk_alloc, DELAUNAY3TETRAHEDRON_CHUNK_SIZE); 
		exit(EXIT_FAILURE);
	}
  // allocate the array of chunks
	tetrahedron_chunks = (Delaunay3Chunk*)malloc(tetrahedron_chunk_alloc*sizeof(Delaunay3Chunk));
	if (tetrahedron_chunks == NULL)
  {
		fprintf(stderr,"ERROR: cannot allocate memory for %d tetrahedron chunks ", tetrahedron_chunk_alloc); 
		exit(EXIT_FAILURE);
	}
  // the initial allocation for the newchunk index array
	tetrahedron_newchunk_alloc = 16;
  // the initial number of elements in the newchunk index array
	tetrahedron_newchunk_number = 0;
  // allocate the newchunk index array of new chunk indices
  tetrahedron_newchunks = (int*)malloc(tetrahedron_newchunk_alloc*sizeof(int));
  // link the chunks into a single-linked list
  int i,k,j;
	for (i = 0; i < tetrahedron_chunk_alloc; i++)
  {
    // make a link from the chunk to the first tetrahedron
    k = i*DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
    tetrahedron_chunks[i].next = k;
    // link the tetrahedra of each chunk into a single-linked list
  	for (j = 0; j < DELAUNAY3TETRAHEDRON_CHUNK_SIZE; j++)
    {
  		tpool[k].RAD2_DEAD = RAD2_DEAD_TRUE;
	  	tpool[k].next = k+1;
      k++;
    }
    tpool[k-1].next = D3_NULL_INDEX;
    // link the chunks
    tetrahedron_chunks[i].available = i+1; 
  }
  tetrahedron_chunks[i-1].available = D3_NULL_INDEX; 
  // initially the first chunk of tetrahedra
  tetrahedron_chunk_current = 0;
  // initially the next available chunk of tetrahedra
  tetrahedron_chunk_available = 1;

	// initialize the point at infinity
	pinf = allocVertex(0,0,0);
  pinf->sq = 1;

	exactinit();
}

Delaunay3::~Delaunay3()
{
#ifdef COLLECT_STATISTICS
  fprintf(stderr, "delayed_deleted_max %d (%d) max_used %d (%d) alloced %d\n", stat_delayed_deleted_tetrahedra_max, stat_delayed_deleted_tetrahedra, tetrahedron_buffer_maxused, tetrahedron_buffer_used, tetrahedron_chunk_alloc*DELAUNAY3TETRAHEDRON_CHUNK_SIZE);
  fprintf(stderr, "new_tetrahedra %d (%.1f) deleted_tetrahedra %d (%.1f)\n", stat_new_tetrahedra,(float)stat_new_tetrahedra/(stat_locate_calls-stat_duplicate_points), stat_deleted_tetrahedra,(float)stat_deleted_tetrahedra/(stat_locate_calls-stat_duplicate_points));
  fprintf(stderr, "locate_calls %d (t_start dead %d no other %d)\n", stat_locate_calls, stat_t_start_was_dead, stat_t_start_no_other);
  if (stat_locate_calls)
  {
    char buffer1[128];
    _i64toa(stat_locate_active_checked, buffer1, 10);
    float buffer2 = (((float)stat_locate_active_checked/stat_locate_calls*10))/10.0f;
    fprintf(stderr, "checked %s (%.1f) walk0 %d walk1 %d maxrun %d\n", buffer1, buffer2, stat_locate_after_0walk, stat_locate_after_1walk,stat_locate_active_maxrun);
    fprintf(stderr, "failed %d (%.1f %%) - steps %d - deleted %d - infinite %d - noexit %d - notsurewhy %d\n", stat_locate_failed_steps+stat_locate_failed_deleted+stat_locate_failed_infinite+stat_locate_failed_noexit+stat_locate_failed_notsurewhy, 100.0f*(stat_locate_failed_steps+stat_locate_failed_deleted+stat_locate_failed_infinite+stat_locate_failed_noexit+stat_locate_failed_notsurewhy)/stat_locate_calls, stat_locate_failed_steps, stat_locate_failed_deleted, stat_locate_failed_infinite, stat_locate_failed_noexit, stat_locate_failed_notsurewhy);
  }
  fprintf(stderr, "stat_locateold_calls %d\n", stat_locateold_calls);
  if (stat_locateold_calls)
  {
    char buffer1[128];
    _i64toa(stat_locateold_active_total, buffer1, 10);
    float buffer2 = (((float)stat_locateold_active_total/stat_locateold_calls*10))/10.0f;
    char buffer3[128];
    _i64toa(stat_locateold_active_checked, buffer3, 10);
    float buffer4 = (((float)stat_locateold_active_checked/stat_locateold_calls*10))/10.0f;
    fprintf(stderr, "total %s (%.1f) searched %s (%.1f)\n", buffer1, buffer2, buffer3, buffer4);
  }
  fprintf(stderr, "duplicate_points %d\n", stat_duplicate_points);
#endif 
}

bool Delaunay3::initialize(Delaunay3Vertex* v0[])
{
  int j, p;
	Delaunay3Vertex* V[5];
	V[0] = pinf;  V[0]->ref_count = 4;
	V[1] = v0[0]; V[1]->ref_count = 4;
	V[2] = v0[1]; V[2]->ref_count = 4;
	V[3] = v0[2]; V[3]->ref_count = 4;
	V[4] = v0[3]; V[4]->ref_count = 4;

  tetrahedron_chunks[tetrahedron_chunk_current].next = 5;
  tetrahedron_chunks[tetrahedron_chunk_current].full = -1;

  tetrahedron_newchunks[tetrahedron_newchunk_number++] = tetrahedron_chunk_current;

#ifdef COLLECT_STATISTICS
  tetrahedron_buffer_used = 5;
  tetrahedron_buffer_maxused = 5;
#endif

  tetrahedron_buffer_maxsize = 5;

  //----------------------------------------------------------------------------------
	//  Create the initial finite tetrahedron that is surrounded by 4 infinite ones                                                       
	//----------------------------------------------------------------------------------
  double d = volume(V[1],V[2],V[3],V[4]); // if d<0, then we need to swap

  if (d == 0.0)
  { 
		fprintf(stderr, "FATAL ERROR: initial tetra must not be flat\n");
    exit(0);
	} 

  for (p=0; p<5; p++)
  {
		tpool[p].next = -1;
		tpool[p].RAD2_DEAD = RAD2_DEAD_FALSE;
    tpool[p].cell = 0;
		
		for (j=0; j<4; j++)
    {
      // pay attention to orientation when assigning vertices
			tpool[p].V[j] = V[initialv[d<0][p][j]];
			tpool[p].N[j] = initialopp[p][j];
    } 

    tpool[p].rad = -1.0f;
	}

	// initialize the seed tetra for walk
	t_start = 0;

	// allow a second try if the first walk fails 
  secondTry = true;

	return true;
}

void Delaunay3::insert(Delaunay3Vertex* pv)
{
  // find an initial tetra whose circumsphere contains the point pv
  int t = locate(pv);

	if (t < 0)
  {
    if (t == DUPLICATE_POINT)
    {
#ifdef COLLECT_STATISTICS
      stat_duplicate_points++;
#endif
     return;
    }
    else if (t == LOCATE_FAIL)
    {
		  fprintf(stderr,"ERROR: fail to locate point %f %f %f\n", pv->x[0],pv->x[1],pv->x[2]);
		  exit(0);
	  }
    else
    {
		  fprintf(stderr,"ERROR: t = locate(pv) returned invalid index %d\n", t);
    }
  }

	// All tetrahedra that contain pv (or whose circumsphere contains pv) are
  // declared "dead" and pushed onto the kill stack. We perform a DFS using
  // a stack pst to find and kill these tetrahedra.
	// All corners along the live-dead boundary are saved onto the stack nhbr,
	// then make new tetras and hook in to live by setting the last opp pointer.
	//
	// Invariants/operations: Tetrahedron p is marked alive or dead on first visit.
	// Corner c is pushed on stack when D3_TETRA(s[c].opp) is marked dead. 
	//
	// On termination, stack nhbr contains dead corners. Each dead corner points 
	// to the negated new corner. Stack kill contains old tetrahedra for final recycling.

	stkINIT(dfs);  // DFS stack holds corners opposite dead tetras 
	stkINIT(nhbr); // stack for dead corners with live nhbr tetras
	stkINIT(kill); // stack of dead tetras to recycle

  // decrement the reference counters of the initial tetra t
	tpool[t].V[0]->ref_count--;
	tpool[t].V[1]->ref_count--;
	tpool[t].V[2]->ref_count--;
	tpool[t].V[3]->ref_count--;
  // before killing we must de-link this cell-linked tet
  if (tpool[t].cell) 
  {
    if (tpool[t].next == t) // is it the last in the list
    {
      assert(tpool[t].cell->data == t);
      tpool[t].cell->data = -1;
    }
    else
    {
      tpool[tpool[t].prev].next = tpool[t].next;
      tpool[tpool[t].next].prev = tpool[t].prev;
      if (tpool[t].cell->data == t)
      {
        tpool[t].cell->data = tpool[t].next;
      }
    }
  }
	// kill the initial tetra t
  PUSH(t, kill);
  tpool[t].RAD2_DEAD = RAD2_DEAD_TRUE;

#ifdef COLLECT_STATISTICS
  stat_deleted_tetrahedra++;
#endif

  //--------------------------------------------------------------------------------
	//          Search the tetra complex in depth first order.
  //--------------------------------------------------------------------------------
	//  After the dfs is completed, we have two stacks
	//   nhbr: dead corners opposite live corners or void corners (due to finalization)
	//   kill: dead tetras 
	//
	//  The corners in the nhbr stack correspond to triangle faces on the ``horizon".
	//  There are three kinds of corners that correspond to a horizon triangle. 
  //    c1: corner from a dead tetra (those in the nhbr)
  //    c2: corner from a live tetra (opposite c1 before insertion), which could be missing because of finalization
  //    c3: corner from a new tetra  (opposite c2 after insertion)
  //
	//  After the dfs is completed, 
	//    c1 points to COMPLEMENT(c2) ( the reason for this is that later we wish 
	//                           to walk inside the dead tetra and needs to know when to stop)
	//    c2 points to c3
	//    c3 points to c2  
	// 
	//  In addition, the vertex of c3 is set to be pv 
  //--------------------------------------------------------------------------------

  int c, ci;
  // initialize the dfs stack by pushing the four faces of t onto the stack
	for (ci = 0; ci < 4; ci++)
  {
		c = tpool[t].N[ci];
		if ( c != D3_NULL_INDEX )
		{
      // put existing neighbors onto the dfs stack
			PUSH(c, dfs);
		}
		else
    {
      // for non-existing (i.e. finalized) neighbors put a horizon face on stack
			PUSH(COMPLEMENT(D3_CORNER(t,ci)), dfs); 
		}
	}

	bool reachedHorizon;

	while (!isEMPTY(dfs))
  {
 		c = POP(dfs);

		// check if the horizon is reached, which happens either
		// 1) the corner c is negative, signaling that the tetrahedron has
		//    been finalized therefore cannot have a sphere with pv inside
		// or,
		// 2) pv is not inside the sphere of the tetra of c 
		if ( c < 0)
    {
			reachedHorizon = true;
    }
		else
    {
			ASSERT(tpool[D3_TETRA(tpool[D3_TETRA(c)].N[D3_INDEX(c)])].RAD2_DEAD == RAD2_DEAD_TRUE, "dfs stack element with non-dead neighbor");
			t = D3_TETRA(c);
			if (tpool[t].RAD2_DEAD == RAD2_DEAD_TRUE)
      {
        continue; // this tetrahedron is already dead 
      }
      ci = D3_INDEX(c);
			reachedHorizon = (inSphere(t, pv) > 0);
		}
 
		if (reachedHorizon == false) // kill and continue dfs if pv is strictly inside this tetrahedron
    {
      // decrement the reference counters of this tetra t
			tpool[t].V[0]->ref_count--;
			tpool[t].V[1]->ref_count--;
			tpool[t].V[2]->ref_count--;
			tpool[t].V[3]->ref_count--;

      // before killing we must de-link this cell-linked tet
    	if (tpool[t].cell) 
      {
        if (tpool[t].next == t) // is it the last in the list
        {
          assert(tpool[t].cell->data == t);
          tpool[t].cell->data = -1;
        }
        else
        {
          tpool[tpool[t].prev].next = tpool[t].next;
          tpool[tpool[t].next].prev = tpool[t].prev;
          if (tpool[t].cell->data == t)
          {
            tpool[t].cell->data = tpool[t].next;
          }
        }
      }
      // kill this tetra t
      PUSH(t, kill);

      tpool[t].RAD2_DEAD = RAD2_DEAD_TRUE; 
#ifdef COLLECT_STATISTICS
      stat_deleted_tetrahedra++;
#endif

      // loop over the three non-entry faces of t
			for (int k = 1; k <= 3; k++)
      {
				// nc is the new corner opposite the kth corner of t
				// if it is valid, we just push it on the dfs stack
				// otherwise, we negate the kth corner of t and push
				// it on the stack

				int nc = tpool[t].N[ci+offset[k][ci]];
				if (nc != D3_NULL_INDEX)
        {
          // put existing neighbors onto the dfs stack
					PUSH(nc, dfs);
				}
				else
        {
          // for non-existing (i.e. finalized) neighbors put a horizon face on stack
					PUSH(COMPLEMENT(c+offset[k][ci]), dfs);
				}
			}
		}
		else
    {
      // retrieve the corner opposite to c and store it in cDead
			// if c is valid (>=0), cDead is the opposite corner
      // otherwise, c is the corner we want after negation.
			int cDead = c < 0 ? COMPLEMENT(c): tpool[t].N[ci]; 

      //push the dead corner on the stack nhbr
      PUSH(cDead, nhbr);
	 
      // create the tetrahedron that will cover this horizon
      int newt = allocTetrahedron();

#ifdef COLLECT_STATISTICS
      stat_new_tetrahedra++;
#endif

      // get last corner of new tetra
			int newc = D3_CORNER(newt,3);
			tpool[newt].V[3] = pv;

			pv->ref_count++;

			if (c < 0) // was this horizon face a "finalization boundary"?
      {
				tpool[newt].N[3] = D3_NULL_INDEX;
			}
			else // otherwise connect the horizon face to the newly created tetrahedron
      {
				tpool[newt].N[3] = c;
				tpool[t].N[ci] = newc;
			}
			tpool[D3_TETRA(cDead)].N[D3_INDEX(cDead)] = COMPLEMENT(newc);
		}
	}

  //--------------------------------------------------------------------------------
	//         Walk in the tetra complex to hook up new neighbors
  //--------------------------------------------------------------------------------
  //  Each new tetra has one old neighbor, which it already knows, and
  //  three new neighbors that we will hook it up with.  
  //  For each corner cDead belonging to the tetra tDead from the nhbr
  //  stack, we consider the other three corners in tDead. For each of these
  //  corners c, we do a walk around the edge that does not use c and cDead
  //  until we find a negative corner. The fact that this is negative means
  //  we have just reached a horizon triangle so the complement of this
  //  negative corner gives the new tetra we want to hook up with.   
  //--------------------------------------------------------------------------------
	Delaunay3Vertex *v0, *v1, *v2; // pointers to vertices

	int newc;
	while (!isEMPTY(nhbr))
  {
		int cDead = POP(nhbr); 
		int jdead = D3_INDEX(cDead); //  dead tetra and index of dropped corner.
		int tDead = D3_TETRA(cDead);

		/*        fprintf(stderr,"--Popped %d(%d)\n", dead, jdead); /**/
		ASSERT(tpool[tDead].RAD2_DEAD == RAD2_DEAD_TRUE, "corner on nhbr stack is not dead!?");
		newc = COMPLEMENT(tpool[tDead].N[jdead]);
    int newi = 0;  //index of the new corner we try to hook up with a neighbor
		int newt = D3_TETRA(newc); 
 
		int t,off,nc, ni,i,j;

		cDead -= jdead; // just use base of dead one.
		// new tetra has 0,1,2,3=pv; 
		// corresponding old indices before jdead is dropped: 
		//   drop[j][0],..,drop[j][3], (no corresp to pv)
		j = jdead; 
		i = drop[jdead][0]; // old index of new corner 0;
		t = tDead;
		v0 = tpool[tDead].V[i];
		nc = tpool[tDead].N[i];

		//in each iteration, we jump from corner specified by
		// (t,i) to the corner nc. We also maintain j, the index
		// of the corner nc from the previous jump. The corners 
		// given by the js are the 1-ring around an edge
		while (nc>=0)
    { 
			ni = D3_INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
			j = ni; i = ni + off; t = D3_TETRA(nc + off);
			nc = tpool[t].N[i];// fix new j, i, c, and try neighbor
		} 
		nc = COMPLEMENT(nc); // go to new tetra
		assert(tpool[D3_TETRA(nc)].V[D3_INDEX(nc)]==pv);

		tpool[newt].V[newi] = v0;
		tpool[newt].N[newi] = nc-3+invdrop[i][j];

		v0->ref_count++; newi++; 

		j = jdead; 
		i = drop[jdead][1]; // old index of new corner 1;
		t = tDead;
		v1 = tpool[tDead].V[i]; // copy vertex v1
		nc = tpool[tDead].N[i];

		//in each iteration, we jump from corner specified by
		// (t,i) to the corner nc. We also maintain j, the index
		// of the corner nc from the previous jump. The corners 
		// given by the js are the 1-ring around an edge
		while (nc>=0)
    { 
			ni = D3_INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
			j = ni; i = ni + off; t = D3_TETRA(nc + off);
			nc = tpool[t].N[i];// fix new j, i, c, and try neighbor
		} 
		nc = COMPLEMENT(nc); // go to new tetra
		
		assert(tpool[D3_TETRA(nc)].V[D3_INDEX(nc)]==pv);

		tpool[newt].V[newi] = v1;
		tpool[newt].N[newi] = nc-3+invdrop[i][j];

		v1->ref_count++; newi++;  

		j = jdead; 
		i = drop[jdead][2]; // old index of new corner 2;
		t = tDead;
		v2 = tpool[tDead].V[i];
		nc = tpool[tDead].N[i]; // go to neighbor

		//in each iteration, we jump from corner specified by
		// (t,i) to the corner nc. We also maintain j, the index
		// of the corner nc from the previous jump. The corners 
		// given by the js are the 1-ring around an edge
		while (nc>=0)
    { 
			ni = D3_INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
			j = ni; i = ni + off; t = D3_TETRA(nc + off);
			nc = tpool[t].N[i];// fix new j, i, c, and try neighbor
		} 
		nc = COMPLEMENT(nc); // go to new tetra
		assert(tpool[D3_TETRA(nc)].V[D3_INDEX(nc)]==pv);

		tpool[newt].V[newi] = v2;
		tpool[newt].N[newi] = nc-3+invdrop[i][j];
		v2->ref_count++;  newi++;  

    tpool[newt].rad = -1.0f;
    
#ifndef NDEBUG
    if (!isInf(newt) && volume(v0,v1,v2,pv)<=0)
    {
			double d= volume(v0,v1,v2,pv);
			fprintf(stderr, "FATAL ERROR: made a negatively oriented tetra: volume %lg\n",d);
			exit(0);
		}
#endif
 	}

  // recycle memory of dead tetrahedra

	while (!isEMPTY(kill))
  {
		deallocTetrahedron(POP(kill));
	}

#ifdef COLLECT_STATISTICS
  if (stat_delayed_deleted_tetrahedra > stat_delayed_deleted_tetrahedra_max) stat_delayed_deleted_tetrahedra_max = stat_delayed_deleted_tetrahedra;
#endif

  assert(newc>=0 && !(tpool[D3_TETRA(newc)].RAD2_DEAD == RAD2_DEAD_TRUE));
	t_start = D3_TETRA(newc);

#ifdef AUDIT_INSERT
	audit();
#endif
}


//Given negatively oriented tetra [V] such that V[pi] and q are on the opposite
//side of the plane through V\V[pi]. 
//This tests whether q is inside the cone formed by shooting rays from V[pi]
//through the triangle V\V[pi]. 
//Another way to say is testing whether the union of tetras [V] and [V\V[pi]\cup{q}] 
//is convex
inline bool IsInCone(Delaunay3Vertex* V[], int pi, Delaunay3Vertex* q)
{
	Delaunay3Vertex* temp; int qi;
	assert(pi>=0 && pi<4);
	qi = pi + offset[1][pi]; 
	SWAP(V[qi],q,temp); 
	if (volume(V[0], V[1], V[2], V[3])>0) return false;
	SWAP(V[qi],q,temp); 

	qi = pi + offset[2][pi]; 
	SWAP(V[qi],q,temp); 
	if (volume(V[0], V[1], V[2], V[3])>0) return false;
	SWAP(V[qi],q,temp); 
 
	qi = pi + offset[3][pi]; 
	SWAP(V[qi],q,temp); 
	if (volume(V[0], V[1], V[2], V[3])>0) return false;
	SWAP(V[qi],q,temp); 
 
	return true;
}

//given a tetra [V[0...3]] indexed by D3_TETRA(c). and a line qp that satisfies the following
//conditions. Return whether pq intersects the triangle t opposite c
// - qp intersects [V]
// - qp does not enter [V] through t 
//for the inExact Delaunay, we can get a little speed up from the previously 
//computed I and sphere equation stored
bool Delaunay3::intersectTriangle(int c, Delaunay3Vertex* q, Delaunay3Vertex* p)
{
	//first check whether pq intersects the plane through the triangle
	Delaunay3Vertex* tV[4]; 
	int i = D3_INDEX(c);
 
	int t = D3_TETRA(c);
	tV[0] = tpool[t].V[0];
	tV[1] = tpool[t].V[1];
	tV[2] = tpool[t].V[2];
	tV[3] = tpool[t].V[3];

	tV[i] = p;
 
	if (volume(tV[0],tV[1],tV[2],tV[3])>0) 	return false;

	//check whether q is on the other side of the plane. This should be true given the precondition
#ifndef NDEBUG 
	Delaunay3Vertex* temp; SWAP(tV[i],q,temp);
	assert(volume(tV[0],tV[1],tV[2],tV[3])>=0); 
	SWAP(tV[i],q,temp); 
#endif

 	// now check if the line intersects in the interior of the triangle
  return(IsInCone(tV,i,q));  
}

bool Delaunay3::DuplicateVertex(Delaunay3Tetrahedron& t, Delaunay3Vertex* p)
{
  return ( ( t.V[0]->x[0] == p->x[0] && t.V[0]->x[1] == p->x[1]  && t.V[0]->x[2] == p->x[2] && t.V[0] != pinf ) ||
					 ( t.V[1]->x[0] == p->x[0] && t.V[1]->x[1] == p->x[1]  && t.V[1]->x[2] == p->x[2] && t.V[1] != pinf ) ||
					 ( t.V[2]->x[0] == p->x[0] && t.V[2]->x[1] == p->x[1]  && t.V[2]->x[2] == p->x[2] ) ||
					 ( t.V[3]->x[0] == p->x[0] && t.V[3]->x[1] == p->x[1]  && t.V[3]->x[2] == p->x[2] ) );
}

int Delaunay3::locate(Delaunay3Vertex* p)
{
#ifdef COLLECT_STATISTICS
  stat_locate_calls++;
#endif

  // find a triangle face f opposite c0 the vertex at c0 and p are on the same 
  // side of f and set q to be in the interior of f
	int t = t_start;
  
  // t_start could be deleted or finalized
  if (tpool[t].RAD2_DEAD == RAD2_DEAD_TRUE)
  {
#ifdef COLLECT_STATISTICS
    stat_t_start_was_dead++;
#endif
    // use the finalization grid to locate a good starting tetrahedra
    t = otherStartForLocate(p);
    // we did not find anything 
    if (t == -1)
    {
#ifdef COLLECT_STATISTICS
      stat_t_start_no_other++;
#endif
      // so look for the last live triangle (not necessarily the last created)
      // in the tpool as the starting triangle,
      for (t = tetrahedron_buffer_maxsize-1; t>=0; t--)
      {
        if (tpool[t].RAD2_DEAD != RAD2_DEAD_TRUE)
        {
          break;
        }
      }
    }
  }

  assert(t>=0);

#ifdef COLLECT_STATISTICS
  stat_locate_active_checked++;
  stat_locate_active_run = 0;
#endif

  // check the circumsphere of the first tet we picked 
  double I = inSphereExact(t,p);

  // does it already contain the point
	if (I < 0)
  {
#ifdef COLLECT_STATISTICS
    stat_locate_after_0walk++;
#endif
    return t;
  }
	else if (I == 0)
  {
    if (DuplicateVertex(tpool[t],p))
    {
      return DUPLICATE_POINT;
    }
  }

  int c = -1;

	if (tpool[t].V[0] == pinf || tpool[t].V[1] == pinf) // if v0 or v1 are infinite
  {
    c = tpool[t].N[(tpool[t].V[0] == pinf ? 0 : 1)];
		if (c == D3_NULL_INDEX)
    {
#ifdef COLLECT_STATISTICS
      stat_locate_failed_notsurewhy++;
#endif
			return bruteLocate(p);
    }
    I = inSphereExact(D3_TETRA(c),p);
		if (I < 0)
    {
#ifdef COLLECT_STATISTICS
      stat_locate_after_1walk++;
#endif
			return D3_TETRA(c);
    }
		else if (I == 0 && DuplicateVertex(tpool[D3_TETRA(c)],p))  
    {
			return DUPLICATE_POINT; 
    }
	}
  else
  {
		if      ( volume(p, tpool[t].V[1], tpool[t].V[2], tpool[t].V[3])>0) c = D3_CORNER(t,0);
		else if ( volume(tpool[t].V[0], p, tpool[t].V[2], tpool[t].V[3])>0) c = D3_CORNER(t,1);
		else if ( volume(tpool[t].V[0], tpool[t].V[1], p, tpool[t].V[3])>0) c = D3_CORNER(t,2);
		else if ( volume(tpool[t].V[0], tpool[t].V[1], tpool[t].V[2], p)>0) c = D3_CORNER(t,3);
	}

  assert(c != -1);

  t = D3_TETRA(c);

  assert(tpool[t].RAD2_DEAD != RAD2_DEAD_TRUE);
  assert(tpool[t].V[0]!=NULL && tpool[t].V[1]!=NULL && tpool[t].V[2]!=NULL && tpool[t].V[3]!=NULL);

  //compute q, the other end of the striaght line we walk.
	Delaunay3Vertex q;  
	int c_index = D3_INDEX(c);
  Average(&q, tpool[t].V[c_index+offset[1][c_index]],
					    tpool[t].V[c_index+offset[2][c_index]],
					    tpool[t].V[c_index+offset[3][c_index]]);

  //walk and maintain the    
#define MAX_LOCATE_STEPS 3000
 	int steps = 0; 
  while((steps++)<MAX_LOCATE_STEPS)
  { 
    //check the three triangle sides for exit
    int k;
    for (k=1; k<=3; k++)
    {
      int c1 = c + offset[k][D3_INDEX(c)];
      if (intersectTriangle(c1, &q, p))
      { 
        if (tpool[D3_TETRA(c1)].N[D3_INDEX(c1)]==D3_NULL_INDEX)
        {
#ifdef COLLECT_STATISTICS
          stat_locate_failed_deleted++;
#endif
          return bruteLocate(p);
        }
        c = tpool[D3_TETRA(c1)].N[D3_INDEX(c1)]; 
				t = D3_TETRA(c);
        break;
      }
    }

    if (k==4)
    {
#ifdef COLLECT_STATISTICS
      stat_locate_failed_noexit++;
#endif
      return bruteLocate(p); //can't find the exit triangle. This could happen because ... why?
    }

#ifdef COLLECT_STATISTICS
  stat_locate_active_checked++;
  stat_locate_active_run++;
#endif

    I = inSphereExact(t,p);

    if (I < 0)
    {
#ifdef COLLECT_STATISTICS
      if (stat_locate_active_run > stat_locate_active_maxrun) stat_locate_active_maxrun = stat_locate_active_run;
      if (stat_locate_active_run == 1) stat_locate_after_1walk++;
#endif
      return t;
    }
    else if (I == 0)
    {
      if (DuplicateVertex(tpool[D3_TETRA(c)],p))
      {
        return DUPLICATE_POINT;
      }
      if (isInf(t))
      {
#ifdef COLLECT_STATISTICS
        stat_locate_failed_infinite++;
#endif
        return bruteLocate(p);
      }
    }
  }
  
#ifdef COLLECT_STATISTICS
  stat_locate_failed_steps++;
#endif

  return bruteLocate(p);  
}

extern int getDataFromLeaf(const float* pos_f);

int Delaunay3::otherStartForLocate(Delaunay3Vertex* v)
{
  float pos_f[3];
  pos_f[0] = (float)v->x[0];
  pos_f[1] = (float)v->x[1];
  pos_f[2] = (float)v->x[2];
  int idx = getDataFromLeaf(pos_f);
  while (idx != -1)
  {
    if (tpool[idx].RAD2_DEAD == RAD2_DEAD_TRUE || idx == t_start)
    {
      idx = tpool[idx].next;
    }
    else
    {
      return idx;
    }
  }
  return -1;
}

int Delaunay3::bruteLocate(Delaunay3Vertex* v)
{
  if (secondTry == true)
  {
    int t = otherStartForLocate(v);
    if (t != -1)
    {
      t_start = t;
      secondTry = false;
      t = locate(v);
      secondTry = true;
      return t;
    }
  }

#ifdef COLLECT_STATISTICS
	stat_locateold_calls++;
  stat_locateold_active_total += tetrahedron_buffer_maxsize;
#endif

  double d;	
  for (int i=0; i<tetrahedron_buffer_maxsize; i++)
  {
    if (tpool[i].RAD2_DEAD == RAD2_DEAD_TRUE) continue; 

#ifdef COLLECT_STATISTICS
    stat_locateold_active_checked++;
#endif

    if ((d=inSphereExact(i,v))<0)
    {
      return i;
    }
		else if (d==0 && DuplicateVertex(tpool[i],v))
    {
	    return DUPLICATE_POINT; 
    }
 	}
	return LOCATE_FAIL;
}

bool Delaunay3::isInf(int i)
{
 	assert(tpool[i].V[2]!=pinf && tpool[i].V[3]!=pinf);
	return(tpool[i].V[0]==pinf || tpool[i].V[1]==pinf);
}

void Delaunay3::cornerPrint(int c)
{ 
	Delaunay3Tetrahedron& tet = tpool[D3_TETRA(c)];
	int c_index = D3_INDEX(c);
	
	fprintf(stderr,"%%%3d(%2d,%1d)%2d=", c, D3_TETRA(c), c_index, tet.V[D3_INDEX(c)]);

	fprintf(stderr,"(%d %5.0f %5.0f %5.0f) opp:%4d(%3d,%2d) \n", 
	   tet.V[c_index]==pinf?0:1, 
		 tet.V[c_index]->x[0], 
		 tet.V[c_index]->x[1], 
		 tet.V[c_index]->x[2], 
		 tet.N[c_index],
		 D3_TETRA(tet.N[c_index]),D3_INDEX(tet.N[c_index])); 
	(void)fflush(stdout);
}

void Delaunay3::cornerPrint4(int c)
{ 
	int t =D3_TETRA(c);
	SSsphere* sp = &tpool[t];

	fprintf(stderr,"disp('Sphere(%d) = <%5.0f %5.0f %5.0f %5.0f>')\n", t, 
				 sp->cen[0],sp->cen[1],sp->cen[2],sp->rad);

	fprintf(stderr,"DetCheckH([");  
	for (int k = 0; k < 4; k++) 
		fprintf(stderr," %d %5.0f %5.0f %5.0f %5.0f;\n", 
					 (tpool[t].V[k]==pinf?0:1), tpool[t].V[k]->x[0], tpool[t].V[k]->x[1],
					                 tpool[t].V[k]->x[2],tpool[t].V[k]->sq);
	fprintf(stderr,"]);\n");
}

void Delaunay3::audit()
{
	int p,b,i,k;
	double d; 
	  
	for (p = 0; p < tetrahedron_buffer_maxsize; p++) { 
		if (tpool[p].RAD2_DEAD == RAD2_DEAD_TRUE) continue; // don't audit tetras on free list
		b = D3_CORNER(p,0);
		for (int j=0; j<4; j++) { // per corner checks
			int c = D3_CORNER(p,j);
			i = tpool[p].N[j]; // check opposite
			if (i==D3_NULL_INDEX) continue;
			if (tpool[D3_TETRA(i)].N[D3_INDEX(i)] != c) {
				fprintf(stderr,"%%AUDIT: wrong opp.opp \n");
				cornerPrint(c); 
				cornerPrint(i); 
			}
			 
			// check sphere opposite corner c
			k = tpool[p].N[j];
			int t1 = D3_TETRA(k);
			if (tpool[D3_TETRA(c)].V[D3_INDEX(c)]==pinf) {
				d = volume(tpool[t1].V[0], tpool[t1].V[1],tpool[t1].V[2], tpool[t1].V[3]);
			} else {
 				d = inSphereExact(t1, tpool[D3_TETRA(c)].V[D3_INDEX(c)]);
			}
   
			if (d  < 0) { 
				fprintf(stderr,"d <0\n");
			} 
		}
	}
}
