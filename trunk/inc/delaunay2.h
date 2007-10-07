/*
===============================================================================

  FILE:  delaunay2.h
  
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

    30 August 2005 - added leo's fix for finding all duplicate points
    21 August 2005 - improved handling of degenerate case for infinite circles.
    19 August 2005 - gotten rid of 'zombie' triangles.
    16 August 2005 - removed unnecessary incircle checks.
    15 August 2005 - improved point location by wwalking at the most recent triangle.
    12 August 2005 - initial prototype completed.  

===============================================================================
*/

#ifndef DELAUNAY2_H
#define DELAUNAY2_H

#include "sscontainer2d.h"

#define TRI int
#define CNR int
#define LIFT(x,y) (x*x+y*y)
#define CTRI(i) ((i)>>2) //the triangle that corner i belongs to
#define CORNER(i, j) (((i)<<2)+j) //the jth corner of triangle i
#define MOD4(a) (a & 3)
#define CIND(i) (MOD4(i))    

#define MAXSTACK 10000
#define D2_NULL_INDEX -19999999

typedef struct Delaunay2Vertex
{
  Delaunay2Vertex* buffer_next; // used for efficient memory management
  double x[3];
  float h;
	int index;
	int use_count;
  bool finalized;
} Delaunay2Vertex;

class Delaunay2Triangle : public SScircle
{
  public:
  int buffer_next; // used for efficient memory management
  int buffer_prev; // used for efficient memory management
	Delaunay2Vertex* V[3];
	CNR N[3];
	bool dead;
};

class Delaunay2
{
public:
	Delaunay2Vertex* ct;
  Delaunay2Vertex* pinf;

  inline bool IsInf(Delaunay2Triangle* t);

	// efficient memory management of vertices

  Delaunay2Vertex* vertex_buffer;
  int vertex_buffer_size;
  int vertex_buffer_maxsize;
  int vertex_buffer_alloc;
	
  Delaunay2Vertex* allocDelaunay2Vertex();
  Delaunay2Vertex* allocDelaunay2Vertex(const float* xyh);
  Delaunay2Vertex* allocDelaunay2Vertex(const double* xyh);
  Delaunay2Vertex* allocDelaunay2Vertex(double x, double y, double h);
  void deallocDelaunay2Vertex(Delaunay2Vertex* v);

	// efficient memory management of triangles

  Delaunay2Triangle* triangle_buffer; 
  TRI triangle_buffer_next;
  int triangle_buffer_size;
  int triangle_buffer_maxsize;
  int triangle_buffer_alloc;

  TRI allocDelaunay2Triangle();
	void deallocDelaunay2Triangle(TRI i);

	// efficient management of active triangles
  int act_t_list;
  int act_t_num;
  void addActiveTriangle(int idx);
  void delActiveTriangle(int idx);
  int youngestActiveTriangle() const;
  int nextYoungestActiveTriangle(int idx) const;

  //global variables for search
	CNR nc[MAXSTACK]; int nnc;
	CNR dc[MAXSTACK];
	TRI dt[MAXSTACK]; int ndt;

	//geometric predicates
  int inSphereExact(Delaunay2Triangle* t, Delaunay2Vertex* p);
	bool inSphere(Delaunay2Triangle* t, Delaunay2Vertex* p);

	//utilities
  void Average(Delaunay2Vertex* a, Delaunay2Vertex* b1, Delaunay2Vertex* b2);
  bool DuplicateVertex(Delaunay2Triangle* t, Delaunay2Vertex* p);
	TRI locate_brute(Delaunay2Vertex* p);
	TRI locate(Delaunay2Vertex* p);
	void search(CNR c, Delaunay2Vertex* p, int root);
  void search(Delaunay2Triangle* t0, Delaunay2Vertex* p);

	int SearchCorner(Delaunay2Triangle* t, Delaunay2Vertex* v);
	void ConnectNeighbors(TRI t1, TRI t2);

	//Delaunay update functions
	void insert(Delaunay2Vertex* p);
	void initialize(Delaunay2Vertex* p1, Delaunay2Vertex* p2, Delaunay2Vertex* p3);

	//audit routines
	bool mutualNeighbor(Delaunay2Triangle* t, int ci);
	void audit();
	void stats();
	void print();

	Delaunay2();
  ~Delaunay2();
};

inline bool Delaunay2::IsInf(Delaunay2Triangle* t)
{
	return (t->V[0]==pinf || t->V[1]==pinf);
}

#endif
