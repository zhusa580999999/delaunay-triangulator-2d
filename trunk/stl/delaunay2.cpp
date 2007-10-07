/*
===============================================================================

  FILE:  delaunay2.cpp
  
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
#include <iostream.h>
#include <algorithm>
#include "delaunay2.h"

#define EXACTPREDICATE
//#undef EXACTPREDICATE

using namespace std;

#define NEXT(i) ((i+1)%3)
#define PREV(i) ((i+2)%3)

#define CTRIP(i) (triangle_buffer+CTRI(i))
#define OPP(c) (CTRIP(c)->N[CIND(c)])
#define MAKEOPP(c, c1) {CTRIP(c)->N[CIND(c)]=c1;} //make c1 the opposite of c
#define COMPLEMENT(i) (-(i+1))

#define DUPLICATE_POINT -1
#define LOCATE_FAIL -2

static inline int SignDet33(double* a, double *b, double* c)
{
	double d = 
         (-a[2]*b[1]*c[0] + a[1]* b[2]* c[0] 
         + a[2]*b[0]*c[1] - a[0]* b[2]* c[1] 
				 - a[1]*b[0]*c[2] + a[0]* b[1]* c[2]);

	if (d==0) return 0;
	return (d<0? -1: 1);
}

//Assuming p is on the line ab, returns -1, 0, 1, if p is inside,
//on, or outside the segment, respectively
#define SIGN(a) ((a)==0? 0 : ((a)<0? -1:1))
static inline int InSegment(Delaunay2Vertex* a, Delaunay2Vertex* b, Delaunay2Vertex* p)
{
	assert(!(a->x[0]==b->x[0]&&a->x[1]==b->x[1]));
	if (a->x[0]!=b->x[0]){
		return (SIGN( (a->x[0]-p->x[0])*(b->x[0]-p->x[0])));
 	}
	else{
		return (SIGN( (a->x[1]-p->x[1])*(b->x[1]-p->x[1])));
	}
}

/* statistics */
#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

static int stat_duplicate_points = 0;

#ifdef COLLECT_STATISTICS

static __int64 stat_new_triangles = 0;
static __int64 stat_deleted_triangles = 0;

static int stat_locate_calls = 0;
static int stat_locate_after_0walk = 0;
static int stat_locate_after_1walk = 0;
static int stat_locate_active_maxrun = 0;
static __int64 stat_locate_active_checked = 0;
static int stat_locate_failed_steps = 0;
static int stat_locate_failed_deleted = 0;

static int stat_locateold_calls = 0;
static __int64 stat_locateold_active_total = 0;
static __int64 stat_locateold_active_checked = 0;
#endif

//Jonathan Shewchuk's code 
void exactinit();
double incircle(double *pa, double* pb, double* pc, double* pd);
double orient2d(double *pa, double* pb, double* pc);

int area_sign(Delaunay2Vertex* v1,Delaunay2Vertex* v2, Delaunay2Vertex* v3 )
{
#ifdef EXACTPREDICATE
	double d=orient2d(v1->x, v2->x, v3->x);
	if (d==0) return 0;
	else return (d<0? -1: 1);
#else
	double a[3], b[3], c[3];
	a[0] = 1; a[1] = v1->x[0]; a[2] = v1->x[1];
	b[0] = 1; b[1] = v2->x[0]; b[2] = v2->x[1];
	c[0] = 1; c[1] = v3->x[0]; c[2] = v3->x[1];
	return SignDet33(a,b,c);
#endif
}

// efficient memory allocation

Delaunay2Vertex* Delaunay2::allocDelaunay2Vertex()
{
  if (vertex_buffer == 0)
  {
    vertex_buffer = (Delaunay2Vertex*)malloc(sizeof(Delaunay2Vertex)*vertex_buffer_alloc);
    if (vertex_buffer == 0)
    {
      fprintf(stderr,"malloc for delaunay vertex buffer failed\n");
      return 0;
    }
    for (int i = 0; i < vertex_buffer_alloc; i++)
    {
      vertex_buffer[i].buffer_next = &(vertex_buffer[i+1]);
    }
    vertex_buffer[vertex_buffer_alloc-1].buffer_next = 0;
    vertex_buffer_alloc = 1.5*vertex_buffer_alloc;
  }
  // get pointer to next available delaunay vertex
  Delaunay2Vertex* vertex = vertex_buffer;
  vertex_buffer = vertex->buffer_next;

  // init delaunay vertex
  vertex->index = -1;
  vertex->use_count = 0;

  vertex_buffer_size++; if (vertex_buffer_size > vertex_buffer_maxsize) vertex_buffer_maxsize = vertex_buffer_size;

  return vertex;
}

Delaunay2Vertex* Delaunay2::allocDelaunay2Vertex(const float* xyh)
{
  Delaunay2Vertex* vertex = allocDelaunay2Vertex();

  vertex->x[0] = xyh[0];
  vertex->x[1] = xyh[1];
  vertex->x[2] = LIFT(xyh[0],xyh[1]);
  vertex->h = xyh[2];

  return vertex;
}

Delaunay2Vertex* Delaunay2::allocDelaunay2Vertex(const double* xyh)
{
  Delaunay2Vertex* vertex = allocDelaunay2Vertex();

  vertex->x[0] = xyh[0];
  vertex->x[1] = xyh[1];
  vertex->x[2] = LIFT(xyh[0],xyh[1]);
  vertex->h = (float)xyh[2];

  return vertex;
}

Delaunay2Vertex* Delaunay2::allocDelaunay2Vertex(double x, double y, double h)
{
  Delaunay2Vertex* vertex = allocDelaunay2Vertex();

  vertex->x[0] = x;
  vertex->x[1] = y;
  vertex->x[2] = LIFT(x,y);
  vertex->h = (float)h;

  return vertex;
}

void Delaunay2::deallocDelaunay2Vertex(Delaunay2Vertex* vertex)
{
  vertex->buffer_next = vertex_buffer;
  vertex_buffer = vertex;
  vertex_buffer_size--;
}

TRI Delaunay2::allocDelaunay2Triangle()
{
  // assure that there is enough memory
  if (triangle_buffer_next == -1)
  {
    if (triangle_buffer)
    {
      triangle_buffer_alloc = 1.5*triangle_buffer_alloc;
      triangle_buffer = (Delaunay2Triangle*)realloc(triangle_buffer,sizeof(Delaunay2Triangle)*triangle_buffer_alloc);
      if (triangle_buffer == 0)
      {
        fprintf(stderr,"realloc %d for delaunay triangle buffer failed\n", triangle_buffer_alloc);
        return 0;
      }
    }
    else
    {
      triangle_buffer = (Delaunay2Triangle*)malloc(sizeof(Delaunay2Triangle)*triangle_buffer_alloc);
      if (triangle_buffer == 0)
      {
        fprintf(stderr,"malloc %d for delaunay triangle buffer failed\n", triangle_buffer_alloc);
        return 0;
      }
    }
    for (int i = triangle_buffer_size; i < triangle_buffer_alloc; i++)
    {
      triangle_buffer[i].buffer_next = i+1;
    }
    triangle_buffer[triangle_buffer_alloc-1].buffer_next = -1;
    triangle_buffer_next = triangle_buffer_size;
  }

  // get next available delaunay triangle
  TRI triangle_idx = triangle_buffer_next;
  Delaunay2Triangle* triangle = &(triangle_buffer[triangle_idx]);
  triangle_buffer_next = triangle->buffer_next;
  triangle_buffer_size++;
  
  if (triangle_buffer_size > triangle_buffer_maxsize) triangle_buffer_maxsize = triangle_buffer_size;

  // init delaunay triangle
  triangle->dead = false;
  triangle->N[0] = triangle->N[1] = triangle->N[2] = 0;
  triangle->rad = -1.0f;

  return triangle_idx;
}

void Delaunay2::deallocDelaunay2Triangle(TRI triangle_idx)
{
  triangle_buffer[triangle_idx].dead = true;
  triangle_buffer[triangle_idx].buffer_next = triangle_buffer_next;
  triangle_buffer_next = triangle_idx;
  triangle_buffer_size--;
}

void Delaunay2::addActiveTriangle(int t_idx)
{
  if (act_t_num == 0)
  {
    triangle_buffer[t_idx].buffer_prev = t_idx;
    triangle_buffer[t_idx].buffer_next = t_idx;
    act_t_list = t_idx;
  }
  else
  {
    Delaunay2Triangle* list = triangle_buffer + act_t_list;
    triangle_buffer[list->buffer_prev].buffer_next = t_idx;
    triangle_buffer[t_idx].buffer_prev = list->buffer_prev;
    triangle_buffer[t_idx].buffer_next = act_t_list;
    list->buffer_prev = t_idx;
  }
  act_t_num++;
}

void Delaunay2::delActiveTriangle(int t_idx)
{
  if (act_t_num > 1)
  {
    Delaunay2Triangle* t = triangle_buffer + t_idx;
    triangle_buffer[t->buffer_prev].buffer_next = t->buffer_next;
    triangle_buffer[t->buffer_next].buffer_prev = t->buffer_prev;
    if (t_idx == act_t_list)
    {
      act_t_list = t->buffer_next;
    }
  }
  act_t_num--;
}

int Delaunay2::nextYoungestActiveTriangle(int t_idx) const
{
  return triangle_buffer[t_idx].buffer_prev;
}

int Delaunay2::youngestActiveTriangle() const
{
  return triangle_buffer[act_t_list].buffer_prev;
}

// a = (b1+b2)/2
void Delaunay2::Average(Delaunay2Vertex* a, Delaunay2Vertex* b1, Delaunay2Vertex* b2)
{
  a->x[0] = (b1->x[0] + b2->x[0])/2.0;
  a->x[1] = (b1->x[1] + b2->x[1])/2.0;
  a->x[2] = LIFT(a->x[0], a->x[1]);
}

//does p have the same coordinates as a vertex in t?
bool Delaunay2::DuplicateVertex(Delaunay2Triangle* t, Delaunay2Vertex* p)
{
  return ( (t->V[0]!=pinf && t->V[0]->x[0]== p->x[0] && t->V[0]->x[1]== p->x[1]  ) ||
           (t->V[1]!=pinf && t->V[1]->x[0]== p->x[0] && t->V[1]->x[1]== p->x[1]  ) || 
	         (t->V[2]!=pinf && t->V[2]->x[0]== p->x[0] && t->V[2]->x[1]== p->x[1]  ) );
}

bool Delaunay2::mutualNeighbor(Delaunay2Triangle* t, int ci)
{
	CNR c = t->N[ci];
	return (CTRIP(c)->N[CIND(c)] == CORNER(t-triangle_buffer,ci));
}

void Delaunay2::audit()
{
//  printf("WARNING ... audit is enabled\n");
    return;

  for (int i=0; i<triangle_buffer_maxsize; i++)
  {
		Delaunay2Triangle* t = triangle_buffer+i;

		if (t->dead)
    {
      continue;
    }

    assert( t->V[0]->use_count > 0 );
    assert( t->V[1]->use_count > 0 );
    assert( t->V[2]->use_count > 0 );

    if (t->V[0] != pinf && t->V[1] != pinf)
    {
      assert( area_sign(t->V[0],t->V[1],t->V[2]) > 0 );
    }

		// commented out because deletion from streaming
		// 		assert(t->V[2]!=pinf);
		for (int j=0; j<3; j++)
    {
			if (t->N[j]!=D2_NULL_INDEX) assert(mutualNeighbor(t, j));
		}
	}
}

void Delaunay2::stats()
{
	int infin_count=0;
	int finit_count=0;

  for (int i=0; i<triangle_buffer_maxsize; i++)
  {
		Delaunay2Triangle* t = &(triangle_buffer[i]);

    if (t->dead)
    {
      continue;
    }

		if (t->V[0]==pinf || t->V[1]==pinf || t->V[2]==pinf)
    {
			infin_count++;
		}
		else
    {
			finit_count++;
      printf("(%g %g) / (%g %g) / (%g %g) \n", t->V[0]->x[0], t->V[0]->x[1], t->V[1]->x[0], t->V[1]->x[1], t->V[2]->x[0], t->V[2]->x[1]);
    }
	}
	printf("There are %3d    finite triangles\n", finit_count);
	printf("There are %3d  infinite triangles\n", infin_count);
}

// geometric predicates

TRI Delaunay2::locate_brute(Delaunay2Vertex* p)
{
  int sign;
  Delaunay2Triangle* t;
  int t_idx = youngestActiveTriangle();

#ifdef COLLECT_STATISTICS
  stat_locateold_calls++;
  stat_locateold_active_total += act_t_num;
#endif

  for (int i = 0; i < act_t_num; i++)
  {
#ifdef COLLECT_STATISTICS
    stat_locateold_active_checked++;
#endif

		t = triangle_buffer + t_idx; 

		assert(t->V[0] != NULL);
    assert(t->V[1] != NULL);
    assert(t->V[2] != NULL);

    sign = inSphereExact(t,p);
    if ( sign < 0 )
    {
      return t_idx;
    }
    else if (sign == 0 && DuplicateVertex(t,p))
    {
      return DUPLICATE_POINT;
    }
    t_idx = nextYoungestActiveTriangle(t_idx);
	}
  // we should never get here
//  assert(false);
  fprintf(stderr, "ERROR: fail to locate point\n");
	return LOCATE_FAIL;
}

int Delaunay2::inSphereExact(Delaunay2Triangle* t, Delaunay2Vertex* p)
{
  int sign;
  int inf=D2_NULL_INDEX;

  assert(t->V[2]!=pinf);

#ifdef EXACTPREDICATE

  double d;
  if (t->V[0]==pinf) inf = 0;
  if (t->V[1]==pinf) inf = 1;
  if (inf == D2_NULL_INDEX) // finite circle
  {
    d = -incircle(t->V[0]->x, t->V[1]->x, t->V[2]->x, p->x);
  }
  else // infinite circle
  {
    d = -orient2d(t->V[NEXT(inf)]->x, t->V[PREV(inf)]->x, p->x);
  }	
  if (d==0) sign = 0;
  else sign = d < 0? -1 :1;	

#else

  double a[3], b[3], c[3];
  a[0] = t->V[0]->x[0];  b[0] = t->V[1]->x[0];  c[0] = p->x[0];
  a[1] = t->V[0]->x[1];  b[1] = t->V[1]->x[1];	c[1] = p->x[1];
  a[2] = t->V[0]->x[2];	 b[2] = t->V[1]->x[2];	c[2] = p->x[2];
  if (t->V[0]!=pinf)
  {
    a[0] -= t->V[2]->x[0];
    a[1] -= t->V[2]->x[1];
    a[2] -= t->V[2]->x[2];
  }
  if (t->V[1]!=pinf)
  {
    b[0] -= t->V[2]->x[0];
    b[1] -= t->V[2]->x[1];
    b[2] -= t->V[2]->x[2];
  }
  if (p!=pinf)
  {
    c[0] -= t->V[2]->x[0];
    c[1] -= t->V[2]->x[1];
    c[2] -= t->V[2]->x[2];
  }	
  sign = SignDet33(a,b,c);

#endif

  return sign;
}

bool Delaunay2::inSphere(Delaunay2Triangle* t, Delaunay2Vertex* p)
{
  int sign;
  int inf=D2_NULL_INDEX;

  assert(t->V[2]!=pinf);
  assert((t->V[0]!=NULL && t->V[1]!=NULL && t->V[2]!=NULL));

#ifdef EXACTPREDICATE

  double d;
  if (t->V[0]==pinf) inf = 0;
  if (t->V[1]==pinf) inf = 1;
  if (inf == D2_NULL_INDEX) //finite circle
  { 
    d = -incircle(t->V[0]->x, t->V[1]->x, t->V[2]->x, p->x);
  }
  else //infinite circle
  {
    d = -orient2d(t->V[NEXT(inf)]->x, t->V[PREV(inf)]->x, p->x);
  }
	
  if (d==0) sign=0;
  else sign = d<0? -1 :1;	

#else

  double a[3], b[3], c[3];

  a[0] = t->V[0]->x[0];  b[0] = t->V[1]->x[0];  c[0] = p->x[0];
  a[1] = t->V[0]->x[1];  b[1] = t->V[1]->x[1];	c[1] = p->x[1];
  a[2] = t->V[0]->x[2];	 b[2] = t->V[1]->x[2];	c[2] = p->x[2];
		
  if (t->V[0]!=pinf)
  {
    a[0] -= t->V[2]->x[0];
    a[1] -= t->V[2]->x[1];
    a[2] -= t->V[2]->x[2];
  }
  else
  {
    inf = 0;
  }
  if (t->V[1]!=pinf)
  {
    b[0] -= t->V[2]->x[0];
    b[1] -= t->V[2]->x[1];
    b[2] -= t->V[2]->x[2];
  }
  else
  {
    inf = 1;
  }
  if (p!=pinf)
  {
    c[0] -= t->V[2]->x[0];
    c[1] -= t->V[2]->x[1];
    c[2] -= t->V[2]->x[2];
  }
		
  sign = SignDet33(a,b,c);
#endif

	if (sign == 0)
  {
	  if (inf != D2_NULL_INDEX)
    {
		  //check if p is in the 1-dimensional sphere (a segment)
		  //   t->V[NEXT(inf)]  t->V[PREV(inf)
		  //since the algebraic degree of InSegment is 1
		  //we do not need special code for numerical computation
 		  return (InSegment(t->V[NEXT(inf)], t->V[PREV(inf)], p)<=0);
    }
    else
    {
      return true;
    }
  }
	else
  {
		return (sign <= 0);
	}
}

#define MAX_LOCATE_STEPS 3000

TRI Delaunay2::locate(Delaunay2Vertex* p)
{
#ifdef COLLECT_STATISTICS
  stat_locate_calls++;
#endif

  //find an edge (t, ci) so that the triangle and p are on the same 
  //side of it and set q to be the mid point of the edge
  int t_idx = youngestActiveTriangle();
  Delaunay2Triangle* t = triangle_buffer + t_idx;
  
#ifdef COLLECT_STATISTICS
  stat_locate_active_checked++;
  int stat_locate_active_run = 1;
#endif

  assert(!t->dead && t->V[0]!=NULL && t->V[1]!=NULL && t->V[2]!=NULL);

  int sign = inSphereExact(t,p);

  if (sign < 0)
  {
#ifdef COLLECT_STATISTICS
    stat_locate_after_0walk++;
#endif
    return t_idx;
  }
  else if (sign == 0 && DuplicateVertex(t,p))
  {
    return DUPLICATE_POINT;
  }

  int ci = -1;
  if      (area_sign(p, t->V[1], t->V[2]) > 0) ci = 0; 
  else if (area_sign(t->V[0], p, t->V[2]) > 0) ci = 1; 
  else if (area_sign(t->V[0], t->V[1], p) > 0) ci = 2; 

  if (ci == -1) // must be an infinite triangle, use the finite edge as the starting edge
  {
    assert(IsInf(t));
    for (int k=0; k<3; k++)
    {
      if (t->V[k]==pinf)
      {
        if (t->N[k] == D2_NULL_INDEX)
        {
          return locate_brute(p);
        }
        ci = CIND(t->N[k]);
        t_idx = CTRI(t->N[k]);
        break;
      }
    }

    t = triangle_buffer + t_idx;

#ifdef COLLECT_STATISTICS
    stat_locate_active_checked++;
    stat_locate_active_run++;
#endif

    assert(!t->dead && t->V[0]!=NULL && t->V[1]!=NULL && t->V[2]!=NULL);

    sign = inSphereExact(t,p);

    if (sign < 0)
    {
      return t_idx;
    }
    else if (sign==0 && DuplicateVertex(t,p)) 
    {
      return DUPLICATE_POINT;
    }
  }

  //compute q, the other end of the striaght line we walk.
  Delaunay2Vertex q;
  Average(&q, t->V[PREV(ci)], t->V[NEXT(ci)]);

  //walk and maintain the edge (t, ci)
  int steps = 0;
  while( (steps++) < MAX_LOCATE_STEPS )
  {
    CNR c1 = t->N[NEXT(ci)];
    CNR c2 = t->N[PREV(ci)];

    if ( area_sign(p, &q, t->V[ci]) <=0 )
    {
      if (c1==D2_NULL_INDEX) break;
      t = triangle_buffer + CTRI(c1);
      ci = CIND(c1);
    }
    else
    {
      if (c2==D2_NULL_INDEX) break;
      t = triangle_buffer + CTRI(c2);
      ci = CIND(c2);
    }     
    
    assert(!t->dead && t->V[0]!=NULL && t->V[1]!=NULL && t->V[2]!=NULL);

#ifdef COLLECT_STATISTICS
    stat_locate_active_checked++;
    stat_locate_active_run++;
#endif

    sign = inSphereExact(t,p);

    if (sign < 0)
    {
#ifdef COLLECT_STATISTICS
      if (stat_locate_active_run > stat_locate_active_maxrun) stat_locate_active_maxrun = stat_locate_active_run;
#endif
      return (t-triangle_buffer);
    }
    else if (sign==0 && DuplicateVertex(t,p)) 
    {
#ifdef COLLECT_STATISTICS
      if (stat_locate_active_run > stat_locate_active_maxrun) stat_locate_active_maxrun = stat_locate_active_run;
#endif
      return DUPLICATE_POINT;
    }
  }

//  printf("brute \n");
	//the walk has failed either because
	//it runs into some deleted region
	//or becasue it takes too many steps, indicating 
	//there is an infinite loop.(there shouldn't be, but let's just be safe)
	//so just use bruteLocate

#ifdef COLLECT_STATISTICS
  if (steps < MAX_LOCATE_STEPS)
    stat_locate_failed_deleted++;
  else
    stat_locate_failed_steps++;
#endif

  return locate_brute(p);
}

//after search
#define POP(stack) (stack##st[stack##sp--])
#define isEMPTY(stack) (stack##sp < 0)
#define stkINIT(stack) {stack##sp = -1; }
#define PUSH(value, stack) { \
	stack##st[++stack##sp] = value; \
	if (stack##max < stack##sp) { stack##max = stack##sp; \
	if (stack##max >= MAXSTACK) { \
	printf("ERROR: overflow stack %x pushing %d",  stack##st, value); exit(EXIT_FAILURE); } } \
}
#define stkDECLARE(stack,stn) static int stack##sp, stack##st[MAXSTACK]; static int stack##max;
 
stkDECLARE(dfs,"dfs")  // DFS stack

#define DEC_REF_COUNT(t) {(t)->V[0]->use_count--;(t)->V[1]->use_count--;(t)->V[2]->use_count--;}
#define INC_REF_COUNT(t) {(t)->V[0]->use_count++;(t)->V[1]->use_count++;(t)->V[2]->use_count++;}

void Delaunay2::search(Delaunay2Triangle* t0, Delaunay2Vertex* p)
{
	stkINIT(dfs);

	TRI t0_idx = t0-triangle_buffer;
	t0->dead = true;
	DEC_REF_COUNT(t0);
	dt[ndt++] = t0_idx;

#ifndef NDEBUG
	if (ndt==MAXSTACK)
  {
		fprintf(stderr, "ERROR: ndt stack overflow\n");
    exit(0);
	}
#endif

  for (int i=0; i<3;i++)
  {
		if (t0->N[i]!=D2_NULL_INDEX){
			PUSH(t0->N[i],dfs);
		}
		else{
			PUSH(COMPLEMENT(CORNER(t0_idx,i)),dfs);
		}
	}

	while(!isEMPTY(dfs)){
		CNR c = POP(dfs);
		Delaunay2Triangle*  t; TRI t_idx;	int  ci; 

		bool reachedHorizon;

		//initialize the live triangle given by t, t_idx, ci 
		if (c<0){
			reachedHorizon= true;
		}
		else{
			t_idx = CTRI(c); 
			t = triangle_buffer + t_idx;
			ci = CIND(c);
			if (t->dead) return; //already visited			
			reachedHorizon = !inSphere(t,p);
		}

		if (!reachedHorizon)
    {
			t->dead = true;
			dt[ndt] = t_idx;
			DEC_REF_COUNT(t);
			ndt++;

#ifndef NDEBUG
	    if (ndt==MAXSTACK)
      {
		    fprintf(stderr, "ERROR: ndt stack overflow\n");
        exit(0);
	    }
#endif

			for (int i=0; i<2;i++)
      {
				t = triangle_buffer + t_idx;
				ci = NEXT(ci);   
				if (t->N[ci]!=D2_NULL_INDEX){
					PUSH(t->N[ci],dfs);
				}
				else{ 
					PUSH(COMPLEMENT(CORNER(t_idx,ci)),dfs);//push the negative corner, indicating 
				}
			}
		}
		else
    {
			Delaunay2Triangle *t0; int c0i; TRI t0_idx;
			CNR c0; 
			
			//intialize the dead corne c0
			if (c<0) c0 = COMPLEMENT(c);
			else c0 = t->N[ci];

			//make new triangle
			TRI tp_idx = allocDelaunay2Triangle();
			Delaunay2Triangle* tp = triangle_buffer+tp_idx;

			//retrieve dead triangle
			t0 = CTRIP(c0); c0i = CIND(c0);
			t0_idx = t0-triangle_buffer;


      tp->V[0] = t0->V[0];  // tp is oriented the same way as that now dead triangle t
	    tp->V[1] = t0->V[1];  
	    tp->V[2] = t0->V[2];  
	    tp->V[c0i] = p;       // but the vertex opposite from the edge shared with t is p

			INC_REF_COUNT(tp); 

			// because pointer t could have changed during realloc
			t =  triangle_buffer + t_idx; //make sense only if c<0, but doesn't hurt
			t0 = triangle_buffer + t0_idx;

			CNR cp = CORNER(tp_idx, c0i);

			//set up opp pointers
			if (c<0){
				tp->N[c0i]=D2_NULL_INDEX;
			}
			else{
				tp->N[c0i]=c; 	
				t->N[ci]=cp;		
			}	
			t0->N[c0i] = COMPLEMENT(cp);

			//stack up the new corners
			nc[nnc] = cp;
			dc[nnc] = c0;
			nnc++;
#ifndef NDEBUG
			if (nnc==MAXSTACK)
      {
				fprintf(stderr, "ERROR: nnc stack overflow\n");
        exit(0);
			}
#endif
		}
	}

  for (i=0; i<ndt; i++)
  {
	  delActiveTriangle(dt[i]);
  }

  for (i=0; i<nnc; i++)
  {
	  addActiveTriangle(CTRI(nc[i]));
  }		
}

//Delaunay update functions
void Delaunay2::insert(Delaunay2Vertex* p)
{
#ifdef VERBOSE
	fprintf(stderr, "insert %g %g %g\n", p->x[0], p->x[1], p->h);
#endif

  // find a triangle whose circumcicle contains p

  TRI t0_idx = locate(p);

  // discard duplicate points 

	if (t0_idx == DUPLICATE_POINT)
  {
    deallocDelaunay2Vertex(p);
//#ifdef COLLECT_STATISTICS
    stat_duplicate_points++;
//#endif
    return;
	}

  nnc = ndt = 0;

	search(triangle_buffer+t0_idx, p); 

  assert(nnc>0);

  int i;

	//connect neighbors among new triangles
  for (i = 0; i < nnc; i++)
  {
		Delaunay2Triangle* tp = CTRIP(nc[i]); // the new triangle 
	  int cp = CIND(nc[i]); //the index of the new corner

		TRI tN; int cN; //neighbor 
 
		//connect tp to its two new neighbors by traversing
		//around around its two horizon vertices and maintains 
    //c and c>=0.  If c<0, it indicates that -c is the new corner 
		//the new corner facing the same horizon edge 
		CNR c;

		c = dc[i];
		while(c >= 0)
    {
			Delaunay2Triangle* s = CTRIP(c);  int si = CIND(c);
			c = s->N[NEXT(si)];
		};
 		assert(c<0);
		tN = CTRI(COMPLEMENT(c)); cN = CIND(COMPLEMENT(c));   
  	tp->N[NEXT(cp)] =  CORNER(tN, PREV(cN)) ;  //pointing cp down, the left neighbor

		c = dc[i];
		while(c >= 0)
    {
			Delaunay2Triangle* s = CTRIP(c);  int si = CIND(c);
			c = s->N[PREV(si)];
		};
 		assert(c<0);
		tN = CTRI(COMPLEMENT(c)); cN = CIND(COMPLEMENT(c));   
  	tp->N[PREV(cp)] =  CORNER(tN, NEXT(cN)) ;  //pointing cp down, the left neighbor    
  }

#ifdef COLLECT_STATISTICS
  stat_new_triangles += nnc;
#endif

  // dealloc the triangles that were removed because of the inserted point
	for (i=0; i<ndt; i++)
  {
		deallocDelaunay2Triangle(dt[i]);
	}

#ifdef COLLECT_STATISTICS
  stat_deleted_triangles += ndt;
#endif

#ifndef NDEBUG		
	audit();
#endif
}

int Delaunay2::SearchCorner(Delaunay2Triangle* t, Delaunay2Vertex* v)
{
	if (t->V[0]==v) return 0;
	if (t->V[1]==v) return 1;
	if (t->V[2]==v) return 2;
	return D2_NULL_INDEX;
}

void Delaunay2::ConnectNeighbors(TRI t1, TRI t2)
{
  Delaunay2Triangle* pt1 = triangle_buffer + t1;
  Delaunay2Triangle* pt2 = triangle_buffer + t2;

  int i;
  int c1=D2_NULL_INDEX;
  int c2=D2_NULL_INDEX;
		
  for (i=0; i<3; i++)
  {
    if (D2_NULL_INDEX==SearchCorner(pt2, pt1->V[i]))	
    {
      assert(c1==D2_NULL_INDEX);
      c1 = i;
    }
  }
  assert(c1!=D2_NULL_INDEX);

  for (i=0; i<3; i++)
  {
    if (D2_NULL_INDEX==SearchCorner(pt1, pt2->V[i]))	
    {
      assert(c2==D2_NULL_INDEX);
      c2 = i;
    }
  }
  assert(c2!=D2_NULL_INDEX);
		
  pt1->N[c1] = CORNER(t2, c2);
  pt2->N[c2] = CORNER(t1, c1);
}
	
void Delaunay2::initialize(Delaunay2Vertex* p1, Delaunay2Vertex* p2, Delaunay2Vertex* p3)
{				
  int i;
  int t_idx[4];
	Delaunay2Triangle* t[4];

  // create the initial four triangles 

  for (i=0; i<4; i++) t_idx[i] = allocDelaunay2Triangle();

  // add them to the linked list of active triangle

  for (i=0; i<4; i++) addActiveTriangle(t_idx[i]);

  // get valid pointers to them

  for (i=0; i<4; i++) t[i] = triangle_buffer + t_idx[i];
		
	// the initial finite triangle ...

  t[0]->V[0]=p1;	 t[0]->V[1]=p2;  t[0]->V[2]=p3;

  // ... should not be degenerate ...

	if (area_sign(t[0]->V[0],t[0]->V[1],t[0]->V[2])==0)
  {
		printf("initial triangle must not be degenerate\n");
		exit(0);
	}

  // ... and be oriented correctly.

  if (area_sign(t[0]->V[0],t[0]->V[1],t[0]->V[2])<0)
  {
		swap(t[0]->V[0],t[0]->V[1]);
  }

	// the other three infinite triangles
		
	for (i=1; i<4; i++)
  {
		t[i]->V[0] = t[0]->V[1];			
		t[i]->V[1] = t[0]->V[0];			
		t[i]->V[2] = t[0]->V[2];			
		t[i]->V[i-1] = pinf;
	}

	// pinf can't be the third vertex so have to fix this

	swap(t[3]->V[2], t[3]->V[1]);
	swap(t[3]->V[1], t[3]->V[0]);

	//set up neighbor pointers

	for (i=0; i<4; i++)
  {
		for (int j=i+1; j<4; j++)
    {
			ConnectNeighbors(t[i]-triangle_buffer,t[j]-triangle_buffer);
		}
	}

  //set the interior point of the triangle
	ct = allocDelaunay2Vertex((p1->x[0]+p2->x[0]+p3->x[0])/3.0, (p1->x[1]+p2->x[1]+p3->x[1])/3.0, 0.0);

  assert( (area_sign(ct, t[0]->V[1], t[0]->V[2]) > 0) &&
          (area_sign(t[0]->V[0], ct, t[0]->V[2]) > 0) &&
      	  (area_sign(t[0]->V[0], t[0]->V[1], ct) > 0) );

  //set up use_count 
	pinf->use_count=3;
	p1->use_count=3;
	p2->use_count=3;
	p3->use_count=3;

#ifndef NDEBUG
	audit();		
#endif		
}

Delaunay2::Delaunay2()
{
	// initialize the vertex pool
	vertex_buffer = 0;
	vertex_buffer_size = 0;
	vertex_buffer_maxsize = 0;
	vertex_buffer_alloc = 1024;

	// initialize the triangle pool
	triangle_buffer = 0;
	triangle_buffer_size = 0;
	triangle_buffer_next = -1;
	triangle_buffer_maxsize = 0;
	triangle_buffer_alloc = 1024;

	// create the point at infinity
	pinf = allocDelaunay2Vertex();
	pinf->x[0]=0;
	pinf->x[1]=0;
	pinf->x[2]=1;

	// make sure that some triangles are allocated
	int idx = allocDelaunay2Triangle();
	deallocDelaunay2Triangle(idx);

	// zero counter and start of active triangle list
	act_t_num = 0;
	act_t_list = 0;

#ifdef EXACTPREDICATE
	fprintf(stderr, "using exact predicates.\n");
	exactinit();
#else
	fprintf(stderr, "using NON-exact predicates.\n");
#endif
}


Delaunay2::~Delaunay2()
{
#ifdef COLLECT_STATISTICS
  char buffer1[128];
  _i64toa(stat_new_triangles, buffer1, 10);
  float buffer2 = (((float)stat_new_triangles/(stat_locate_calls-stat_duplicate_points)*10))/10.0f;
  char buffer3[128];
  _i64toa(stat_deleted_triangles, buffer3, 10);
  float buffer4 = (((float)stat_deleted_triangles/(stat_locate_calls-stat_duplicate_points)*10))/10.0f;
  fprintf(stderr, "stat_new_triangles %s (%.1f) stat_deleted_triangles %s (%.1f)\n", buffer1, buffer2, buffer3, buffer4);
  fprintf(stderr, "stat_locate_calls %d\n", stat_locate_calls);
  _i64toa(stat_locate_active_checked, buffer1, 10);
  buffer2 = (((float)stat_locate_active_checked/stat_locate_calls*10))/10.0f;
  fprintf(stderr, "checked %s (%.1f) walk0 %d walk1 %d maxrun %d\n", buffer1, buffer2, stat_locate_after_0walk, stat_locate_after_1walk,stat_locate_active_maxrun);
  fprintf(stderr, "failed %d (%.1f %%) - steps %d - deleted %d\n", stat_locate_failed_steps+stat_locate_failed_deleted, 100.0f*(stat_locate_failed_steps+stat_locate_failed_deleted)/stat_locate_calls, stat_locate_failed_steps, stat_locate_failed_deleted);

  fprintf(stderr, "stat_locateold_calls %d\n", stat_locateold_calls);
  _i64toa(stat_locateold_active_total, buffer1, 10);
  buffer2 = (((float)stat_locateold_active_total/stat_locateold_calls*10))/10.0f;
  _i64toa(stat_locateold_active_checked, buffer3, 10);
  buffer4 = (((float)stat_locateold_active_checked/stat_locateold_calls*10))/10.0f;
  float buffer5 = (((float)(stat_locateold_active_checked+stat_locate_active_checked)/stat_locate_calls*10))/10.0f;
  fprintf(stderr, "total %s (%.1f) searched %s (%.1f) totavg:(%.1f)\n", buffer1, buffer2, buffer3, buffer4, buffer5);
#endif
  fprintf(stderr, "duplicate_points %d\n", stat_duplicate_points);
}
