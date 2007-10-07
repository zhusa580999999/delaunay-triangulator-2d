/*
===============================================================================

  FILE:  SScontainer2D.cpp
  
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
#include "sscontainer2d.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

#ifdef COLLECT_STATISTICS
static int stat_finalize_active_justtested = 0;
static int stat_finalize_active_justfinalized = 0;
static int stat_finalize_active_easyintersect = 0;
static int stat_finalize_active_nowfinalized = 0;
static int stat_finalize_active_stillintersect = 0;
#endif

#include <ext/hash_map>

typedef struct SScell
{
  SScell* buffer_next;     // used for efficient memory management
  SScell* parent;
  int idx;
  SScell* child[4];
  int num_children;
  float r_min[2];
  float r_mid[2];
  float r_max[2];
} SScell;

typedef __gnu_cxx::hash_map<int, SScell*> my_cell_hash;

static my_cell_hash* cell_hash = 0;
static my_cell_hash::iterator hash_iterator;

// statistics

static int cell_buffer_size;
static int cell_buffer_maxsize;

// efficient memory management

static int cell_buffer_alloc = 16;
static SScell* cell_buffer_next = 0;

static SScell* allocCell()
{
  if (cell_buffer_next == 0)
  {
    cell_buffer_next = (SScell*)malloc(sizeof(SScell)*cell_buffer_alloc);
    if (cell_buffer_next == 0)
    {
      fprintf(stderr,"ERROR: malloc for cell buffer failed\n");
      return 0;
    }
    for (int i = 0; i < cell_buffer_alloc; i++)
    {
      cell_buffer_next[i].buffer_next = &(cell_buffer_next[i+1]);
    }
    cell_buffer_next[cell_buffer_alloc-1].buffer_next = 0;
    cell_buffer_alloc = 2*cell_buffer_alloc;
  }
  // get index of next available cell
  SScell* cell = cell_buffer_next;
  cell_buffer_next = cell->buffer_next;
 
  // clean cell
  cell->num_children = 0;
  cell->parent = 0;

  cell_buffer_size++; if (cell_buffer_size > cell_buffer_maxsize) cell_buffer_maxsize = cell_buffer_size;

  return cell;
}

static void deallocCell(SScell* cell)
{
  cell->buffer_next = cell_buffer_next;
  cell_buffer_next = cell;
  cell_buffer_size--;
}

static float single_epsilon;
static double double_epsilon;
static float sqrt_epsilon;

void SScircle::initialize(const double* v0, const double* v1, const double* v2)
{
#define LIFT(x,y) (x*x+y*y)
  assert(v0[2] == LIFT(v0[0],v0[1]));
 	assert(v1[2] == LIFT(v1[0],v1[1]));
	assert(v2[2] == LIFT(v2[0],v2[1]));
#undef LIFT

  double x1, y1, q1;
	double x2, y2, q2;

	x1 = v0[0] - v2[0];
	y1 = v0[1] - v2[1];
	q1 = v0[2] - v2[2];

	x2 = v1[0] - v2[0];
	y2 = v1[1] - v2[1];
	q2 = v1[2] - v2[2];

	double Mx, nx;
	double My, ny;
	double Mq, nq ;

#define DET2_AND_ERR(a,b,p,q,i,j) {\
             double t1 = (i##p)*(j##q); \
	           double t2 = (j##p)*(i##q); \
             a = t1 - t2; \
	           b = fabs(t1) + fabs(t2);}
	DET2_AND_ERR(Mx, nx, 1,2,y,q);
	DET2_AND_ERR(My, ny, 1,2,x,q); My*=-1;
	DET2_AND_ERR(Mq, nq, 1,2,x,y);
#undef DET2_AND_ERR

  assert(nx>=0 && ny>=0 && nq>=0);

  // center of circle

  cen[0] = (float)(-Mx/(2.0*Mq));
  cen[1] = (float)(-My/(2.0*Mq));

 	double M1 = -(v2[0]*Mx +v2[1]*My + v2[2]*Mq);

	double Mx2 = Mx*Mx;
	double My2 = My*My;
	double Mq2 = Mq*Mq;

  // squared radius of circle

 	double rad2_tmp = fabs(((Mx2 + My2)/(4*Mq2)) - M1/Mq);

	Mx = fabs(Mx);
	My = fabs(My);
	Mq = fabs(Mq);
	M1 = fabs(M1);

	double topx = Mx*Mq+2*(Mq*nx+Mx*nq);
	double topy = My*Mq+2*(Mq*ny+My*nq);

	double bot = fabs(Mq*(Mq-4*double_epsilon*(Mq+4*nq)));

  if (bot<=0)
  {					
		printf("Zero is divided in the error bound. Use epsilon\n");
		bot = double_epsilon;
	}

  assert(bot<=Mq2); // should be lower bound to maxmize top/bot

	double E1 =  2*(5*M1*Mq + 4*M1*nq)/bot;
	double Ex =  5*topx/bot;
	double Ey =  5*topy/bot;

	float cen_err_xy[2];
  cen_err_xy[0] = (float)(Ex*double_epsilon);
	cen_err_xy[1] = (float)(Ey*double_epsilon);

  // the error we need to add to linear operations (additions and substraction
  // involving the center coodinates
	cen_err = cen_err_xy[0] > cen_err_xy[1] ?  cen_err_xy[0] : cen_err_xy[1];

  double m1 = 4*(cen[0])*(cen[0]) + 4*(cen[1])*(cen[1]);
	double m2 = Ex*(fabs(-cen[0]*2)) + Ey*(fabs(-cen[1]*2));
	double m3 = Ex*Ex+Ey*Ey;

  // the error we need to add to the squared radius
  float rad2_err = (float)(2*double_epsilon*(4*rad2_tmp + E1 + m1/2+m2 +2*double_epsilon*(m1/2+3*m2+m3)));
	assert(rad2_err>=0);
  // the square radius of the circle with error
	rad2 = (float)(rad2_tmp + rad2_err);
  // the radius of the circle with error
  rad = (float)((1+sqrt_epsilon)*sqrt(rad2_tmp + rad2_err) + cen_err)*(1+single_epsilon);
}

int SScontainer2D::get_idx(const float* pos_f, int level)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  int level_idx = 0;
  int l = level;

  while (l)
  {
    level_idx <<= 2;

    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    if (pos_f[0] < r_mid[0])
    {
      r_max[0] = r_mid[0];
      level_idx |= 1;
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (pos_f[1] < r_mid[1])
    {
      r_max[1] = r_mid[1];
    }
    else
    {
      r_min[1] = r_mid[1];
      level_idx |= 2;
    }
    l--;
  }
  return level_offset[level]+level_idx;
}

int SScontainer2D::get_idx(const double* pos_d, int level)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  int level_idx = 0;
  int l = level;

  while (l)
  {
    level_idx <<= 2;

    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    if (pos_d[0] < r_mid[0])
    {
      r_max[0] = r_mid[0];
      level_idx |= 1;
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (pos_d[1] < r_mid[1])
    {
      r_max[1] = r_mid[1];
    }
    else
    {
      r_min[1] = r_mid[1];
      level_idx |= 2;
    }
    l--;
  }
  return level_offset[level]+level_idx;
}

void SScontainer2D::get_mid(float* mid_f, const float* pos_f, int idx)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  while (idx)
  {
    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    if (pos_f[0] < r_mid[0])
    {
      r_max[0] = r_mid[0];
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (pos_f[1] < r_mid[1])
    {
      r_max[1] = r_mid[1];
    }
    else
    {
      r_min[1] = r_mid[1];
    }
    idx = (idx-1) / 4;
  }

  mid_f[0] = (r_min[0] + r_max[0])/2;
  mid_f[1] = (r_min[1] + r_max[1])/2;
}

void SScontainer2D::get_min_max(float* min_f, float* max_f, const float* pos_f, int idx)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  int parent_idx = 0;
  while (parent_idx != idx)
  {

    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    if (pos_f[0] < r_mid[0])
    {
      r_max[0] = r_mid[0];
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (pos_f[1] < r_mid[1])
    {
      r_max[1] = r_mid[1];
    }
    else
    {
      r_min[1] = r_mid[1];
    }
    idx = (idx-1) / 4;
  }

  min_f[0] = r_min[0];
  min_f[1] = r_min[1];

  max_f[0] = r_max[0];
  max_f[1] = r_max[1];
}

void SScontainer2D::get_min_max(float* min_f, float* max_f, int idx)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  int parent_idx = 0;
  while (parent_idx != idx)
  {
    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    int child_idx = idx;
    while ((child_idx-1)/4 != parent_idx)
    {
      child_idx = (child_idx-1) / 4;
    }

    int child = child_idx - (4*parent_idx) - 1;

    if (child & 1)
    {
      r_max[0] = r_mid[0];
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (child & 2)
    {
      r_min[1] = r_mid[1];
    }
    else
    {
      r_max[1] = r_mid[1];
    }
    parent_idx = child_idx;
  }

  min_f[0] = r_min[0];
  min_f[1] = r_min[1];

  max_f[0] = r_max[0];
  max_f[1] = r_max[1];
}

void SScontainer2D::get_mid(float* mid_f, const double* pos_d, int idx)
{
  float r_min[2];
  float r_max[2];
  float r_mid[2];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];

  while (idx)
  {
    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;

    if (pos_d[0] < r_mid[0])
    {
      r_max[0] = r_mid[0];
    }
    else
    {
      r_min[0] = r_mid[0];
    }
    if (pos_d[1] < r_mid[1])
    {
      r_max[1] = r_mid[1];
    }
    else
    {
      r_min[1] = r_mid[1];
    }
    idx = (idx-1) / 4;
  }

  mid_f[0] = (r_min[0] + r_max[0])/2;
  mid_f[1] = (r_min[1] + r_max[1])/2;
}

/* 
Returns TRUE iff rectangle is inside of half_plane

static bool check_intersect(const float* r_min, const float* r_max, const float* p_nor,  const float* p_pnt)
{
  if (p_nor[0]*r_min[0]+p_nor[1]*r_min[1] >= p_pnt[0]) return true;
  if (p_nor[0]*r_min[0]+p_nor[1]*r_max[1] >= p_pnt[0]) return true;
  if (p_nor[0]*r_max[0]+p_nor[1]*r_min[1] >= p_pnt[0]) return true;
  if (p_nor[0]*r_max[0]+p_nor[1]*r_max[1] >= p_pnt[0]) return true;
  return false;
}

static bool intersect(SScell* cell, const float* p_nor, const float* p_pnt)
{
}
*/

static bool intersect_point(SScell* cell, const float* p_pos)
{
  if (cell->num_children) // check if cell has children
  {
    if (p_pos[1] < cell->r_mid[1]) // check if i need to check child 0 or child 1
    {
      if (p_pos[0] < cell->r_mid[0]) // do i need to check child 1
      {
        if (cell->child[1] && intersect_point(cell->child[1], p_pos)) return true;
      }
      else // i need to check child 0
      {
        if (cell->child[0] && intersect_point(cell->child[0], p_pos)) return true;
      }
    }
    else // i need to check child 2 or child 3
    {
      if (p_pos[0] < cell->r_mid[0]) // do i need to check child 3
      {
        if (cell->child[3] && intersect_point(cell->child[3], p_pos)) return true;
      }
      else // i need to check child 2
      {
        if (cell->child[2] && intersect_point(cell->child[2], p_pos)) return true;
      }
    }
    return false;
  }
  else
  {
    return true;
  }
}

static bool intersect_point(SScell* cell, const double* p_pos)
{
  if (cell->num_children) // check if cell has children
  {
    if (p_pos[1] < cell->r_mid[1]) // check if i need to check child 0 or child 1
    {
      if (p_pos[0] < cell->r_mid[0]) // do i need to check child 1
      {
        if (cell->child[1] && intersect_point(cell->child[1], p_pos)) return true;
      }
      else // i need to check child 0
      {
        if (cell->child[0] && intersect_point(cell->child[0], p_pos)) return true;
      }
    }
    else // i need to check child 2 or child 3
    {
      if (p_pos[0] < cell->r_mid[0]) // do i need to check child 3
      {
        if (cell->child[3] && intersect_point(cell->child[3], p_pos)) return true;
      }
      else // i need to check child 2
      {
        if (cell->child[2] && intersect_point(cell->child[2], p_pos)) return true;
      }
    }
    return false;
  }
  else
  {
    return true;
  }
}

bool SScontainer2D::is_finalized(const float* p_pos)
{
  return (root == 0) || (intersect_point((SScell*)root, p_pos) == false);
}

/*
bool SScontainer2D::is_finalized(const double* p_pos)
{
  if (root == 0) return true;
  
  int cell_idx = get_idx(p_pos, 5);

  my_cell_hash::iterator hash_element = cell_hash->find((cell_idx-1)/4);
  while(hash_element == cell_hash->end())
  {
    cell_idx = (cell_idx-1)/4;
    hash_element = cell_hash->find((cell_idx-1)/4);
  }
  SScell* cell = (*hash_element).second;
  if (cell->num_children == 0)
  {
    return false;
  }
  if (cell->child[cell_idx - (cell->idx*4)-1] == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}
*/


bool SScontainer2D::is_finalized(const double* p_pos)
{
  return (root == 0) || (intersect_point((SScell*)root, p_pos) == false);
}

/* 
Adapted from
Fast Circle-Rectangle Intersection Checking
by Clifford A. Shaffer
from "Graphics Gems", Academic Press, 1990

Returns TRUE iff rectangle r_min/r_max intersects
circle with centerpoint c_cen and radius c_rad.
*/

bool check_intersect_circle_no_error(const float* r_min, const float* r_max, const float* c_cen, float c_rad)
{
  float c_rad_squared = c_rad * c_rad;
  float r_minx = r_min[0] - c_cen[0];
  float r_miny = r_min[1] - c_cen[1];
  float r_maxx = r_max[0] - c_cen[0];
  float r_maxy = r_max[1] - c_cen[1];

  if (r_maxx < 0) 			// R to left of circle center
    if (r_maxy < 0) 		// R in lower left corner
      return ((r_maxx * r_maxx + r_maxy * r_maxy) < c_rad_squared);
    else if (r_miny > 0) 	// R in upper left corner
      return ((r_maxx * r_maxx + r_miny * r_miny) < c_rad_squared);
    else 					// R due West of circle
      return((-r_maxx) < c_rad);
  else if (r_minx > 0)  	// R to right of circle center
   	if (r_maxy < 0) 	// R in lower right corner
     	return ((r_minx * r_minx + r_maxy * r_maxy) < c_rad_squared);
    else if (r_miny > 0)  	// R in upper right corner
     	return ((r_minx * r_minx + r_miny * r_miny) < c_rad_squared);
    else 				// R due East of circle
     	return (r_minx < c_rad);
  else				// R on circle vertical centerline
   	if (r_maxy < 0) 	// R due South of circle
     	return ((-r_maxy) < c_rad);
    else if (r_miny > 0)  	// R due North of circle
     	return (r_miny < c_rad);
    else 				// R contains circle centerpoint
     	return true;
}

/* 
Same as above but with error treatment for floating-point round-off.
*/
bool check_intersect_circle(const float* r_min, const float* r_max, SScircle* circle)
{
  float r_minx = r_min[0] - circle->cen[0] + circle->cen_err;
  float r_miny = r_min[1] - circle->cen[1] + circle->cen_err;
  float r_maxx = r_max[0] - circle->cen[0] - circle->cen_err;
  float r_maxy = r_max[1] - circle->cen[1] - circle->cen_err;

  if (r_maxx < 0) 			// R to left of circle center
    if (r_maxy < 0) 		// R in lower left corner
      return ((r_maxx * r_maxx + r_maxy * r_maxy) < circle->rad2);
    else if (r_miny > 0) 	// R in upper left corner
      return ((r_maxx * r_maxx + r_miny * r_miny) < circle->rad2);
    else 					// R due West of circle
      return((-r_maxx) < circle->rad);
  else if (r_minx > 0)  	// R to right of circle center
   	if (r_maxy < 0) 	// R in lower right corner
     	return ((r_minx * r_minx + r_maxy * r_maxy) < circle->rad2);
    else if (r_miny > 0)  	// R in upper right corner
     	return ((r_minx * r_minx + r_miny * r_miny) < circle->rad2);
    else 				// R due East of circle
     	return (r_minx < circle->rad);
  else				// R on circle vertical centerline
   	if (r_maxy < 0) 	// R due South of circle
     	return ((-r_maxy) < circle->rad);
    else if (r_miny > 0)  	// R due North of circle
     	return (r_miny < circle->rad);
    else 				// R contains circle centerpoint
     	return true;
}

#define CHILD0    ((1<<0))
#define CHILD1    ((1<<1))
#define CHILD2    ((1<<2))
#define CHILD3    ((1<<3))
#define CHILD02   ((1<<2)|(1<<0)) // large X
#define CHILD13   ((1<<3)|(1<<1)) // small X
#define CHILD23   ((1<<3)|(1<<2)) // large Y
#define CHILD01   ((1<<1)|(1<<0)) // small Y

const static bool intersect_circle_unrolled(SScell* cell, SScircle* circle)
{
  int num_children = cell->num_children;

  if (num_children) // check if cell has children
  {
    if (num_children & CHILD02) // do children 02 with large X exist
    {
      if ((circle->cen[0] + circle->rad) > cell->r_mid[0]) /// does circle intersect large X
      {
        if (num_children & CHILD2) // does child 2 with large X and large Y exist
        {
          if ((circle->cen[1] + circle->rad) > cell->r_mid[1]) // does circle intersect large Y
          {
            if (intersect_circle_unrolled(cell->child[2], circle)) return true;
            if (num_children & CHILD0) // does child 0 with large X and small Y also exist
            {
              if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // does circle also intersect small Y
              {
                if (intersect_circle_unrolled(cell->child[0], circle)) return true;
              }
              else // remove children with small Y from further consideration
              {
                num_children &= (~CHILD01);
              }
            }
          }
          else // remove children with large Y from further consideration
          {
            num_children &= (~CHILD23); // that means circle must intersect small Y
            if (num_children & CHILD0) // does child 0 with large X and small Y exist
            {
              if (intersect_circle_unrolled(cell->child[0], circle)) return true;
            }
          }
        }
        else // if child 2 (with large Y) does not exist then child 0 (with small Y) must exist
        {
          assert(num_children & CHILD0);
          if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // does circle intersect small Y
          {
            if (intersect_circle_unrolled(cell->child[0], circle)) return true;
          }
          else // remove children with small Y from further consideration
          {
            num_children &= (~CHILD01);
          }
        }
      }
      else // remove children with large X from further consideration
      {
        num_children &= (~CHILD02); // that means circle must intersect small X
      }
    }
    if (num_children & CHILD13) // do children 13 with small X exist
    {
      if ((circle->cen[0] - circle->rad) < cell->r_mid[0]) /// does circle intersect small X
      {
        if (num_children & CHILD3) // does child 3 with small X and large Y exist
        {
          if ((circle->cen[1] + circle->rad) > cell->r_mid[1]) // does circle intersect large Y
          {
            if (intersect_circle_unrolled(cell->child[3], circle)) return true;
            if (num_children & CHILD1) // does child 1 with small X and small Y also exist
            {
              if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // does circle also intersect small Y
              {
                if (intersect_circle_unrolled(cell->child[1], circle)) return true;
              }
              else // remove children with small Y from further consideration
              {
                num_children &= (~CHILD01);
              }
            }
          }
          else // remove children with large Y from further consideration
          {
            num_children &= (~CHILD23); // that means circle must intersect small Y
            if (num_children & CHILD1) // does child 0 with small X and small Y exist
            {
              if (intersect_circle_unrolled(cell->child[1], circle)) return true;
            }
          }
        }
        else // if child 3 (with large Y) does not exist then child 1 (with small Y) must exist
        {
          assert(num_children & CHILD1);
          if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // does circle intersect small Y
          {
            if (intersect_circle_unrolled(cell->child[1], circle)) return true;
          }
          else // remove children with small Y from further consideration
          {
            num_children &= (~CHILD01);
          }
        }
      }
      else // remove children with small X from further consideration
      {
        num_children &= (~CHILD13); // that means circle must intersect small X
      }
    }
    return false;
  }
  else
  {
    if (check_intersect_circle(cell->r_min, cell->r_max, circle))
    {
      circle->cell_idx = cell->idx;
      return true;
    }
    else
    {
      return false;
    }
  }
}

static bool intersect_circle_new(SScell* cell, SScircle* circle)
{
  int num_children = cell->num_children;

  if (num_children) // check if cell has children
  {
    if (num_children & CHILD01) // check if i need to check child 0 and/or child 1
    {
      if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // can the circle intersect child 0 and/or child 1
      {
        if (num_children & CHILD0) // check if i need to check child 0
        {
          if ((circle->cen[0] + circle->rad) > cell->r_mid[0]) // can the circle intersect child 0
          {
            if (intersect_circle_new(cell->child[0], circle)) return true; // check it.
          }
          else
          {
            num_children &= (~CHILD02);
          }
        }
        if (num_children & CHILD1) // check if i need to check child 1
        {
          if ((circle->cen[0] - circle->rad) < cell->r_mid[0]) // can the circle intersect child 1
          {
            if (intersect_circle_new(cell->child[1], circle)) return true; // check it.
          }
          else
          {
            num_children &= (~CHILD13);
          }
        }
      }
    }
    if (num_children & CHILD23) // check if i need to check child 2 and/or child 3
    {
      if ((circle->cen[1] + circle->rad) > cell->r_mid[1]) // can the circle intersect child 2 and/or child 3
      {
        if (num_children & CHILD2) // check if i need to check child 2
        {
          if ((circle->cen[0] + circle->rad) > cell->r_mid[0]) // can the circle intersect child 2
          {
            if (intersect_circle_new(cell->child[2], circle)) return true; // check it.
          }
        }
        if (num_children & CHILD3) // check if i need to check child 3
        {
          if ((circle->cen[0] - circle->rad) < cell->r_mid[0]) // can the circle intersect child 3
          {
            if (intersect_circle_new(cell->child[3], circle)) return true; // check it.
          }
        }
      }
    }
    return false;
  }
  else
  {
    if (check_intersect_circle(cell->r_min, cell->r_max, circle))
    {
      circle->cell_idx = cell->idx;
      return true;
    }
    else
    {
      return false;
    }
  }
}

static bool intersect_circle(SScell* cell, SScircle* circle)
{
  if (cell->num_children) // check if cell has children
  {
    if ((circle->cen[1] - circle->rad) < cell->r_mid[1]) // check if i need to check child 0 and/or child 1
    {
      if ((circle->cen[0] + circle->rad) > cell->r_mid[0]) // check if i need to check child 0
      {
        if (cell->child[0] && intersect_circle(cell->child[0], circle)) return true; // check if i have child 0. if yes -> recursion
      }
      if ((circle->cen[0] - circle->rad) < cell->r_mid[0]) // check if i need to check child 1
      {
        if (cell->child[1] && intersect_circle(cell->child[1], circle)) return true; // check if i have child 1. if yes -> recursion
      }
    }
    if ((circle->cen[1] + circle->rad) > cell->r_mid[1]) // check if i need to check child 2 and/or child 3
    {
      if ((circle->cen[0] + circle->rad) > cell->r_mid[0]) // check if i need to check child 2
      {
        if (cell->child[2] && intersect_circle(cell->child[2], circle)) return true; // check if i have child 2. if yes -> recursion
      }
      if ((circle->cen[0] - circle->rad) < cell->r_mid[0]) // check if i need to check child 3
      {
        if (cell->child[3] && intersect_circle(cell->child[3], circle)) return true; // check if i have child 1. if yes -> recursion
      }
    }
    return false;
  }
  else
  {
    if (check_intersect_circle(cell->r_min, cell->r_max, circle))
    {
      circle->cell_idx = cell->idx;
      return true;
    }
    else
    {
      return false;
    }
  }
}

static bool fin_cell_enlarged = false;
static float fin_cell_min_f[2];
static float fin_cell_max_f[2];

static bool just_finalized_circle(const float* c_cen, float c_rad)
{
  if (c_cen[0] - c_rad < fin_cell_min_f[0]) return false;
  if (c_cen[0] + c_rad > fin_cell_max_f[0]) return false;
  if (c_cen[1] - c_rad < fin_cell_min_f[1]) return false;
  if (c_cen[1] + c_rad > fin_cell_max_f[1]) return false;
  return true;
}

bool SScontainer2D::is_finalized(SScircle* circle)
{
  // for efficiency we should not call this function once everything is finalized
  assert (root != 0);
  // also ... when we call this function the radius should be computed
  assert (circle->rad > 0);
#ifdef COLLECT_STATISTICS
    stat_finalize_active_justtested++;
#endif
  // often these previously unsees circles are entirely in the just finalized cell
  if (just_finalized_circle(circle->cen, circle->rad))
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_justfinalized++;
#endif
    return true;
  }
  return (intersect_circle((SScell*)root, circle) == false);
}

bool SScontainer2D::is_still_finalized(SScircle* circle)
{
  // for efficiency we should not call this function once everything is finalized
  assert (root != 0);
  // also ... when we call this function the radius should be computed
  assert (circle->rad > 0);
  // does the cell that prohibited finalization last time still exist?
  my_cell_hash::iterator hash_element = cell_hash->find(circle->cell_idx);
  if ((hash_element != cell_hash->end()) && (((*hash_element).second->num_children) == 0))
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_easyintersect++;
#endif
    return false;
  }
  // previously seen circumcircles are typically only in enlarged finalized cells
#ifdef ENLARGE_FINALIZED_CELLS
  if (fin_cell_enlarged)
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_justtested++;
#endif
    if (0 && just_finalized_circle(c_cen, c_rad))
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_active_justfinalized++;
#endif
      return true;
    }
  }
#endif // ENLARGE_FINALIZED_CELLS
  if (intersect_circle((SScell*)root, circle) == false)
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_nowfinalized++;
#endif
    return true;
  }
  else
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_stillintersect++;
#endif
    return false;
  }
}

/* 
Adapted from
Fast Circle-Rectangle Intersection Checking
by Clifford A. Shaffer
from "Graphics Gems", Academic Press, 1990

Returns TRUE iff rectangle r_min/r_max intersects
circle with centerpoint c_cen and squared radius c_rad_squared.

bool check_intersect_circle_squared(const float* r_min, const float* r_max, const float* c_cen, float c_rad_squared)
{
  float r_minx = r_min[0] - c_cen[0];
  float r_miny = r_min[1] - c_cen[1];
  float r_maxx = r_max[0] - c_cen[0];
  float r_maxy = r_max[1] - c_cen[1];

  if (r_maxx < 0) 			// R to left of circle center
    if (r_maxy < 0) 		// R in lower left corner
      return ((r_maxx * r_maxx + r_maxy * r_maxy) < c_rad_squared);
    else if (r_miny > 0) 	// R in upper left corner
      return ((r_maxx * r_maxx + r_miny * r_miny) < c_rad_squared);
    else 					// R due West of circle
      return ((r_maxx * r_maxx) < c_rad_squared);
  else if (r_minx > 0)  	// R to right of circle center
   	if (r_maxy < 0) 	// R in lower right corner
     	return ((r_minx * r_minx + r_maxy * r_maxy) < c_rad_squared);
    else if (r_miny > 0)  	// R in upper right corner
     	return ((r_minx * r_minx + r_miny * r_miny) < c_rad_squared);
    else 				// R due East of circle
     	return ((r_minx * r_minx) < c_rad_squared);
  else				// R on circle vertical centerline
   	if (r_maxy < 0) 	// R due South of circle
     	return ((r_maxy * r_maxy) < c_rad_squared);
    else if (r_miny > 0)  	// R due North of circle
     	return ((r_miny * r_miny) < c_rad_squared);
    else 				// R contains circle centerpoint
     	return true;
}
*/

static void finalize_subtree(SScell* cell)
{
  // dealloc children cells
  if (cell->num_children)
  {
    for (int i = 0; i < 4; i++)
    {
      if (cell->child[i])
      {
        cell_hash->erase(cell->child[i]->idx);
        finalize_subtree(cell->child[i]);
      }
    }
  }
  deallocCell(cell);
}

static SScell* create_parent_and_siblings(int idx, bool finalize)
{
  SScell* parent;
  int parent_idx = (idx - 1) / 4;
  my_cell_hash::iterator hash_element = cell_hash->find(parent_idx);

  if (hash_element == cell_hash->end())
  {
    // did not find parent ... create the parent's parent and all siblings
    parent = create_parent_and_siblings(parent_idx, false);
  }
  else
  {
    parent = (*hash_element).second;
  }

  // compute mid point of parent cell
  parent->r_mid[0] = (parent->r_min[0] + parent->r_max[0])/2;
  parent->r_mid[1] = (parent->r_min[1] + parent->r_max[1])/2;

  // create siblings
  SScell* sibling;
  SScell* return_parent = 0;
  int sibling_idx = 4*parent_idx;

  for (int i = 0; i < 4; i++)
  {
    // compute index of siblings
    sibling_idx++;
    // is it a parental or an actual sibling (but not the cell itself)
    if (sibling_idx == idx)
    {
      if (finalize)
      {
#ifdef PRINT_DEBUG_OUTPUT
        fprintf(stderr, "finalize %d (%g,%g)/(%g,%g)\n", idx, (i&1?parent->r_min[0]:parent->r_mid[0]), (i&2?r_mid[1]:parent->r_min[1]), (i&1?r_mid[0]:parent->r_max[0]), (i&2?parent->r_max[1]:r_mid[1]));
#endif
        // copy the bounding box of the finalized cell
        if (i&1)
        {
          fin_cell_min_f[0] = parent->r_min[0];
          fin_cell_max_f[0] = parent->r_mid[0];
        }
        else
        {
          fin_cell_min_f[0] = parent->r_mid[0];
          fin_cell_max_f[0] = parent->r_max[0];
        }
        if (i&2)
        {
          fin_cell_min_f[1] = parent->r_mid[1];
          fin_cell_max_f[1] = parent->r_max[1];
        }
        else
        {
          fin_cell_min_f[1] = parent->r_min[1];
          fin_cell_max_f[1] = parent->r_mid[1];
        }
        parent->child[i] = 0;
        continue;
      }
      // allocate sibling (and return_parent)
      sibling = allocCell();
      return_parent = sibling;
    }
    else
    {
      // allocate sibling cell
      sibling = allocCell();
    }

    // give sibling its idx
    sibling->idx = sibling_idx;
    // give sibling its bounding box
    if (i&1)
    {
      sibling->r_min[0] = parent->r_min[0];
      sibling->r_max[0] = parent->r_mid[0];
    }
    else
    {
      sibling->r_min[0] = parent->r_mid[0];
      sibling->r_max[0] = parent->r_max[0];
    }
    if (i&2)
    {
      sibling->r_min[1] = parent->r_mid[1];
      sibling->r_max[1] = parent->r_max[1];
    }
    else
    {
      sibling->r_min[1] = parent->r_min[1];
      sibling->r_max[1] = parent->r_mid[1];
    }
    // insert it into the hash
    cell_hash->insert(my_cell_hash::value_type(sibling_idx, sibling));
    // attach sibling to parent
    parent->child[i] = sibling;
    parent->num_children |= (1<<i);
    // make sibling point to parent
    sibling->parent = parent;
  }
  assert(finalize || return_parent);
  return return_parent;
}

void SScontainer2D::finalize_cell(int cell_idx)
{
  my_cell_hash::iterator hash_element;
  
  hash_element = cell_hash->find(cell_idx);

  if (hash_element == cell_hash->end())
  {
    // did not find cell ... create the subtree surrounding it
    create_parent_and_siblings(cell_idx, true);
    // this cell clearly will not be enlarged
#ifdef ENLARGE_FINALIZED_CELLS
    fin_cell_enlarged = false;
#endif
//    fprintf(stderr, "finalize_cell non-existing %d = %g/%g %g/%g %d\n", cell_idx, fin_cell_min_f[0], fin_cell_min_f[1], fin_cell_max_f[0], fin_cell_max_f[1], fin_cell_enlarged);
  }
  else
  {
    // found cell
    SScell* cell = (*hash_element).second;
    // make sure this is true
    assert(cell->idx == cell_idx);
    // remove from hash
    cell_hash->erase(hash_element);
    // get parent of cell
    SScell* parent = cell->parent;
    // does parent exist?
    if (parent)
    {
      // which of its children is it
      int i = cell_idx-(4*parent->idx+1);
      // make sure this child exists
      assert(parent->child[i]);
      // it will exist no longer
      parent->child[i] = 0;
      parent->num_children &= (~(1<<i));
      // we must still have children
      assert(parent->num_children);
      // copy the bounding box of the finalized cell
      fin_cell_min_f[0] = cell->r_min[0];
      fin_cell_min_f[1] = cell->r_min[1];
      fin_cell_max_f[0] = cell->r_max[0];
      fin_cell_max_f[1] = cell->r_max[1];
      // can we enlarge it?
#ifdef ENLARGE_FINALIZED_CELLS
      if (parent->child[i^1] == 0) // in x direction?
      {
        fin_cell_min_f[0] = parent->r_min[0];
        fin_cell_max_f[0] = parent->r_max[0];
        fin_cell_enlarged = true;
      }
      else if (parent->child[i^2] == 0) // in y direction?
      {
        fin_cell_min_f[1] = parent->r_min[1];
        fin_cell_max_f[1] = parent->r_max[1];
        fin_cell_enlarged = true;
      }
      else
      {
        // the cell can be enlarged by finalizing a parent 
        if (cell->num_children)
        {
          fin_cell_enlarged = true;
        }
        else
        {
          fin_cell_enlarged = false;
        }
      }
#endif
    }
    else
    {
      // this must be the root cell
      assert(cell == (SScell*)root);
      // which will now longer exist
      root = 0;
    }
    // finalize the subtree rooted in cell
    finalize_subtree(cell);

//    fprintf(stderr, "finalize_cell existing %d = %g/%g %g/%g %d\n", cell_idx, fin_cell_min_f[0], fin_cell_min_f[1], fin_cell_max_f[0], fin_cell_max_f[1], fin_cell_enlarged);
  }
}


void SScontainer2D::iterateInit()
{
  hash_iterator = cell_hash->begin();
}

bool SScontainer2D::iterateNext()
{
  SScell* cell;

  while (hash_iterator != cell_hash->end())
  {
    cell = (*hash_iterator).second;
    hash_iterator++;
    if (cell->num_children == 0)
    {
      r_min_f = cell->r_min;
      r_max_f = cell->r_max;
      return true;
    }
  }
  return false;
}

void SScontainer2D::open(const float* bb_min_f, const float* bb_max_f)
{
  this->bb_min_f[0] = bb_min_f[0];
  this->bb_min_f[1] = bb_min_f[1];
  this->bb_max_f[0] = bb_max_f[0];
  this->bb_max_f[1] = bb_max_f[1];
  root = allocCell();
  ((SScell*)root)->idx = 0;
  ((SScell*)root)->r_min[0] = this->bb_min_f[0];
  ((SScell*)root)->r_min[1] = this->bb_min_f[1];
  ((SScell*)root)->r_max[0] = this->bb_max_f[0];
  ((SScell*)root)->r_max[1] = this->bb_max_f[1];
  cell_hash->insert(my_cell_hash::value_type(0, (SScell*)root));
}

void SScontainer2D::open(const double* bb_min_d, const double* bb_max_d)
{
  this->bb_min_f[0] = (float)bb_min_d[0];
  this->bb_min_f[1] = (float)bb_min_d[1];
  this->bb_max_f[0] = (float)bb_max_d[0];
  this->bb_max_f[1] = (float)bb_max_d[1];

  // make sure the min values are *small* enough
  if (bb_min_d[0] < this->bb_min_f[0])
  {
    float epsilon = single_epsilon;
    float result = this->bb_min_f[0] - epsilon;
    fprintf(stderr,"WARNING: bb_min_f[0] %f > %f ",this->bb_min_f[0],bb_min_d[0]);
    while (bb_min_d[0] < result)
    {
      epsilon = epsilon*2;
      result = this->bb_min_f[0] - epsilon;
    }
    this->bb_min_f[0] = result;
    fprintf(stderr,"corrected to %f 0==%d \n",this->bb_min_f[0],bb_min_d[0] < this->bb_min_f[0]);
  }
  if (bb_min_d[1] < this->bb_min_f[1])
  {
    float epsilon = single_epsilon;
    float result = this->bb_min_f[1] - epsilon;
    fprintf(stderr,"WARNING: bb_min_f[1] %g > %f ",this->bb_min_f[1],bb_min_d[1]);
    while (bb_min_d[1] < result)
    {
      epsilon = epsilon*2;
      result = this->bb_min_f[1] - epsilon;
    }
    this->bb_min_f[1] = result;
    fprintf(stderr,"corrected to %g\n",this->bb_min_f[1]);
  }
  // make sure the max values are *big* enough
  if (bb_max_d[0] > this->bb_max_f[0])
  {
    float epsilon = single_epsilon;
    float result = this->bb_max_f[0] + epsilon;
    fprintf(stderr,"WARNING: bb_max_f[0] %g < %f ",this->bb_max_f[0],bb_max_d[0]);
    while (bb_max_d[0] > result)
    {
      epsilon = epsilon*2;
      result = this->bb_max_f[0] + epsilon;
    }
    this->bb_max_f[0] = result;
    fprintf(stderr,"corrected to %g\n",this->bb_max_f[0]);
  }
  if (bb_max_d[1] > this->bb_max_f[1])
  {
    float epsilon = single_epsilon;
    float result = this->bb_max_f[1] - epsilon;
    fprintf(stderr,"WARNING: bb_max_f[1] %g < %f ",this->bb_max_f[1],bb_max_d[1]);
    while (bb_max_d[1] > result)
    {
      epsilon = epsilon*2;
      result = this->bb_max_f[1] + epsilon;
    }
    this->bb_max_f[1] = result;
    fprintf(stderr,"corrected to %g\n",this->bb_max_f[1]);
  }

  root = allocCell();
  ((SScell*)root)->idx = 0;
  ((SScell*)root)->r_min[0] = this->bb_min_f[0];
  ((SScell*)root)->r_min[1] = this->bb_min_f[1];
  ((SScell*)root)->r_max[0] = this->bb_max_f[0];
  ((SScell*)root)->r_max[1] = this->bb_max_f[1];
  cell_hash->insert(my_cell_hash::value_type(0, (SScell*)root));
}

void SScontainer2D::close()
{
  if (root) finalize_subtree((SScell*)root);
  cell_hash->clear();
}

SScontainer2D::SScontainer2D()
{
  cell_hash = new my_cell_hash;
  root = 0;

  // set hard-coded epsilons
	single_epsilon = 1.0;
	{ for (int i=0; i<23; i++){ single_epsilon*=0.5;}}
	double_epsilon = 1.0;
	{ for (int i=0; i<52; i++){ double_epsilon*=0.5;}}
	sqrt_epsilon = single_epsilon;

  // compute level offsets
  level_offset[0] = 0;
  for (int l = 0; l < 19; l++) level_offset[l+1] = level_offset[l] + ((1<<l)*(1<<l));
}

SScontainer2D::~SScontainer2D()
{
#ifdef COLLECT_STATISTICS
  fprintf(stderr, "justtested %d justfinalized %d fulltested %d\n", stat_finalize_active_justtested, stat_finalize_active_justfinalized, stat_finalize_active_justtested-stat_finalize_active_justfinalized);
  fprintf(stderr, "easyintersect %d nowfinalized %d stillintersect %d\n", stat_finalize_active_easyintersect, stat_finalize_active_nowfinalized, stat_finalize_active_stillintersect);
#endif
  delete cell_hash;
}
