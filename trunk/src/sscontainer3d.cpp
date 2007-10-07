/*
===============================================================================

  FILE:  SScontainer3D.cpp
  
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
#include "sscontainer3d.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

#ifdef COLLECT_STATISTICS
static int stat_sphere_initialized = 0;
static int stat_finalize_justtested = 0;
static int stat_finalize_justfinalized = 0;
static int stat_box_overlap_true = 0;
static int stat_box_overlap_false = 0;
static int stat_box_overlap1_true = 0;
static int stat_box_overlap1_false = 0;
#endif

#include <ext/hash_map>

typedef __gnu_cxx::hash_map<int, SScell3D*> my_cell_hash;

// local variables
static my_cell_hash* cell_hash;
static my_cell_hash::iterator hash_iterator;
static float grid_cell_diagonal[16];
static int deepest_level = -1;

// statistics

static int cell_buffer_size;
static int cell_buffer_maxsize;

// for external query

SScontainer3D* this_container = 0;

int getDataFromLeaf(const float* pos_f)
{
  return this_container->get_leaf_data(pos_f);
}

// efficient memory management

static int cell_buffer_alloc = 16;
static SScell3D* cell_buffer_next = (SScell3D*)-1;

static SScell3D* allocCell()
{
  if (cell_buffer_next == (SScell3D*)-1)
  {
    cell_buffer_next = (SScell3D*)malloc(sizeof(SScell3D)*cell_buffer_alloc);
    if (cell_buffer_next == 0)
    {
      fprintf(stderr,"ERROR: malloc for cell buffer failed\n");
      return 0;
    }
    for (int i = 0; i < cell_buffer_alloc; i++)
    {
      cell_buffer_next[i].buffer_next = &(cell_buffer_next[i+1]);
    }
    cell_buffer_next[cell_buffer_alloc-1].buffer_next = (SScell3D*)-1;
    cell_buffer_alloc = 2*cell_buffer_alloc;
  }
  // get index of next available cell
  SScell3D* cell = cell_buffer_next;
  cell_buffer_next = cell->buffer_next;
 
  // clean cell
  cell->buffer_next = 0;
  cell->num_children = 0;
  cell->parent = 0;
  cell->data = -1;

  cell_buffer_size++; if (cell_buffer_size > cell_buffer_maxsize) cell_buffer_maxsize = cell_buffer_size;

  return cell;
}

static void deallocCell(SScell3D* cell)
{
  cell->buffer_next = cell_buffer_next;
  cell_buffer_next = cell;
  cell_buffer_size--;
}

/*
	//============== Description of the error bounds=================================/
	// Two kinds of error bounds are computed: the center errors in each coordinate, and 
	//   the errors in the squared radius. Instead of storing the errors, they are added 
	//   into rad and rad2 to minimize code changes.
	// The errors are added into _rad_ and _rad2_ differently:
	//   _rad2_ is added with only the squared radius error
	//   _rad_  is (before taking the sqrt) added with both the squared radius error AND the center error
	// The reason is that if sphere collision test involves only linear operations 
	//   the radius and center (such as bounding box test), then the radius and center errors
	//   can be simply added. 
	// However, if the test is any more complicated, such as computing distance
	//   from point to center, then the center error must be known.
	// The function isOutSideSphere(float x, float y, float z) demonstrates the use

	// Safe test for testing if the point (x,y,z) is outside the sphere:
	//   if the true answer is yes, the test always returns yes, but not always vice versa.
	bool isOutsideSphere(float x, float y, float z)
	{
		double dx =  fabs(x -  cen[0] ) - center_err;  if (dx<0) dx=0;
		double dy =  fabs(y -  cen[1] ) - center_err;  if (dy<0) dy=0;
		double dz =  fabs(z -  cen[2] ) - center_err;  if (dz<0) dz=0;
		double d = (dx*dx+dy*dy+dz*dz) -  rad2;
		return d>=0;
	}
*/

static float single_epsilon;
static double double_epsilon;
static float sqrt_epsilon;

void SSsphere::initialize(const double* v0, const double* v1, const double* v2, const double* v3)
{
#define LIFT(x,y,z) (x*x+y*y+z*z)
  assert(v0[3] == LIFT(v0[0],v0[1],v0[2]));
 	assert(v1[3] == LIFT(v1[0],v1[1],v1[2]));
	assert(v2[3] == LIFT(v2[0],v2[1],v2[2]));
	assert(v3[3] == LIFT(v3[0],v3[1],v3[2]));
#undef LIFT

#ifdef COLLECT_STATISTICS
  stat_sphere_initialized++;
#endif

  double x0, y0, z0, sq0, x1, y1, z1, sq1, x2, y2, z2, sq2;
  x0 = v0[0] - v3[0]; y0 = v0[1] - v3[1]; z0 = v0[2] - v3[2]; sq0 = v0[3] - v3[3];
  x1 = v1[0] - v3[0]; y1 = v1[1] - v3[1]; z1 = v1[2] - v3[2]; sq1 = v1[3] - v3[3];
  x2 = v2[0] - v3[0]; y2 = v2[1] - v3[1]; z2 = v2[2] - v3[2]; sq2 = v2[3] - v3[3];

  double xy, xz, xs, yz, ys, zs; // 2x2 minors
	double xyE, xzE, xsE, yzE, ysE, zsE; // 2x2 minors errors
#define DET2_AND_ERR(a,b,p,q,i,j) {double t1 = (i##p)*(j##q); double t2 = (j##p)*(i##q); a=t1-t2; b=fabs(t1)+fabs(t2);}
	DET2_AND_ERR(xy,xyE,0,1,x,y); 
	DET2_AND_ERR(xz,xzE,0,1,x,z); 
	DET2_AND_ERR(yz,yzE,0,1,y,z); 
	DET2_AND_ERR(xs,xsE,0,1,x,sq); 
	DET2_AND_ERR(ys,ysE,0,1,y,sq);
	DET2_AND_ERR(zs,zsE,0,1,z,sq);
#undef DET2_AND_ERR

  double Mx,My,Mz,Mq;
	double nx, ny, nz, nq; //error terms
	Mx = -y2*zs +z2*ys -sq2*yz;  nx = fabs(y2*zsE) + fabs(z2*ysE) + fabs(sq2*yzE);
	My =  x2*zs -z2*xs +sq2*xz;  ny = fabs(x2*zsE) + fabs(z2*xsE) + fabs(sq2*xzE);  
	Mz = -x2*ys +y2*xs -sq2*xy;  nz = fabs(x2*ysE) + fabs(y2*xsE) + fabs(sq2*xyE); 
	Mq =  x2*yz -y2*xz + z2*xy;  nq = fabs(x2*yzE) + fabs(y2*xzE) + fabs( z2*xyE);

  assert(Mq!=0);

  // center of circle
	double cen_double[3];

  cen_double[0] = (-Mx/(2.0*Mq));
	cen_double[1] = (-My/(2.0*Mq));
	cen_double[2] = (-Mz/(2.0*Mq));
	
	this->cen[0] = cen_double[0];
	this->cen[1] = cen_double[1];
	this->cen[2] = cen_double[2];


	cen_double[0] = fabs(cen_double[0]);
	cen_double[1] = fabs(cen_double[1]);
	cen_double[2] = fabs(cen_double[2]);

	double M1 = -Mx*v3[0]- My*v3[1]- Mz*v3[2]- Mq*v3[3];

  double Mx2 = Mx*Mx;
	double My2 = My*My;
	double Mz2 = Mz*Mz;
	double Mq2 = Mq*Mq;

  // squared radius of circle

  double rad2_tmp = fabs((Mx2 + My2 + Mz2)/(4*Mq2) -M1/Mq);

	M1 = fabs(M1);
	Mx = fabs(Mx);
	My = fabs(My);
	Mz = fabs(Mz);
	Mq = fabs(Mq);

	double topx = Mx*Mq+4*(Mq*nx+ Mx*nq);
	double topy = My*Mq+4*(Mq*ny+ My*nq);
	double topz = Mz*Mq+4*(Mq*nz+ Mz*nq);

  double bot = Mq*(Mq-2*double_epsilon*(Mq+8*double_epsilon*nq));

  if (bot<=0)
  {					
		fprintf(stderr, "The divisor in the error bound computation is zero. Using double_epsilon.\n");
		bot = double_epsilon;
	}

	assert(bot<=Mq2); //should be lower bound to maxmize top/bot

	double E1 =  2*(11*M1*Mq + 8*M1*nq)/bot;          
	double Ex =  5*topx/bot;
	double Ey =  5*topy/bot;
	double Ez =  5*topz/bot;

	float cen_err_xyz[3]; 
	cen_err_xyz[0] = (float)((cen_double[0]+Ex)*double_epsilon + single_epsilon*fabs(cen[0])); 
	cen_err_xyz[1] = (float)((cen_double[1]+Ey)*double_epsilon + single_epsilon*fabs(cen[1]));
	cen_err_xyz[2] = (float)((cen_double[2]+Ez)*double_epsilon + single_epsilon*fabs(cen[2]));

	double cen_err_temp; 
  // the error we need to add to linear operations (additions and substraction
  // involving the center coodinates
	cen_err_temp = cen_err_xyz[0] > cen_err_xyz[1] ? cen_err_xyz[0] : cen_err_xyz[1];
	cen_err_temp = cen_err_temp > cen_err_xyz[2] ? cen_err_temp : cen_err_xyz[2];
	cen_err = cen_err_temp*(1+single_epsilon);

  double m1 = Ex*(cen_double[0]*2) + Ey*(cen_double[1]*2) + Ez*(cen_double[2]*2); 
	double m2 = double_epsilon*(Ex*Ex+Ey*Ey+Ez*Ez);
 
  // the error we need to add to the squared radius
	double rad2_err = double_epsilon*(10*rad2_tmp + 4*E1 + 2*m1 + m2);
	assert(rad2_err>=0);
  // the square radius of the circle with error
	rad2 = (float)(rad2_tmp + rad2_err);
  // the radius of the circle with error
	rad = (float)((1+sqrt_epsilon) * sqrt(rad2) + cen_err)*(1+single_epsilon); 
}

// what should happen to points that lie *exactly* on the
// boundary between two cells.  we simply declare those
// point that are exactly on the bounday consistently as
int SScontainer3D::get_idx(const float* pos_f, int level)
{
  float r_min[3];
  float r_max[3];
  float r_mid[3];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_min[2] = bb_min_f[2];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];
  r_max[2] = bb_max_f[2];

  int level_idx = 0;
  int l = level;

  while (l)
  {
    level_idx <<= 3;

    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;
    r_mid[2] = (r_min[2] + r_max[2])/2;

    if (pos_f[0] < r_mid[0]) // only if strictly less
    {
      r_max[0] = r_mid[0];
      level_idx |= 1;
    }
    else  // if equal or larger
    {
      r_min[0] = r_mid[0];
    }
    if (pos_f[1] < r_mid[1]) // only if strictly less
    {
      r_max[1] = r_mid[1];
    }
    else // if equal or larger
    {
      r_min[1] = r_mid[1];
      level_idx |= 2;
    }
    if (pos_f[2] < r_mid[2]) // only if strictly less
    {
      r_max[2] = r_mid[2];
      level_idx |= 4;
    }
    else // if equal or larger
    {
      r_min[2] = r_mid[2];
    }
    l--;
  }
  return level_offset[level]+level_idx;
}

int SScontainer3D::get_idx(const double* pos_d, int level)
{
  float r_min[3];
  float r_max[3];
  float r_mid[3];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_min[2] = bb_min_f[2];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];
  r_max[2] = bb_max_f[2];

  int level_idx = 0;
  int l = level;

  while (l)
  {
    level_idx <<= 3;

    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;
    r_mid[2] = (r_min[2] + r_max[2])/2;

    if (pos_d[0] < r_mid[0]) // only if strictly less
    {
      r_max[0] = r_mid[0];
      level_idx |= 1;
    }
    else  // if equal or larger
    {
      r_min[0] = r_mid[0];
    }
    if (pos_d[1] < r_mid[1]) // only if strictly less
    {
      r_max[1] = r_mid[1];
    }
    else // if equal or larger
    {
      r_min[1] = r_mid[1];
      level_idx |= 2;
    }
    if (pos_d[2] < r_mid[2]) // only if strictly less
    {
      r_max[2] = r_mid[2];
      level_idx |= 4;
    }
    else // if equal or larger
    {
      r_min[2] = r_mid[2];
    }
    l--;
  }
  return level_offset[level]+level_idx;
}

void SScontainer3D::get_mid(float* mid_f, const float* pos_f, int idx)
{
  float r_min[3];
  float r_max[3];
  float r_mid[3];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_min[2] = bb_min_f[2];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];
  r_max[2] = bb_max_f[2];

  while (idx)
  {
    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;
    r_mid[2] = (r_min[2] + r_max[2])/2;

    if (pos_f[0] < r_mid[0]) // only if strictly less
    {
      r_max[0] = r_mid[0];
    }
    else  // if equal or larger
    {
      r_min[0] = r_mid[0];
    }
    if (pos_f[1] < r_mid[1]) // only if strictly less
    {
      r_max[1] = r_mid[1];
    }
    else // if equal or larger
    {
      r_min[1] = r_mid[1];
    }
    if (pos_f[2] < r_mid[2]) // only if strictly less
    {
      r_max[2] = r_mid[2];
    }
    else // if equal or larger
    {
      r_min[2] = r_mid[2];
    }
    idx = (idx-1) / 8;
  }

  mid_f[0] = (r_min[0] + r_max[0])/2;
  mid_f[1] = (r_min[1] + r_max[1])/2;
  mid_f[2] = (r_min[2] + r_max[2])/2;
}

void SScontainer3D::get_mid(float* mid_f, const double* pos_d, int idx)
{
  float r_min[3];
  float r_max[3];
  float r_mid[3];
  
  r_min[0] = bb_min_f[0];
  r_min[1] = bb_min_f[1];
  r_min[2] = bb_min_f[2];
  r_max[0] = bb_max_f[0];
  r_max[1] = bb_max_f[1];
  r_max[2] = bb_max_f[2];

  while (idx)
  {
    r_mid[0] = (r_min[0] + r_max[0])/2;
    r_mid[1] = (r_min[1] + r_max[1])/2;
    r_mid[2] = (r_min[2] + r_max[2])/2;

    if (pos_d[0] < r_mid[0]) // only if strictly less
    {
      r_max[0] = r_mid[0];
    }
    else  // if equal or larger
    {
      r_min[0] = r_mid[0];
    }
    if (pos_d[1] < r_mid[1]) // only if strictly less
    {
      r_max[1] = r_mid[1];
    }
    else // if equal or larger
    {
      r_min[1] = r_mid[1];
    }
    if (pos_d[2] < r_mid[2]) // only if strictly less
    {
      r_max[2] = r_mid[2];
    }
    else // if equal or larger
    {
      r_min[2] = r_mid[2];
    }
    idx = (idx-1) / 8;
  }

  mid_f[0] = (r_min[0] + r_max[0])/2;
  mid_f[1] = (r_min[1] + r_max[1])/2;
  mid_f[2] = (r_min[2] + r_max[2])/2;
}

int SScontainer3D::get_leaf_data(const float* pos_f)
{
  SScell3D* leaf = root;
  while (leaf && leaf->num_children)
  {
    int child = 0;
    if (pos_f[0] < leaf->r_mid[0]) // only if strictly less
    {
      child |= 1;
    }
    if (!(pos_f[1] < leaf->r_mid[1])) // only if not strictly less
    {
      child |= 2;
    }
    if (pos_f[2] < leaf->r_mid[2]) // only if strictly less
    {
      child |= 4;
    }
    leaf = leaf->child[child];
  }
  assert(leaf);
  return leaf->data;
}

const static bool check_intersect_point(const float* r_min, const float* r_max, const float* p_pos)
{
  if (p_pos[0] < r_min[0]) return false;
  if (p_pos[0] > r_max[0]) return false;
  if (p_pos[1] < r_min[1]) return false;
  if (p_pos[1] > r_max[1]) return false;
  if (p_pos[2] < r_min[2]) return false;
  if (p_pos[2] > r_max[2]) return false;
  return true;
}

const static bool intersect_point(const SScell3D* cell, const float* p_pos)
{
  if (cell->num_children) // check if cell has children
  {
    if (p_pos[0] < cell->r_mid[0]) // do i need to check child 1/3/5/7
    {
      if (p_pos[1] < cell->r_mid[1]) // do i need to check child 1 or 5
      {
	      if (p_pos[2] < cell->r_mid[2]) // do i need to check child 5
		    {
		      return (cell->child[5] && intersect_point(cell->child[5], p_pos));
        }
        else // need to check child 1
        {
		      return (cell->child[1] && intersect_point(cell->child[1], p_pos));
        }
      }
      else // do i need to check child 3 or 7
      {
	      if (p_pos[2] < cell->r_mid[2]) // do i need to check child 7
		    {
		      return (cell->child[7] && intersect_point(cell->child[7], p_pos));
        }
        else // need to check child 3
        {
		      return (cell->child[3] && intersect_point(cell->child[3], p_pos));
        }
      }
    }
    else // need to check child 0/2/4/6
    {
      if (p_pos[1] < cell->r_mid[1]) // do i need to check child 0 or 4
      {
	      if (p_pos[2] < cell->r_mid[2]) // do i need to check child 4
		    {
		      return (cell->child[4] && intersect_point(cell->child[4], p_pos));
        }
        else // need to check child 0
        {
		      return (cell->child[0] && intersect_point(cell->child[0], p_pos));
        }
      }
      else // need to check child 2 or 6
      {
	      if (p_pos[2] < cell->r_mid[2]) // do i need to check child 6
		    {
		      return (cell->child[6] && intersect_point(cell->child[6], p_pos));
        }
        else // need to check child 2
        {
		      return (cell->child[2] && intersect_point(cell->child[2], p_pos));
        }
      }
    }
  }
  else
  {
    return check_intersect_point(cell->r_min, cell->r_max, p_pos);
  }
}

bool SScontainer3D::is_point_finalized(const float* p_pos) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root); 
  // intersect the point with the finalization grid
  return (intersect_point(root, p_pos) == false);
}

bool SScontainer3D::is_point_finalized(const double* p_pos) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root);
  // cast the double-precision float to single-precision
  float pos[3];
  pos[0] = (float)p_pos[0];
  pos[1] = (float)p_pos[1];
  pos[2] = (float)p_pos[2];
  // intersect the point with the finalization grid
  return (intersect_point(root, pos) == false);
}

static const bool check_intersect_halfplane(const float* r_min, const float* r_max, const float* p_nor, const float* p_pos)
{
  float temp[3];
  VecSubtract3fv(temp, p_pos, r_min[0], r_min[1], r_min[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_max[0], r_max[1], r_max[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_min[0], r_max[1], r_min[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_max[0], r_min[1], r_max[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_max[0], r_min[1], r_min[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_min[0], r_max[1], r_max[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_min[0], r_max[1], r_min[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  VecSubtract3fv(temp, p_pos, r_max[0], r_max[1], r_min[2]);
  if (VecDotProd3fv(temp, p_nor) > 0) return true;
  return false;
}

static const bool intersect_halfplane(const SScell3D* cell, const float* p_nor, const float* p_pos)
{
  if (cell->num_children) // check if cell has children
  {
    for (int c = 0; c < 8; c++)
    {
      if (cell->child[c] && intersect_halfplane(cell->child[c], p_nor, p_pos)) return true;
    }
    return false;
  }
  else
  {
    return check_intersect_halfplane(cell->r_min, cell->r_max, p_nor, p_pos);
  }
}

bool SScontainer3D::is_halfplane_finalized(const float* p0_pos, const float* p1_pos, const float* p2_pos) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root);
  float normal[3];
  VecCcwNormal3fv(normal, p0_pos, p1_pos, p2_pos);
  return (intersect_halfplane(root, normal, p0_pos) == false);
}

const static bool cell_overlaps_sphere1(const SScell3D* cell, const SSsphere* sphere) 
{ 
  double s, d = 0;
  // make sure we are not testing useless boxes
  assert((sphere->cen[0] + sphere->rad >= cell->r_min[0]) && (sphere->cen[0] - sphere->rad <= cell->r_max[0]));
  assert((sphere->cen[1] + sphere->rad >= cell->r_min[1]) && (sphere->cen[1] - sphere->rad <= cell->r_max[1]));
  assert((sphere->cen[2] + sphere->rad >= cell->r_min[2]) && (sphere->cen[2] - sphere->rad <= cell->r_max[2]));  

	double r_min, r_max ;

  // compute the square of the distance from the sphere to the box
  for (int i = 0; i < 3; i++) 
  { 
		r_min = (cell->r_min[i] - sphere->cen_err);
		r_max = (cell->r_max[i] + sphere->cen_err);


    if ( sphere->cen[i] <  r_min  )
    {
      s = sphere->cen[i] -  r_min ;
      d += s * s;
    }
    else if ( sphere->cen[i] >  r_max  )
    {
      s = sphere->cen[i] - r_max  ;
      d += s * s; 
    }
  } 

	//CAUTION: this is fragile code
	// and relies on the assumption that the double
	// precision has 53 bit code which is just few more
	// bits more than 2*23 that d can be computed exactly
	// otherwise we need to set d = d * (1 - n*double_epsilons), n<3
 
  if (d <= sphere->rad2)
  {
#ifdef COLLECT_STATISTICS
    stat_box_overlap1_true++;
#endif
    return true;
  }
  else
  {
#ifdef COLLECT_STATISTICS
    stat_box_overlap1_false++;
#endif
    return false;
  }
}

const static bool cell_overlaps_sphere(const SScell3D* cell, const SSsphere* sphere) 
{ 
  double  s, d = 0;
  // make sure we are not testing useless boxes
  assert((sphere->cen[0] + sphere->rad >= cell->r_min[0]) && (sphere->cen[0] - sphere->rad <= cell->r_max[0]));
  assert((sphere->cen[1] + sphere->rad >= cell->r_min[1]) && (sphere->cen[1] - sphere->rad <= cell->r_max[1]));
  assert((sphere->cen[2] + sphere->rad >= cell->r_min[2]) && (sphere->cen[2] - sphere->rad <= cell->r_max[2]));  

	double r_min, r_max;

  // compute the square of the distance from the sphere to the box
  for (int i = 0; i < 3; i++) 
  { 
		r_min  = cell->r_min[i] - sphere->cen_err;
		r_max  = cell->r_max[i] + sphere->cen_err;

    if ( sphere->cen[i] <  r_min  )
    {
      s = sphere->cen[i] -  r_min ;
      d += s * s;
    }
    else if ( sphere->cen[i] >  r_max  )
    {
      s = sphere->cen[i] - r_max ;
      d += s * s; 
    }
  } 

	//CAUTION: this is fragile code
	// and relies on the assumption that the double
	// precision has 53 bit code which is just few more
	// bits more than 2*23 that d can be computed exactly
	// otherwise we need to set d = d * (1 - n*double_epsilons), n<3

  return (d <= sphere->rad2);
}

// this function takes a cell and a sphere as input and check if there is some leaf cell that
// overlaps with the sphere. if yes, it uses the last three arguments cell_idx, tet_next, and
// tet_idx to link the corresponding tetrahedron into a linked list maintained by that gridcell

#define CHILD0123 ((1<<3)|(1<<2)|(1<<1)|(1<<0)) // large Z
#define CHILD4567 ((1<<7)|(1<<6)|(1<<5)|(1<<4)) // small Z
#define CHILD2367 ((1<<7)|(1<<6)|(1<<3)|(1<<2)) // large Y
#define CHILD0145 ((1<<5)|(1<<4)|(1<<1)|(1<<0)) // small Y

#define CHILD0246 ((1<<6)|(1<<4)|(1<<2)|(1<<0)) // large X
#define CHILD26   ((1<<6)|(1<<2))
#define CHILD2    ((1<<2))
#define CHILD6    ((1<<6))
#define CHILD04   ((1<<4)|(1<<0))
#define CHILD0    ((1<<0))
#define CHILD4    ((1<<4))
#define CHILD1357 ((1<<7)|(1<<5)|(1<<3)|(1<<1)) // small X
#define CHILD37   ((1<<7)|(1<<3))
#define CHILD3    ((1<<3))
#define CHILD7    ((1<<7))
#define CHILD15   ((1<<5)|(1<<1))
#define CHILD1    ((1<<1))
#define CHILD5    ((1<<5))

/*
const static bool intersect_sphere_new_notunrolled(SScell3D* cell, SSsphere* sphere, int* tet_next, int tet_idx) 
{
  int num_children = cell->num_children;

  if (num_children)
  {
    if ((grid_cell_diagonal[cell->level] > sphere->rad) || box_overlaps_sphere1(cell->r_min, cell->r_max, s_cen, sphere->rad))
    {
      if (num_children & CHILD0246) // do children 0246 with large X exist
      {
        if ((sphere->cen[0] + sphere->rad) > cell->r_mid[0]) // does sphere intersect large X
        {
          if (num_children & CHILD26) // do children 26 with large X and large Y exist
          {
            if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // does sphere intersect large Y
            {
              if (num_children & CHILD2) // does child 2 with large X and large Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[2], sphere, s_idx)) return true;
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123);
                }
              }
              if (num_children & CHILD6) // does child 6 with large X and large Y and small Z exist
              {
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[6], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with large Y from further consideration
            {
              num_children &= (~CHILD2367); 
            }
          }
          if (num_children & CHILD04) // do children 04 with large X and small Y exist
          {
            if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // does sphere intersect small Y
            {
              if (num_children & CHILD0) // does child 0 with large X and small Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[0], sphere, s_idx)) return true;
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123);
                }
              }
              if (num_children & CHILD4) // does child 4 with large X and small Y and small Z exist
              {
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[4], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with small Y from further consideration
            {
              num_children &= (~CHILD0145); 
            }
          }
        }
        else // remove children with large X from further consideration
        {
          num_children &= (~CHILD0246);
        }
      }

      if (num_children & CHILD1357) // do children 1357 with small X exist
      {
        if ((sphere->cen[0] - sphere->rad) < cell->r_mid[0]) // does sphere intersect small X
        {
          if (num_children & CHILD37) // do children 37 with small X and large Y exist
          {
            if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // does sphere intersect large Y
            {
              if (num_children & CHILD3) // does child 3 with small X and large Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[3], sphere, s_idx)) return true;
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123);
                }
              }
              if (num_children & CHILD7) // does child 7 with small X and large Y and small Z exist
              {
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[7], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with large Y from further consideration
            {
              num_children &= (~CHILD2367); 
            }
          }
          if (num_children & CHILD15) // do children 15 with small X and small Y exist
          {
            if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // does sphere intersect small Y
            {
              if (num_children & CHILD1) // does child 1 with small X and small Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[1], sphere, s_idx)) return true;
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123);
                }
              }
              if (num_children & CHILD5) // does child 5 with small X and small Y and small Z exist
              {
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new_notunrolled(cell->child[5], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with small Y from further consideration
            {
              num_children &= (~CHILD0145); 
            }
          }
        }
        else // remove children with small X from further consideration
        {
          num_children &= (~CHILD1357);
        }
      }
    }
    return false;
  }
  else
  {
    if (cell_overlaps_sphere(cell, sphere))
    {
      sphere->cell_idx = cell->idx;
      *tet_next = cell->data;
      cell->data = tet_idx;
#ifdef COLLECT_STATISTICS
      stat_box_overlap_true++;
#endif
      return true;
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_box_overlap_false++;
#endif
      return false;
    }
  }
}
*/

const static bool intersect_sphere_new(SScell3D* cell, SSsphere* sphere, int s_idx) 
{
  int num_children = cell->num_children;

  if (num_children)
  {
    if ((grid_cell_diagonal[cell->level] > sphere->rad) || cell_overlaps_sphere1(cell, sphere))
    {
      if (num_children & CHILD0246) // do children 0246 with large X exist
      {
        if ((sphere->cen[0] + sphere->rad) > cell->r_mid[0]) // does sphere intersect large X
        {
          if (num_children & CHILD26) // do children 26 with large X and large Y exist
          {
            if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // does sphere intersect large Y
            {
              if (num_children & CHILD2) // does child 2 with large X and large Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new(cell->child[2], sphere, s_idx)) return true;
                  if (num_children & CHILD6) // does child 6 with large X and large Y and small Z also exist
                  {
                    if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere also intersect small Z
                    {
                      if (intersect_sphere_new(cell->child[6], sphere, s_idx)) return true;
                    }
                    else // remove children with small Z from further consideration
                    {
                      num_children &= (~CHILD4567);
                    }
                  }
                }
                else // remove children with large Z from further consideration 
                {
                  num_children &= (~CHILD0123); // that means sphere must intersect small Z
                  if (num_children & CHILD6) // does child 6 with large X and large Y and small Z exist
                  {
                    if (intersect_sphere_new(cell->child[6], sphere, s_idx)) return true;
                  }
                }
              }
              else // if child 2 (with large Z) does not exist then child 6 (with small Z) must exist
              {
                assert(num_children & CHILD6);
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new(cell->child[6], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with large Y from further consideration
            {
              num_children &= (~CHILD2367); // that means sphere must intersect small Y
            }
          }
          if (num_children & CHILD04) // do children 04 with large X and small Y exist
          {
            if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // does sphere intersect small Y
            {
              if (num_children & CHILD0) // does child 0 with large X and small Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new(cell->child[0], sphere, s_idx)) return true;
                  if (num_children & CHILD4) // does child 4 with large X and small Y and small Z also exist
                  {
                    if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere also intersect small Z
                    {
                      if (intersect_sphere_new(cell->child[4], sphere, s_idx)) return true;
                    }
                    else // remove children with small Z from further consideration
                    {
                      num_children &= (~CHILD4567);
                    }
                  }
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123); // that means sphere must intersect small Z
                  if (num_children & CHILD4) // does child 4 with large X and small Y and small Z exist
                  {
                    if (intersect_sphere_new(cell->child[4], sphere, s_idx)) return true;
                  }
                }
              }
              else // if child 0 (with large Z) does not exist then child 4 (with small Z) must exist
              {
                assert(num_children & CHILD4);
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new(cell->child[4], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with small Y from further consideration
            {
              num_children &= (~CHILD0145); 
            }
          }
        }
        else // remove children with large X from further consideration
        {
          num_children &= (~CHILD0246);
        }
      }

      if (num_children & CHILD1357) // do children 1357 with small X exist
      {
        if ((sphere->cen[0] - sphere->rad) < cell->r_mid[0]) // does sphere intersect small X
        {
          if (num_children & CHILD37) // do children 37 with small X and large Y exist
          {
            if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // does sphere intersect large Y
            {
              if (num_children & CHILD3) // does child 3 with small X and large Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new(cell->child[3], sphere, s_idx)) return true;
                  if (num_children & CHILD7) // does child 7 with small X and large Y and small Z also exist
                  {
                    if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere also intersect small Z
                    {
                      if (intersect_sphere_new(cell->child[7], sphere, s_idx)) return true;
                    }
                    else // remove children with small Z from further consideration
                    {
                      num_children &= (~CHILD4567);
                    }
                  }
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123); // that means sphere must intersect small Z
                  if (num_children & CHILD7) // does child 7 with small X and large Y and small Z exist
                  {
                    if (intersect_sphere_new(cell->child[7], sphere, s_idx)) return true;
                  }
                }
              }
              else // if child 3 (with large Z) does not exist then child 7 (with small Z) must exist
              {
                assert(num_children & CHILD7);
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new(cell->child[7], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with large Y from further consideration
            {
              num_children &= (~CHILD2367); 
            }
          }
          if (num_children & CHILD15) // do children 15 with small X and small Y exist
          {
            if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // does sphere intersect small Y
            {
              if (num_children & CHILD1) // does child 1 with small X and small Y and large Z exist
              {
                if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // does sphere intersect large Z
                {
                  if (intersect_sphere_new(cell->child[1], sphere, s_idx)) return true;
                  if (num_children & CHILD5) // does child 5 with small X and small Y and small Z also exist
                  {
                    if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere also intersect small Z
                    {
                      if (intersect_sphere_new(cell->child[5], sphere, s_idx)) return true;
                    }
                    else // remove children with small Z from further consideration
                    {
                      num_children &= (~CHILD4567);
                    }
                  }
                }
                else // remove children with large Z from further consideration
                {
                  num_children &= (~CHILD0123); // that means sphere must intersect small Z
                  if (num_children & CHILD5) // does child 5 with small X and small Y and small Z exist
                  {
                    if (intersect_sphere_new(cell->child[5], sphere, s_idx)) return true;
                  }
                }
              }
              else // if child 1 (with large Z) does not exist then child 5 (with small Z) must exist
              {
                assert(num_children & CHILD5);
                if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // does sphere intersect small Z
                {
                  if (intersect_sphere_new(cell->child[5], sphere, s_idx)) return true;
                }
                else // remove children with small Z from further consideration
                {
                  num_children &= (~CHILD4567);
                }
              }
            }
            else // remove children with small Y from further consideration
            {
              num_children &= (~CHILD0145); 
            }
          }
        }
        else // remove children with small X from further consideration
        {
          num_children &= (~CHILD1357);
        }
      }
    }
    return false;
  }
  else
  {
    if (cell_overlaps_sphere(cell, sphere))
    {
      sphere->cell = cell;
      sphere->next = cell->data;
      cell->data = s_idx;
#ifdef COLLECT_STATISTICS
      stat_box_overlap_true++;
#endif
      return true;
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_box_overlap_false++;
#endif
      return false;
    }
  }
}

/*
const static bool intersect_sphere(SScell3D* cell, SSsphere* sphere, int* tet_next, int tet_idx) 
{
  if (cell->num_children) // check if cell has children
  {
    if ((sphere->cen[0] + sphere->rad) > cell->r_mid[0]) 
    {
      if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // do i need to check 2 and 6
      {
        if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // do i need to check 2
        {
          if (cell->child[2] && intersect_sphere(cell->child[2], sphere, s_idx)) return true;
        }
        if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // do i need to check 6
        {
          if (cell->child[6] && intersect_sphere(cell->child[6], sphere, s_idx)) return true;
        }
      }
      if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // do i need to check 0 and 4
      {
        if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // do i need to check 0
        {
          if (cell->child[0] && intersect_sphere(cell->child[0], sphere, s_idx)) return true;
        }
        if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // do i need to check 4
        {
          if (cell->child[4] && intersect_sphere(cell->child[4], sphere, s_idx)) return true;
        }
      }
    }
    if ((sphere->cen[0] - sphere->rad) < cell->r_mid[0]) // do i need to check 1, 3, 5, and 7
    {
      if ((sphere->cen[1] + sphere->rad) > cell->r_mid[1]) // do i need to check 3 and 7
      {
        if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // do i need to check 3
        {
          if (cell->child[3] && intersect_sphere(cell->child[3], sphere, s_idx)) return true;
        }
        if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // do i need to check 7
        {
          if (cell->child[7] && intersect_sphere(cell->child[7], sphere, s_idx)) return true;
        }
      }
      if ((sphere->cen[1] - sphere->rad) < cell->r_mid[1]) // do i need to check 1 and 5
      {
        if ((sphere->cen[2] + sphere->rad) > cell->r_mid[2]) // do i need to check 1
        {
          if (cell->child[1] && intersect_sphere(cell->child[1], sphere, s_idx)) return true;
        }
        if ((sphere->cen[2] - sphere->rad) < cell->r_mid[2]) // do i need to check 5
        {
          if (cell->child[5] && intersect_sphere(cell->child[5], sphere, s_idx)) return true;
        }
      }
    }
    return false;
  }
  else
  {
    if (cell_overlaps_sphere(cell, sphere))
    {
      sphere->cell_idx = cell->idx;
      *tet_next = cell->data;
      cell->data = tet_idx;
      return true;
    }
    else
    {
      return false;
    }
  }
};
*/

const int children_mask[] = {
  ((1<<7)|(1<<5)|(1<<4)|(1<<3)|(1<<1)),
  ((1<<5)),
  ((1<<7)|(1<<6)|(1<<5)|(1<<4)|(1<<3)|(1<<1)|(1<<0)),
  ((1<<7)|(1<<5)|(1<<1)),
  ((1<<7)|(1<<5)|(1<<3)|(1<<1)),
  ((0)),
  ((1<<7)|(1<<5)|(1<<4)|(1<<3)|(1<<1)|(1<<0)),
  ((1<<5)|(1<<1))
};

const static bool intersect_sphere_remaining(SScell3D* cell, SSsphere* sphere, int s_idx) 
{
  // loop over all ancestors of the finalized cell and continue the search
  while (cell->idx)
  {
    // get the parent cell
    SScell3D* parent = cell->parent;
    // let's see which child we were
    int i = cell->idx - (8 * parent->idx + 1);
    // get the remaining siblings after bistmasking the parent
    int remaining_siblings = parent->num_children & children_mask[i];
    // we only need to check the parent if there are remaining siblings
    if (remaining_siblings)
    {
      // temporarily store the children bistmask of the parent
      int temp = parent->num_children;
      // this temporarily "disables" me and all "if-dependent" siblings in the cell
      parent->num_children = remaining_siblings;
      // check the remaining siblings
      if (intersect_sphere_new(parent, sphere, s_idx))
      {
        // restore the num_children field
        parent->num_children = temp;
        return true;
      }
      // restore the num_children field
      parent->num_children = temp;
    }
    // move one level up
    cell = parent;
  }
  return false;
}

static float fin_cell_min_f[3];
static float fin_cell_max_f[3];

bool SScontainer3D::was_sphere_just_finalized(const SSsphere* sphere) const
{
#ifdef COLLECT_STATISTICS
  stat_finalize_justtested++;
#endif
  // for efficiency we should never use this function if the root is zero
  assert(root); 
  // does the sphere poke out of the cell in any direction 
  if (sphere->cen[0] - sphere->rad <= fin_cell_min_f[0]) return false;
  if (sphere->cen[0] + sphere->rad >= fin_cell_max_f[0]) return false;
  if (sphere->cen[1] - sphere->rad <= fin_cell_min_f[1]) return false;
  if (sphere->cen[1] + sphere->rad >= fin_cell_max_f[1]) return false;
  if (sphere->cen[2] - sphere->rad <= fin_cell_min_f[2]) return false;
  if (sphere->cen[2] + sphere->rad >= fin_cell_max_f[2]) return false;
#ifdef COLLECT_STATISTICS
  stat_finalize_justfinalized++;
#endif
  return true;
}

bool SScontainer3D::is_sphere_finalized_parent(SSsphere* sphere, int s_idx) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root); 
  // intersect the sphere with the sub-tree rooted at the parent of the "finalized" cell 
  if (intersect_sphere_new(finalized, sphere, s_idx))
  {
    return false;
  }
  // intersect the sphere with the remaining tree starting from the parent of the "finalized" cell 
  return (intersect_sphere_remaining(finalized, sphere, s_idx) == false);
}

bool SScontainer3D::is_sphere_finalized_leaf(SSsphere* sphere, int s_idx) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root); 
  // intersect the sphere with the remaining tree starting from the parent of the "finalized" leaf 
  return (intersect_sphere_remaining(finalized, sphere, s_idx) == false);
}

bool SScontainer3D::is_sphere_finalized(SSsphere* sphere, int s_idx) const
{
  // for efficiency we should never use this function if the root is zero
  assert(root); 
  // intersect the sphere with the finalization tree starting from the root in DFS order
  return (intersect_sphere_new(root, sphere, s_idx) == false);
}

// this function deallocates a subtree. it starts at cell. deallocation
// includes that cell. it returns a linked list of leaf cells (but only
// those actually contain data) that are not deallocated yet.
static SScell3D* finalize_subtree(SScell3D* cell, SScell3D* finalized)
{
  // is this cell a parent cell
  if (cell->num_children)
  {
    // parents cell that are deallocated should not have data
    assert(cell->data == -1);
    // finalize all existing children
    for (int i = 0; i < 8; i++)
    {
      if (cell->child[i])
      {
        cell_hash->erase(cell->child[i]->idx);
        finalized = finalize_subtree(cell->child[i], finalized);
      }
    }
    // all parents of this cell are finalized
    cell->num_children = 0;
    // it is safe to deallocate this parent cell already
    deallocCell(cell);
    return finalized;
  }
  else
  {
    // the leaf cell may have data
    if (cell->data != -1)
    {
      // in which case we link it into finalized (but no dealloc)
      cell->buffer_next = finalized;
      return cell;
    }
    else
    {
      // it is safe to deallocate this leaf cell already
      deallocCell(cell);
      return finalized;
    }
  }
}

// this function creates the subtree for the leaf idx starting bottom
// up all the way to an existing parent. for example, the very first
// finalize call (when only the root node exists) creates the unfinalized
// ancestors of the leaf, including the seven sibling leaves. it does
// not create the finalized cell itself.
// that pre-existing parent of the finalized cell is only stored in the 
// field "finalized_parent" if it contains any data.
SScell3D* SScontainer3D::create_parent_and_siblings(int idx, bool finalize)
{
  SScell3D* parent;
  int parent_idx = (idx - 1) / 8;
  my_cell_hash::iterator hash_element = cell_hash->find(parent_idx);

  if (hash_element == cell_hash->end())
  {
    // this could not have been the root node 
    assert(parent_idx > 0);
    // did not find parent ... create that parent, its siblings, and their parent
    parent = create_parent_and_siblings(parent_idx, false);
  }
  else
  {
    // pre-existing parent was found
    parent = (*hash_element).second;
    // does this parent have data
    if (parent->data != -1)
    {
      finalized_parent = parent;
    }
    else
    {
      finalized_parent = 0;
    }
  }

  // just make sure
  assert(parent_idx == parent->idx);

  // compute mid point of parent cell
  parent->r_mid[0] = (parent->r_min[0] + parent->r_max[0])/2;
  parent->r_mid[1] = (parent->r_min[1] + parent->r_max[1])/2;
  parent->r_mid[2] = (parent->r_min[2] + parent->r_max[2])/2;

  // create siblings
  SScell3D* sibling;
  SScell3D* return_parent = 0;
  int sibling_idx = 8*parent_idx;
  int sibling_level = parent->level + 1;

  for (int i = 0; i < 8; i++)
  {
    // compute index of siblings
    sibling_idx++;
    // is it a parental or an actual sibling (but not the cell itself)
    if (sibling_idx == idx)
    {
      if (finalize)
      {
#ifdef PRINT_DEBUG_OUTPUT
//        fprintf(stderr, "finalize %d (%g,%g,%g)/(%g,%g,%g)\n", idx, (i&1?parent->r_min[0]:parent->r_mid[0]), (i&2?parent->r_mid[1]:parent->r_min[1]), (i&4?parent->r_min[2]:parent->r_mid[2]), (i&1?parent->r_mid[0]:parent->r_max[0]), (i&2?parent->r_max[1]:parent->r_mid[1]), (i&4?parent->r_mid[2]:parent->r_max[2]));
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
        if (i&4)
        {
          fin_cell_min_f[2] = parent->r_min[2];
          fin_cell_max_f[2] = parent->r_mid[2];
        }
        else
        {
          fin_cell_min_f[2] = parent->r_mid[2];
          fin_cell_max_f[2] = parent->r_max[2];
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
    // give sibling its level
    sibling->level = sibling_level;
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
    if (i&4)
    {
      sibling->r_min[2] = parent->r_min[2];
      sibling->r_max[2] = parent->r_mid[2];
    }
    else
    {
      sibling->r_min[2] = parent->r_mid[2];
      sibling->r_max[2] = parent->r_max[2];
    }
    // insert it into the hash
    cell_hash->insert(my_cell_hash::value_type(sibling_idx, sibling));
    // attach sibling to parent
    parent->child[i] = sibling;
    parent->num_children |= (1<<i);
    // make sibling point to parent
    sibling->parent = parent;
  }

  // check if we need to compute the cell_grid_diagonal for this level
  if (deepest_level < sibling_level)
  {
    assert(deepest_level + 1 == sibling_level);
    deepest_level = sibling_level;
    assert(grid_cell_diagonal[deepest_level] == -1);
    grid_cell_diagonal[deepest_level] = VecDistance3fv(sibling->r_min, sibling->r_max);
  }

  assert(finalize || return_parent);
  return return_parent;
}

int SScontainer3D::finalize_cell(int cell_idx)
{
  my_cell_hash::iterator hash_element;
  
  hash_element = cell_hash->find(cell_idx);

  // only for debugging
  r_min_f = fin_cell_min_f;
  r_max_f = fin_cell_max_f;

  if (hash_element == cell_hash->end())
  {
//    fprintf(stderr, "finalize non-existing cell %d \n",cell_idx);
    // did not find cell ... create the subtree surrounding it
    create_parent_and_siblings(cell_idx, true);
    // this will place the parent of the cell in the finalized
    // field (but only in case it actually has data)
    if (finalized_parent)
    {
      return 1; // means the finalized cell is a parent
    }
  }
  else
  {
//    fprintf(stderr, "finalize existing cell %d \n",cell_idx);
    // found cell
    SScell3D* cell = (*hash_element).second;
    // make sure this is true
    assert(cell->idx == cell_idx);
    // remove from hash
    cell_hash->erase(hash_element);
    // get parent of cell
    SScell3D* parent = cell->parent;
    // does parent exist?
    if (parent)
    {
      // which of its children is it
      int i = cell_idx-(8*parent->idx+1);
      // make sure this child exists
      assert(parent->child[i]);
      // it will exist no longer
      parent->child[i] = 0;
      parent->num_children &= (~(1<<i));
      // we must still have children
      assert(parent->num_children);
    }
    else
    {
      // this must be the root cell
      assert(cell == root);
      // which will now longer exist
      root = 0;
    }
    // here we could check if any immediate neighbors
    // of the cell are already finalized and then expand
    // the size of this finalized region
//  code goes here for doing this 
    // copy the bounding box of the finalized cell
    fin_cell_min_f[0] = cell->r_min[0];
    fin_cell_min_f[1] = cell->r_min[1];
    fin_cell_min_f[2] = cell->r_min[2];
    fin_cell_max_f[0] = cell->r_max[0];
    fin_cell_max_f[1] = cell->r_max[1];
    fin_cell_max_f[2] = cell->r_max[2];
    // finalize the subtree rooted in cell
    finalized_leaf = finalize_subtree(cell, 0);
    if (finalized_leaf)
    {
      return 2; // means the finalized cell is a leaf
    }
  }
  // means the finalized cell had no data
  return 0;
}

bool SScontainer3D::write_data_to_cell(int idx, int data)
{
  my_cell_hash::iterator hash_element = cell_hash->find(idx);

  if (hash_element == cell_hash->end())
  {
    return false;
  }
  ((*hash_element).second)->data = data;
  return true;
}

bool SScontainer3D::get_data_from_cell(int idx)
{
  my_cell_hash::iterator hash_element = cell_hash->find(idx);

  if (hash_element == cell_hash->end())
  {
    return false;
  }
  data = ((*hash_element).second)->data;
  return true;
}

int SScontainer3D::prepareParent()
{
  assert(finalized_parent);
  finalized = finalized_parent;
  assert(finalized->buffer_next == 0);
  int data = finalized->data;
  finalized->data = -1;
  return data;
}

int SScontainer3D::prepareLeaf()
{
  if (finalized_leaf)
  {
    finalized = finalized_leaf;
    assert(finalized->num_children == 0);
    finalized_leaf = finalized->buffer_next;
    int data = finalized->data;
    deallocCell(finalized);
    return data;
  }
  else
  {
    return -1;
  }
}

void SScontainer3D::iterateInit()
{
  hash_iterator = cell_hash->begin();
}

bool SScontainer3D::iterateNext()
{
  SScell3D* cell;

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

bool SScontainer3D::noiterate(int cell_idx)
{
  my_cell_hash::iterator hash_element = cell_hash->find(cell_idx);
  
  if (hash_element == cell_hash->end())
  {
    return false;
  }
  else
  {
    SScell3D* cell = (*hash_element).second;
    if (cell->num_children == 0)
    {
      r_min_f = cell->r_min;
      r_max_f = cell->r_max;
      return true;
    }
    else
    {
      return false;
    }
  }
}

void SScontainer3D::open(const float* bb_min_f, const float* bb_max_f)
{
  // this is a hack to provide access for a local function
  this_container = this;
  this->bb_min_f[0] = bb_min_f[0];
  this->bb_min_f[1] = bb_min_f[1];
  this->bb_min_f[2] = bb_min_f[2];
  this->bb_max_f[0] = bb_max_f[0];
  this->bb_max_f[1] = bb_max_f[1];
  this->bb_max_f[2] = bb_max_f[2];

  root = allocCell();
  root->idx = 0;
  root->level = 0;
  root->r_min[0] = this->bb_min_f[0];
  root->r_min[1] = this->bb_min_f[1];
  root->r_min[2] = this->bb_min_f[2];
  root->r_max[0] = this->bb_max_f[0];
  root->r_max[1] = this->bb_max_f[1];
  root->r_max[2] = this->bb_max_f[2];
  grid_cell_diagonal[0] = VecDistance3fv(bb_min_f, bb_max_f);
  deepest_level = 0;
  for (int i = 1; i < 16; i++) grid_cell_diagonal[i] = -1;
  cell_hash->insert(my_cell_hash::value_type(0, root));
}

void SScontainer3D::open(const double* bb_min_d, const double* bb_max_d)
{
  // this is a hack to provide access for a local function
  this_container = this;
  this->bb_min_f[0] = (float)bb_min_d[0];
  this->bb_min_f[1] = (float)bb_min_d[1];
  this->bb_min_f[2] = (float)bb_min_d[2];
  this->bb_max_f[0] = (float)bb_max_d[0];
  this->bb_max_f[1] = (float)bb_max_d[1];
  this->bb_max_f[2] = (float)bb_max_d[2];

  // make sure the min values are *small* enough
  while (bb_min_d[0] < this->bb_min_f[0])
  {
    fprintf(stderr,"WARNING: making bb_min_f[0] %g smaller than %f\n",this->bb_min_f[0],bb_min_d[0]);
    this->bb_min_f[0] -= single_epsilon;
  }
  while (bb_min_d[1] < this->bb_min_f[1])
  {
    fprintf(stderr,"WARNING: making bb_min_f[1] %g smaller than %f\n",this->bb_min_f[1],bb_min_d[1]);
    this->bb_min_f[1] -= single_epsilon;
  }
  while (bb_min_d[2] < this->bb_min_f[2])
  {
    fprintf(stderr,"WARNING: making bb_min_f[2] %g smaller than %f\n",this->bb_min_f[2],bb_min_d[2]);
    this->bb_min_f[2] -= single_epsilon;
  }
  // make sure the max values are *big* enough
  while (bb_max_d[0] > this->bb_max_f[0])
  {
    fprintf(stderr,"WARNING: making bb_max_f[0] %g bigger than %f\n",this->bb_max_f[0],bb_max_d[0]);
    this->bb_max_f[0] += single_epsilon;
  }
  while (bb_max_d[1] > this->bb_max_f[1])
  {
    fprintf(stderr,"WARNING: making bb_max_f[1] %g bigger than %f\n",this->bb_max_f[1],bb_max_d[1]);
    this->bb_max_f[1] += single_epsilon;
  }
  while (bb_max_d[2] > this->bb_max_f[2])
  {
    fprintf(stderr,"WARNING: making bb_max_f[2] %g bigger than %f\n",this->bb_max_f[2],bb_max_d[2]);
    this->bb_max_f[2] += single_epsilon;
  }

  root = allocCell();
  root->idx = 0;
  root->level = 0;
  root->r_min[0] = this->bb_min_f[0];
  root->r_min[1] = this->bb_min_f[1];
  root->r_min[2] = this->bb_min_f[2];
  root->r_max[0] = this->bb_max_f[0];
  root->r_max[1] = this->bb_max_f[1];
  root->r_max[2] = this->bb_max_f[2];
  grid_cell_diagonal[0] = VecDistance3fv(bb_min_f, bb_max_f);
  deepest_level = 0;
  for (int i = 1; i < 16; i++) grid_cell_diagonal[i] = -1;
  cell_hash->insert(my_cell_hash::value_type(0, root));
}

void SScontainer3D::close()
{
  if (root) fprintf(stderr, "WARNING: closing container containing unfinalized\n");
  cell_hash->clear();
}

SScontainer3D::SScontainer3D()
{
  cell_hash = new my_cell_hash;
  root = 0;
  finalized = 0;
  finalized_parent = 0;
  finalized_leaf = 0;

  // set hard-coded epsilons
	single_epsilon = 1.0;
	{ for (int i=0; i<23; i++){ single_epsilon*=0.5;}}
	double_epsilon = 1.0;
	{ for (int i=0; i<52; i++){ double_epsilon*=0.5;}}
	sqrt_epsilon = single_epsilon;

  // compute level offsets
  level_offset[0] = 0;
  for (int l = 0; l < 19; l++) level_offset[l+1] = level_offset[l] + ((1<<l)*(1<<l)*(1<<l));
}

SScontainer3D::~SScontainer3D()
{
#ifdef COLLECT_STATISTICS
  fprintf(stderr, "sphere_initialized %d\n", stat_sphere_initialized);
  fprintf(stderr, "justtested %d justfinalized %d fulltested %d\n", stat_finalize_justtested, stat_finalize_justfinalized, stat_finalize_justtested-stat_finalize_justfinalized);
  fprintf(stderr, "stat_box_overlap: true %d false %d\n", stat_box_overlap_true, stat_box_overlap_false);
  fprintf(stderr, "stat_box_overlap1: true %d false %d\n", stat_box_overlap1_true, stat_box_overlap1_false);
#endif
  delete cell_hash;
}
