/*
===============================================================================

  FILE:  spdelaunay2d.cpp
  
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
#include "spdelaunay2d.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "vec3fv.h"

#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

#ifdef COLLECT_STATISTICS
static int stat_finalize_calls = 0;
static int stat_finalize_active_checked = 0;
static int stat_finalize_active_infinite = 0;
static int stat_finalize_active_circled = 0;
static int stat_finalize_active_notcircled = 0;
#endif


/*
static bool isFinalized(Delaunay2Vertex* v0, Delaunay2Vertex* v1, SScontainer2D* ss2d)
{
  float p_nor[2];
  p_nor[0] = (float)(v1->x[1] - v0->x[1]);
  p_nor[1] = (float)(v0->x[0] - v1->x[0]);
  float p_pnt[2];
  p_pnt[0] = (float)(p_nor[0]*v1->x[0] + p_nor[1]*v1->x[1]) - 0.5f;
  return ss2d->is_finalized(p_nor, p_pnt);
}
*/

int area_sign(Delaunay2Vertex* a, Delaunay2Vertex* b, Delaunay2Vertex* c);

static bool isFinalized(Delaunay2Triangle* t, SScontainer2D* ss2d)
{
	assert( area_sign(t->V[0],t->V[1],t->V[2]) );
	
  if (t->rad < 0) // does not yet have a circle
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_circled++;
#endif
    t->initialize(t->V[0]->x,t->V[1]->x,t->V[2]->x);
    return ss2d->is_finalized(t);
  }
  else
  {
#ifdef COLLECT_STATISTICS
    stat_finalize_active_notcircled++;
#endif
    return ss2d->is_still_finalized(t);
  }
}

static void CheckFinalizeAndWrite(Delaunay2* dt, SScontainer2D* ss2d, SMwriter* smwriter)
{
  Delaunay2Triangle* t;
  int act_t_num = dt->act_t_num;
  int idx = dt->youngestActiveTriangle();

#ifdef COLLECT_STATISTICS
  stat_finalize_calls++;
  stat_finalize_active_checked += act_t_num;
#endif

	while (act_t_num)
  {
    t = dt->triangle_buffer + idx;
    act_t_num--;

    // is this an infinite triangle with v0 being infinite vertex
    if ( t->V[0] == dt->pinf ) 
    {
      // and is its finite edge already finalized
      if (t->N[0] == D2_NULL_INDEX)
      {
        // and it's other two vertices have the x coordinate in common
        if (t->V[1]->x[0] == t->V[2]->x[0])
        {
          // and this coordinate is on the bounding box
          if (t->V[1]->x[0] == ss2d->bb_min_f[0] || t->V[1]->x[0] == ss2d->bb_max_f[0])
          {
            // then you can finalize the infinite triangle
            for (int j = 0; j < 3; j++)
            {
              Delaunay2Vertex* v = t->V[j];
              v->use_count--;
              if (v->use_count == 0)
              {
  			        dt->deallocDelaunay2Vertex(v);	
              }
            }
	          //make all its neighbor corners pointing to null
            if (t->N[0]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[0]))->N[CIND(t->N[0])] = D2_NULL_INDEX;
            if (t->N[1]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[1]))->N[CIND(t->N[1])] = D2_NULL_INDEX;
            if (t->N[2]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[2]))->N[CIND(t->N[2])] = D2_NULL_INDEX;
 
            // take it out of the active triangle list
            dt->delActiveTriangle(idx);

            // free its memory
            dt->deallocDelaunay2Triangle(idx);
          }
        }
        // or else the y coordinate in common
        else if (t->V[1]->x[1] == t->V[2]->x[1])
        {
          // and this coordinate is on the bounding box
          if (t->V[1]->x[1] == ss2d->bb_min_f[1] || t->V[1]->x[1] == ss2d->bb_max_f[1])
          {
            // then you can finalize the infinite triangle
            for (int j = 0; j < 3; j++)
            {
              Delaunay2Vertex* v = t->V[j];
              v->use_count--;
              if (v->use_count == 0)
              {
  			        dt->deallocDelaunay2Vertex(v);	
              }
            }
	          //make all its neighbor corners pointing to null
            if (t->N[1]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[1]))->N[CIND(t->N[1])] = D2_NULL_INDEX;
            if (t->N[2]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[2]))->N[CIND(t->N[2])] = D2_NULL_INDEX;
 
            // take it out of the active triangle list
            dt->delActiveTriangle(idx);

            // free its memory
            dt->deallocDelaunay2Triangle(idx);
          }
        }
      }
    }
    else if ( t->V[1] == dt->pinf )
    {
      // and is its finite edge already finalized
      if (t->N[1] == D2_NULL_INDEX)
      {
        // and it's other two vertices have the x coordinate in common
        if (t->V[0]->x[0] == t->V[2]->x[0])
        {
          // and this coordinate is on the bounding box
          if (t->V[0]->x[0] == ss2d->bb_min_f[0] || t->V[0]->x[0] == ss2d->bb_max_f[0])
          {
            // then you can finalize the infinite triangle
            for (int j = 0; j < 3; j++)
            {
              Delaunay2Vertex* v = t->V[j];
              v->use_count--;
              if (v->use_count == 0)
              {
  			        dt->deallocDelaunay2Vertex(v);	
              }
            }
	          //make all its neighbor corners pointing to null
            if (t->N[0]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[0]))->N[CIND(t->N[0])] = D2_NULL_INDEX;
            if (t->N[2]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[2]))->N[CIND(t->N[2])] = D2_NULL_INDEX;
 
            // take it out of the active triangle list
            dt->delActiveTriangle(idx);

            // free its memory
            dt->deallocDelaunay2Triangle(idx);
          }
        }
        // or else the y coordinate in common
        else if (t->V[0]->x[1] == t->V[2]->x[1])
        {
          // and this coordinate is on the bounding box
          if (t->V[0]->x[1] == ss2d->bb_min_f[1] || t->V[0]->x[1] == ss2d->bb_max_f[1])
          {
            // then you can finalize the infinite triangle
            for (int j = 0; j < 3; j++)
            {
              Delaunay2Vertex* v = t->V[j];
              v->use_count--;
              if (v->use_count == 0)
              {
  			        dt->deallocDelaunay2Vertex(v);	
              }
            }
	          //make all its neighbor corners pointing to null
            if (t->N[0]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[0]))->N[CIND(t->N[0])] = D2_NULL_INDEX;
            if (t->N[2]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[2]))->N[CIND(t->N[2])] = D2_NULL_INDEX;
 
            // take it out of the active triangle list
            dt->delActiveTriangle(idx);

            // free its memory
            dt->deallocDelaunay2Triangle(idx);
          }
        }
      }
    }
    else if (t->V[0] != dt->pinf && t->V[1] != dt->pinf)
    {
      assert (t->V[2] != dt->pinf);
      // here we check whether the triangle's circle is completely inside the finalized area or
      // rather whether the active triangle's circles still intersects some non-finalized space
      if (isFinalized(t, ss2d))
      {
	      int t_idx[3];
	      bool t_final[3];
        float v_pos[3];
        for (int j = 0; j < 3; j++)
        {
          Delaunay2Vertex* v = t->V[j];
          if (v->index == -1)
          {
            v->index = smwriter->v_count;
            v_pos[0] = (float)v->x[0];
            v_pos[1] = (float)v->x[1];
            v_pos[2] = v->h;
            smwriter->write_vertex(v_pos);
          }
          t_idx[j] = v->index;
          v->use_count--;
          if (v->use_count == 0)
          {
            t_final[j] = true;
  			    dt->deallocDelaunay2Vertex(v);	
          }
          else
          {
            t_final[j] = false;
            assert(v->use_count > 0);
          }
        }
        smwriter->write_triangle(t_idx, t_final);

//        if (output) printf("%d: %f %f %f\n",idx,t->c_cen[0],t->c_cen[1],t->c_rad);

	      //make all its neighbor corners pointing to null
        if (t->N[0]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[0]))->N[CIND(t->N[0])] = D2_NULL_INDEX;
        if (t->N[1]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[1]))->N[CIND(t->N[1])] = D2_NULL_INDEX;
        if (t->N[2]!=D2_NULL_INDEX) (dt->triangle_buffer + CTRI(t->N[2]))->N[CIND(t->N[2])] = D2_NULL_INDEX;
 
        // take it out of the active triangle list
        dt->delActiveTriangle(idx);

        // free its memory
        dt->deallocDelaunay2Triangle(idx);
      }
      else
      {
#ifdef COLLECT_STATISTICS
        // check how often we could have a cheaper reject.
#endif
      }
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_active_infinite++;
#endif
    }
    idx = dt->nextYoungestActiveTriangle(idx);
  }
}

static void AllFinalizeAndWrite(Delaunay2* dt, SMwriter* smwriter)
{
  Delaunay2Triangle* t;
  int act_t_num = dt->act_t_num;
  int idx = dt->youngestActiveTriangle();

#ifdef COLLECT_STATISTICS
  stat_finalize_calls++;
  stat_finalize_active_checked += act_t_num;
#endif

	while (act_t_num)
  {
    t = dt->triangle_buffer + idx;
    act_t_num--;

    if (t->V[0] != dt->pinf && t->V[1] != dt->pinf)
    {
      // if the triangle is finite we can finalize and write it
      assert (t->V[2] != dt->pinf);
	    int t_idx[3];
	    bool t_final[3];
      float v_pos[3];
      for (int j = 0; j < 3; j++)
      {
        Delaunay2Vertex* v = t->V[j];
        if (v->index == -1)
        {
          v->index = smwriter->v_count;
          v_pos[0] = (float)v->x[0];
          v_pos[1] = (float)v->x[1];
          v_pos[2] = v->h;
          smwriter->write_vertex(v_pos);
        }
        t_idx[j] = v->index;
        v->use_count--;
        if (v->use_count == 0)
        {
          t_final[j] = true;
  			  dt->deallocDelaunay2Vertex(v);	
        }
        else
        {
          t_final[j] = false;
          assert(v->use_count > 0);
        }
      }
      smwriter->write_triangle(t_idx, t_final);
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_active_infinite++;
#endif
    }

    // take them all out of the active triangle list
    dt->delActiveTriangle(idx);

    // free their memory
    dt->deallocDelaunay2Triangle(idx);

    // get the next one
    idx = dt->nextYoungestActiveTriangle(idx);
  }
  assert(act_t_num == 0);
}

bool SPdelaunay2D::open(SMwriter* smwriter)
{
  this->smwriter = smwriter;
  this->ss2d = new SScontainer2D();
  this->dt = new Delaunay2();
  p_count = 0;
  p_count_init = 0;
  return true;
}

void SPdelaunay2D::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  ss2d->open(bb_min_f, bb_max_f);

  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];

  VecCopy3fv(this->bb_min_f, bb_min_f);
  VecCopy3fv(this->bb_max_f, bb_max_f);

  smwriter->set_boundingbox(this->bb_min_f, this->bb_max_f);
}

void SPdelaunay2D::set_boundingbox(const double* bb_min_d, const double* bb_max_d)
{
  ss2d->open(bb_min_d, bb_max_d);

  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];

  VecCopy3fv(this->bb_min_f, bb_min_d);
  VecCopy3fv(this->bb_max_f, bb_max_d);

  smwriter->set_boundingbox(this->bb_min_f, this->bb_max_f);
}

void SPdelaunay2D::set_finalizemethod(SPfinalizemethod finalizemethod)
{
}

static Delaunay2Vertex* p012[3];

void SPdelaunay2D::write_point(const float* p_pos_f)
{
  Delaunay2Vertex* p = dt->allocDelaunay2Vertex(p_pos_f);
  if (p_count_init == 3)
  {
  	dt->insert(p);
  }
  else
  {
    p012[p_count_init] = p;

    if (p_count_init == 2)
    { 
      if (area_sign(p012[0], p012[1], p012[2]) == 0)
      {
		    fprintf(stderr,"WARNING: first three points are in degenerate position.\n");
        exit(0);
      }
  	  dt->initialize(p012[0], p012[1], p012[2]);
    }
    p_count_init++;
  }
  p_count++;
}

void SPdelaunay2D::write_point(const double* p_pos_d)
{
  Delaunay2Vertex* p = dt->allocDelaunay2Vertex(p_pos_d);
  if (p_count_init == 3)
  {
  	dt->insert(p);
  }
  else
  {
    p012[p_count_init] = p;

    if (p_count_init == 2)
    { 
      if (area_sign(p012[0], p012[1], p012[2]) == 0)
      {
		    fprintf(stderr,"WARNING: first three points are in degenerate position.\n");
        exit(0);
      }
  	  dt->initialize(p012[0], p012[1], p012[2]);
    }
    p_count_init++;
  }
  p_count++;
}

void SPdelaunay2D::write_finalize_cell(int cell_idx)
{
  ss2d->finalize_cell(cell_idx);
  if (cell_idx != 0)
  {
    CheckFinalizeAndWrite(dt, ss2d, smwriter);
  }
  else
  {
    AllFinalizeAndWrite(dt, smwriter);
  }
#ifndef NDEBUG		
	dt->audit();
#endif
}

void SPdelaunay2D::close()
{
  // this finalizes whatever space is remaining
  ss2d->close();
  // write whatever is remaining
//  AllFinalizeAndWrite(dt, smwriter);
#ifndef NDEBUG		
  dt->audit();
#endif

//	fprintf(stderr,"max number of allocated quad tree cells: %d\n", ss2d->triangle_buffer_maxsize);
  // delete spatial finalizer

#ifdef COLLECT_STATISTICS
  fprintf(stderr, "stat_finalize_calls %d\n", stat_finalize_calls);
  fprintf(stderr, "stat_finalize_active_checked %d\n", stat_finalize_active_checked);
  fprintf(stderr, "infinite %d circled %d notcircled %d\n", stat_finalize_active_infinite, stat_finalize_active_circled, stat_finalize_active_notcircled);
#endif
  delete ss2d;

  fprintf(stderr,"max number of used vertices/triangles: %d/%d\n", dt->vertex_buffer_maxsize, dt->triangle_buffer_maxsize, dt->vertex_buffer_size, dt->triangle_buffer_size);

  delete dt;
}

SPdelaunay2D::SPdelaunay2D()
{
  // init of SPwriter interface
  ncomments = 0;
  comments = 0;

  npoints = -1;
  p_count = -1;

  bb_min_f = 0;
  bb_max_f = 0;

  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPdelaunay2D interface
  smwriter = 0;
  ss2d = 0;
  dt = 0;
}

void SPdelaunay2D::getActiveTrianglesInit()
{
  active_triangle_next = 0;
}

int SPdelaunay2D::getActiveTrianglesNext(float* v0, float* v1, float* v2)
{
  while (active_triangle_next < dt->triangle_buffer_maxsize)
  {    
    Delaunay2Triangle* t = dt->triangle_buffer + active_triangle_next;

    if (t->dead)
    {
      active_triangle_next++;
      continue;
    }

    if (t->V[0] == dt->pinf || t->V[1] == dt->pinf)
    {
      active_triangle_next++;
      continue;
    }

    v0[0] = (float)t->V[0]->x[0];
    v0[1] = (float)t->V[0]->x[1];
    v0[2] = t->V[0]->h;;

    v1[0] = (float)t->V[1]->x[0];
    v1[1] = (float)t->V[1]->x[1];
    v1[2] = t->V[1]->h;

    v2[0] = (float)t->V[2]->x[0];
    v2[1] = (float)t->V[2]->x[1];
    v2[2] = t->V[2]->h;

    return active_triangle_next++;
  }
  return -1;
}

bool SPdelaunay2D::getActiveTrianglesIdx(int idx, float* v0, float* v1, float* v2)
{
  Delaunay2Triangle* t = dt->triangle_buffer + idx;

  if (t->dead)
  {
    return false;
  }

  if (t->V[0] == dt->pinf || t->V[1] == dt->pinf)
  {
    return false;
  }

  v0[0] = (float)t->V[0]->x[0];
  v0[1] = (float)t->V[0]->x[1];
  v0[2] = (float)t->V[0]->h;
  v1[0] = (float)t->V[1]->x[0];
  v1[1] = (float)t->V[1]->x[1];
  v1[2] = (float)t->V[1]->h;
  v2[0] = (float)t->V[2]->x[0];
  v2[1] = (float)t->V[2]->x[1];
  v2[2] = (float)t->V[2]->h;

  return true;
}

bool  SPdelaunay2D::getActiveTrianglesIdx(int idx, float* cen_rad)
{
  Delaunay2Triangle* t = dt->triangle_buffer + idx;

  if (t->dead)
  {
    return false;
  }

  if (t->V[0] == dt->pinf || t->V[1] == dt->pinf)
  {
    return false;
  }

  if (t->rad < 0)
  {
    fprintf(stderr, "WARNING: circumsphere of triangle %d not computed yet\n", idx);
    return false;
  }

  cen_rad[0] = t->cen[0];
  cen_rad[1] = t->cen[1];
  cen_rad[2] = (float)(t->V[0]->h+t->V[1]->h+t->V[0]->h)/3;
  cen_rad[3] = t->rad;

  return true;
}

void SPdelaunay2D::getInfiniteTrianglesInit()
{
  infinite_triangle_next = 0;
}

bool SPdelaunay2D::getInfiniteTrianglesNext(float* v0, float* v1, float* v2)
{
  while (infinite_triangle_next < dt->triangle_buffer_maxsize)
  {
		Delaunay2Triangle* t = dt->triangle_buffer + infinite_triangle_next;
    
    infinite_triangle_next++;

    if (!t->dead)
    {
      if (t->V[0] == dt->pinf)
      {
        float p_nor[2];
        p_nor[0] = (float)(t->V[2]->x[1] - t->V[1]->x[1]);
        p_nor[1] = (float)(t->V[1]->x[0] - t->V[2]->x[0]);

        if (p_nor[0] > 0) v0[0] = 2*bb_min_f[0] - bb_max_f[0];
        else v0[0] = 2*bb_max_f[0] - bb_min_f[0];

        if (p_nor[1] > 0) v0[1] = 2*bb_min_f[1] - bb_max_f[1];
        else v0[1] = 2*bb_max_f[1] - bb_min_f[1];
        
        v0[2] = bb_min_f[2];

        v1[0] = (float)t->V[1]->x[0];
        v1[1] = (float)t->V[1]->x[1];
        v1[2] = t->V[1]->h;

        v2[0] = (float)t->V[2]->x[0];
        v2[1] = (float)t->V[2]->x[1];
        v2[2] = t->V[2]->h;

        return true;
      }
      else if (t->V[1] == dt->pinf)
      {
        v0[0] = (float)t->V[0]->x[0];
        v0[1] = (float)t->V[0]->x[1];
        v0[2] = t->V[0]->h;;

        float p_nor[2];
        p_nor[0] = (float)(t->V[0]->x[1] - t->V[2]->x[1]);
        p_nor[1] = (float)(t->V[2]->x[0] - t->V[0]->x[0]);

        if (p_nor[0] > 0) v1[0] = 2*bb_min_f[0] - bb_max_f[0];
        else v1[0] = 2*bb_max_f[0] - bb_min_f[0];

        if (p_nor[1] > 0) v1[1] = 2*bb_min_f[1] - bb_max_f[1];
        else v1[1] = 2*bb_max_f[1] - bb_min_f[1];
        
        v1[2] = bb_min_f[2];

        v2[0] = (float)t->V[2]->x[0];
        v2[1] = (float)t->V[2]->x[1];
        v2[2] = t->V[2]->h;

        return true;
      }
    }
  }
  return false;
}
