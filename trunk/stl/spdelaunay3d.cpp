/*
===============================================================================

  FILE:  spdelaunay3d.cpp
  
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
#include "spdelaunay3d.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "vec3fv.h"

#define COLLECT_STATISTICS
#undef COLLECT_STATISTICS

//#define FINALIZE_CHECK_WITHOUT_SQRT

int stat_delayed_deleted_tetrahedra = 0;

#ifdef COLLECT_STATISTICS
static int stat_finalize_cell = 0;
static int stat_finalize_parent_cell = 0;
static int stat_finalize_leaf_cell = 0;
static int stat_finalize_cell_dead = 0;
static int stat_finalize_cell_checked = 0;
static int stat_finalize_cell_relinked = 0;
static int stat_finalize_new = 0;
static int stat_finalize_new_dead = 0;
static int stat_finalize_new_infinite = 0;
static int stat_finalize_new_infinite_finalized = 0;
static int stat_finalize_new_with_cell = 0;
static int stat_finalize_new_checked = 0;
static int stat_finalize_new_linked = 0;
static int stat_finalize_first = 0;
static int stat_finalize_again_parent = 0;
static int stat_finalize_again_leaf = 0;
#endif

// for debug only

static bool check_octree(Delaunay3* dt, SScontainer3D* ss3d)
{
  for (int idx = 0; idx < 256; idx++)
  {
    if (ss3d->get_data_from_cell(idx))
    {
      int t_idx = ss3d->data;
      while (t_idx != -1)
      {
        assert(dt->tpool[t_idx].cell->idx == idx);
        t_idx = dt->tpool[t_idx].next;
      }
    }
  }
  return true;
}

// this function writes out a tetrahedron

static void finalize_tetrahedron(int tet_idx, Delaunay3* dt, SVwriter* svwriter)
{
  int t_idx[4];
  bool t_final[4];
  float v_pos[3];
  Delaunay3Vertex* v;
  for (int j = 0; j < 4; j++)
  {
    v = dt->tpool[tet_idx].V[j];
    if (v->index == -1)
    {
      v->index = svwriter->v_count;
      v_pos[0] = (float)v->x[0];
      v_pos[1] = (float)v->x[1];
      v_pos[2] = (float)v->x[2];
      svwriter->write_vertex(v_pos);
    }
    t_idx[j] = v->index;
    v->ref_count--;
    if (v->ref_count == 0)
    {
      t_final[j] = true;
  	  dt->deallocVertex(v);	
    }
    else
    {
      t_final[j] = false;
      assert(v->ref_count > 0);
    }
  }
  svwriter->write_tetrahedron(t_idx, t_final);

  // make all its neighbor corners pointing to null

  int t1, i1; //indices of tetra and corner position
  if (dt->tpool[tet_idx].N[0]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[0]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[0]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[1]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[1]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[1]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[2]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[2]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[2]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[3]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[3]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[3]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }

  // free its memory
  dt->deallocTetrahedron(tet_idx);
}

static void finalize_infinite_tetrahedron(int tet_idx, Delaunay3* dt)
{
  Delaunay3Vertex* v;

  // make sure it's infinite
  assert(dt->tpool[tet_idx].V[0] == dt->pinf || dt->tpool[tet_idx].V[1] == dt->pinf);

  for (int j = 0; j < 4; j++)
  {
    v = dt->tpool[tet_idx].V[j];
    v->ref_count--;
    if (v->ref_count == 0)
    {
  	  dt->deallocVertex(v);	
    }
    else
    {
      assert(v->ref_count > 0);
    }
  }

  // make all its neighbor corners pointing to null

  int t1, i1; //indices of tetra and corner position
  if (dt->tpool[tet_idx].N[0]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[0]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[0]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[1]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[1]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[1]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[2]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[2]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[2]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }
  if (dt->tpool[tet_idx].N[3]!=D3_NULL_INDEX){
	  t1 = D3_TETRA(dt->tpool[tet_idx].N[3]);
	  i1 = D3_INDEX(dt->tpool[tet_idx].N[3]);
	  dt->tpool[t1].N[i1] = D3_NULL_INDEX;
  }

  // free its memory
  dt->deallocTetrahedron(tet_idx);
}

static bool isFinalizedAgainParent(int tet_idx, Delaunay3Tetrahedron* tet, const SScontainer3D* ss3d)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_again_parent++;
#endif
  // this tet was tested before for finalization. it's radius should be initialized.
  assert(tet->rad != -1.0f);
  // let's test this tet again. maybe it is now finalized. otherwise ... re-link it.
  return ss3d->is_sphere_finalized_parent(tet, tet_idx);
}

static bool isFinalizedAgainLeaf(int tet_idx, Delaunay3Tetrahedron* tet, const SScontainer3D* ss3d)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_again_leaf++;
#endif
  // this tet was tested before for finalization. it's radius should be initialized.
  assert(tet->rad != -1.0f);
  // let's test this tet again. maybe it is now finalized. otherwise ... re-link it.
  return ss3d->is_sphere_finalized_leaf(tet, tet_idx);
}

static bool isFinalizedFirst(int tet_idx, Delaunay3Tetrahedron* tet, SScontainer3D* ss3d)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_first++;
#endif
  // this is the first time that we check this tetrahedron for finalization.
  assert(tet->rad == -1.0f);
  // we need to compute its circumsphere
  tet->initialize(tet->V[0]->x, tet->V[1]->x, tet->V[2]->x, tet->V[3]->x);
  // often this circumsphere is completely inside the cell that was just finalized
  if (ss3d->was_sphere_just_finalized(tet)) return true;
  // otherwise me must perform a full check (and link the tet to a cell it may intersect)
  return ss3d->is_sphere_finalized(tet, tet_idx);
}

static void CheckFinalizeAndWriteParent(Delaunay3* dt, int tet_idx, SScontainer3D* ss3d, SVwriter* svwriter)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_cell++;
  stat_finalize_parent_cell++;
#endif

  // do not call for nothing
  assert(tet_idx != -1);

  int tet_next;
  // make sure the previous element points to me
  assert(dt->tpool[dt->tpool[tet_idx].prev].next == tet_idx);
  // mark the end of the doubly linked list
  dt->tpool[dt->tpool[tet_idx].prev].next = -1;

  while (tet_idx != -1)
  {
    // get the next tet in the linked list
    tet_next = dt->tpool[tet_idx].next;

    // do not process tetrahedra that no longer exist (but are still in the linked list)
    if (dt->tpool[tet_idx].RAD2_DEAD == RAD2_DEAD_TRUE)
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_cell_dead++;
      stat_delayed_deleted_tetrahedra--;
#endif
      // these tetrahedra were destroyed during an earlier point insertion but could not
      // be dealloced because they were still in this list. dealloc them now.
      dt->deallocTetrahedron(tet_idx);
      tet_idx = tet_next;
      continue;
    }

#ifdef COLLECT_STATISTICS
  stat_finalize_cell_checked++;
#endif

    // also ... currently we do not insert infinite tets in these lists
    assert(dt->isInf(tet_idx) == false);
    // check whether the tetrahedron's sphere is now completely inside
    // the finalized area or rather whether it still intersects other
    // unfinalized space elsewhere
    if (isFinalizedAgainParent(tet_idx, &(dt->tpool[tet_idx]), ss3d))
    {
      finalize_tetrahedron(tet_idx, dt, svwriter);
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_cell_relinked++;
#endif
      // make sure the tet was linked correctly
      assert(dt->tpool[tet_idx].cell->data == tet_idx);
      // complete link
      if (dt->tpool[tet_idx].next == -1) // was i the first element?
      {
        dt->tpool[tet_idx].prev = tet_idx;
        dt->tpool[tet_idx].next = tet_idx;
      }
      else
      {
        dt->tpool[tet_idx].prev = dt->tpool[dt->tpool[tet_idx].next].prev;
        dt->tpool[dt->tpool[tet_idx].prev].next = tet_idx;
        dt->tpool[dt->tpool[tet_idx].next].prev = tet_idx;
      }
    }
    // get next tet in linked list
    tet_idx = tet_next;
  }
}

static void FinalizeAndWriteLeaf(Delaunay3* dt, int tet_idx, SVwriter* svwriter)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_cell++;
  stat_finalize_leaf_cell++;
#endif

  // do not call for nothing
  assert(tet_idx != -1);
  
  int tet_next;
  // make sure the previous element points to me
  assert(dt->tpool[dt->tpool[tet_idx].prev].next == tet_idx);
  // mark the end of the doubly linked list
  dt->tpool[dt->tpool[tet_idx].prev].next = -1;


  while (tet_idx != -1)
  {
    // get the next tet in the linked list
    tet_next = dt->tpool[tet_idx].next;
    // do not process tetrahedra that no longer exist (but are still in the linked list)
    if (dt->tpool[tet_idx].RAD2_DEAD == RAD2_DEAD_TRUE)
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_cell_dead++;
      stat_delayed_deleted_tetrahedra--;
#endif
      // these tetrahedra were destroyed during an earlier point insertion but could not
      // be dealloced because of them being in this list. dealloc them now.
      dt->deallocTetrahedron(tet_idx);
    }
    else
    {
      // also ... currently we do not insert infinite tets in these lists
      assert(dt->isInf(tet_idx) == false);
      // finalize and write this tetrahedron
      finalize_tetrahedron(tet_idx, dt, svwriter);
    }
    tet_idx = tet_next;
  }
}

static void CheckFinalizeAndWriteLeaf(Delaunay3* dt, int tet_idx, SScontainer3D* ss3d, SVwriter* svwriter)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_cell++;
  stat_finalize_leaf_cell++;
#endif

  // do not call for nothing
  assert(tet_idx != -1);

  int tet_next;
  // make sure the previous element points to me
  assert(dt->tpool[dt->tpool[tet_idx].prev].next == tet_idx);
  // mark the end of the doubly linked list
  dt->tpool[dt->tpool[tet_idx].prev].next = -1;

  while (tet_idx != -1)
  {
    // get the next tet in the linked list
    tet_next = dt->tpool[tet_idx].next;

    // do not process tetrahedra that no longer exist (but are still in the linked list)
    if (dt->tpool[tet_idx].RAD2_DEAD == RAD2_DEAD_TRUE)
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_cell_dead++;
      stat_delayed_deleted_tetrahedra--;
#endif
      // these tetrahedra were destroyed during an earlier point insertion but could not
      // be dealloced because of them being in this list. dealloc them now.
      dt->deallocTetrahedron(tet_idx);
      tet_idx = tet_next;
      continue;
    }

#ifdef COLLECT_STATISTICS
    stat_finalize_cell_checked++;
#endif

    // also ... currently we do not insert infinite tets in these lists
    assert(dt->isInf(tet_idx) == false);
    // check whether the tetrahedron's sphere is now completely inside
    // the finalized area or rather whether it still intersects other
    // unfinalized space elsewhere
    if (isFinalizedAgainLeaf(tet_idx, &(dt->tpool[tet_idx]), ss3d))
    {
      finalize_tetrahedron(tet_idx, dt, svwriter);
    }
    else
    {
#ifdef COLLECT_STATISTICS
      stat_finalize_cell_relinked++;
#endif
      // make sure the tet was linked correctly
      assert(dt->tpool[tet_idx].cell->data == tet_idx);
      // complete link
      if (dt->tpool[tet_idx].next == -1) // was i the first element?
      {
        dt->tpool[tet_idx].prev = tet_idx;
        dt->tpool[tet_idx].next = tet_idx;
      }
      else
      {
        dt->tpool[tet_idx].prev = dt->tpool[dt->tpool[tet_idx].next].prev;
        dt->tpool[dt->tpool[tet_idx].prev].next = tet_idx;
        dt->tpool[dt->tpool[tet_idx].next].prev = tet_idx;
      }
    }
    // get next tet in linked list
    tet_idx = tet_next;
  }
}

static void FinalizeAndWriteNewTets(Delaunay3* dt, SVwriter* svwriter)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_new++;
#endif
  int i;
  for (int tetrahedron_newchunk = 0; tetrahedron_newchunk < dt->tetrahedron_newchunk_number; tetrahedron_newchunk++)
  {
    i = dt->tetrahedron_newchunks[tetrahedron_newchunk] * DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
	  for (int newchunk_elements = 0; newchunk_elements < DELAUNAY3TETRAHEDRON_CHUNK_SIZE; newchunk_elements++, i++)
    {
      // do not check dead tetrahedra
 		  if (dt->tpool[i].RAD2_DEAD == RAD2_DEAD_TRUE)
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_dead++;
  #endif
        continue;
      }

      // do not check infinite tetrahedra
      if (dt->isInf(i))
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_infinite++;
  #endif
        continue;
      }

      // finalize and write this tetrahedron
      finalize_tetrahedron(i, dt, svwriter);
    }
  }
}

static void CheckFinalizeAndWriteNewTets(Delaunay3* dt, SScontainer3D* ss3d, SVwriter* svwriter)
{
#ifdef COLLECT_STATISTICS
  stat_finalize_new++;
#endif

  int i;
  for (int tetrahedron_newchunk = 0; tetrahedron_newchunk < dt->tetrahedron_newchunk_number; tetrahedron_newchunk++)
  {
    i = dt->tetrahedron_newchunks[tetrahedron_newchunk] * DELAUNAY3TETRAHEDRON_CHUNK_SIZE;
	  for (int newchunk_elements = 0; newchunk_elements < DELAUNAY3TETRAHEDRON_CHUNK_SIZE; newchunk_elements++, i++)
    {
      // do not check dead tetrahedra
 		  if (dt->tpool[i].RAD2_DEAD == RAD2_DEAD_TRUE)
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_dead++;
  #endif
        continue;
      }

#ifdef DONT_CHECK_INFINITE_TETS_FOR_BB_ALIGNGMENT
      // do not check infinite tetrahedra
      if (dt->isInf(i))
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_infinite++;
  #endif
        continue;
      }
#else // ---------------------------- DONT_CHECK_INFINITE_TETS_FOR_BB_ALIGNGMENT ------------------------------ 
      // if this is an infinite tetrahedron with V[0] being infinite
      if (dt->tpool[i].V[0] == dt->pinf)
      {
        // if the tet neighbor adjacent at the finite triangle was finalized
        if (dt->tpool[i].N[0] == D3_NULL_INDEX)
        {
          // if the three finite vertices have the same x coordinate
          if (dt->tpool[i].V[1]->x[0] == dt->tpool[i].V[2]->x[0] && dt->tpool[i].V[2]->x[0] == dt->tpool[i].V[3]->x[0])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[0] == (double)ss3d->bb_min_f[0] || dt->tpool[i].V[2]->x[0] == (double)ss3d->bb_max_f[0])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }
          // or the same y coordinate
          else if (dt->tpool[i].V[1]->x[1] == dt->tpool[i].V[2]->x[1] && dt->tpool[i].V[2]->x[1] == dt->tpool[i].V[3]->x[1])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[1] == (double)ss3d->bb_min_f[1] || dt->tpool[i].V[2]->x[1] == (double)ss3d->bb_max_f[1])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }
          // or the same z coordinate
          else if (dt->tpool[i].V[1]->x[2] == dt->tpool[i].V[2]->x[2] && dt->tpool[i].V[2]->x[2] == dt->tpool[i].V[3]->x[2])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[2] == (double)ss3d->bb_min_f[2] || dt->tpool[i].V[2]->x[2] == (double)ss3d->bb_max_f[2])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }        
        }
#ifdef COLLECT_STATISTICS
        stat_finalize_new_infinite++;
#endif
        continue;
      }
      // if this is an infinite tetrahedron with V[1] being infinite
      else if (dt->tpool[i].V[1] == dt->pinf)
      {
        // if the tet neighbor adjacent at the finite triangle was finalized
        if (dt->tpool[i].N[1] == D3_NULL_INDEX)
        {
          // if the three finite vertices have the same x coordinate
          if (dt->tpool[i].V[0]->x[0] == dt->tpool[i].V[2]->x[0] && dt->tpool[i].V[2]->x[0] == dt->tpool[i].V[3]->x[0])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[0] == (double)ss3d->bb_min_f[0] || dt->tpool[i].V[2]->x[0] == (double)ss3d->bb_max_f[0])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }
          // or the same y coordinate
          else if (dt->tpool[i].V[0]->x[1] == dt->tpool[i].V[2]->x[1] && dt->tpool[i].V[2]->x[1] == dt->tpool[i].V[3]->x[1])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[1] == (double)ss3d->bb_min_f[1] || dt->tpool[i].V[2]->x[1] == (double)ss3d->bb_max_f[1])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }
          // or the same z coordinate
          else if (dt->tpool[i].V[0]->x[2] == dt->tpool[i].V[2]->x[2] && dt->tpool[i].V[2]->x[2] == dt->tpool[i].V[3]->x[2])
          {
            // and if this coordinate is on the bounding box
            if (dt->tpool[i].V[2]->x[2] == (double)ss3d->bb_min_f[2] || dt->tpool[i].V[2]->x[2] == (double)ss3d->bb_max_f[2])
            {
              finalize_infinite_tetrahedron(i, dt);
#ifdef COLLECT_STATISTICS
              stat_finalize_new_infinite_finalized++;
#endif
              continue;
            }
          }
        }
#ifdef COLLECT_STATISTICS
        stat_finalize_new_infinite++;
#endif
        continue;
      }
#endif  // ---------------------------- DONT_CHECK_INFINITE_TETS_FOR_BB_ALIGNGMENT ------------------------------ 

      // do not check tetrahedra that are linked by a cell
 		  if (dt->tpool[i].cell)
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_with_cell++;
  #endif
        continue;
      }

  #ifdef COLLECT_STATISTICS
      stat_finalize_new_checked++;
  #endif

      // it seems we have an un-dead, finite, and brand-new tetrahedron
      // check whether the tetrahedron's sphere is completely inside
      // the finalized area
      if (isFinalizedFirst(i, &(dt->tpool[i]), ss3d))
      {
        finalize_tetrahedron(i, dt, svwriter);
      }
      else
      {
  #ifdef COLLECT_STATISTICS
        stat_finalize_new_linked++;
  #endif
        // make sure the tet was linked correctly
        assert(dt->tpool[i].cell->data == i);
        // complete link
        if (dt->tpool[i].next == -1) // was i the first element?
        {
          dt->tpool[i].prev = i;
          dt->tpool[i].next = i;
        }
        else
        {
          dt->tpool[i].prev = dt->tpool[dt->tpool[i].next].prev;
          dt->tpool[dt->tpool[i].prev].next = i;
          dt->tpool[dt->tpool[i].next].prev = i;
        }
      }
    }
  }
}

// find the subset of the finalized triangles whose vertices
// are all finalized and free up their memory

bool SPdelaunay3D::open(SVwriter* svwriter)
{
  this->svwriter = svwriter;
  this->ss3d = new SScontainer3D();
  this->dt = new Delaunay3();
  p_count = 0;
  return true;
}

void SPdelaunay3D::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  ss3d->open(bb_min_f, bb_max_f);

  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];

  VecCopy3fv(this->bb_min_f, bb_min_f);
  VecCopy3fv(this->bb_max_f, bb_max_f);

  svwriter->set_v_pnum_f(3);
  svwriter->set_v_pbox_f(bb_min_f, bb_max_f);
}

void SPdelaunay3D::set_boundingbox(const double* bb_min_d, const double* bb_max_d)
{
  ss3d->open(bb_min_d, bb_max_d);

  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];

  VecCopy3fv(this->bb_min_f, bb_min_d);
  VecCopy3fv(this->bb_max_f, bb_max_d);

  svwriter->set_v_pnum_f(3);
  svwriter->set_v_pbox_f(this->bb_min_f, this->bb_max_f);
}

static Delaunay3Vertex* p0123[4];

void SPdelaunay3D::write_point(const float* p_pos_f)
{
  Delaunay3Vertex* p = dt->allocVertex(p_pos_f);
  if (p_count >= 4)
  {
  	dt->insert(p);
  }
  else
  {
    p0123[p_count] = p;

    if (p_count == 3)
    {
  	  dt->initialize(p0123);
    }
  }
  p_count++;
}

void SPdelaunay3D::write_finalize_cell(int cell_idx)
{
  // first we finalize the cell (and all its children)
  int type = ss3d->finalize_cell(cell_idx);
  // if p_count is still zero then this was only a finalize empty call
  if (p_count)
  {
    // maybe the root was finalized 
    if (cell_idx == 0)
    {
      // everything is finalized so we write out everything
      assert(type != 1);
      // write out tetrahedra stored with finalized leafs
      while ((type = ss3d->prepareLeaf()) != -1)
      {
        FinalizeAndWriteLeaf(dt, type, svwriter);
      }
      // then write out all new tetrahedra
      FinalizeAndWriteNewTets(dt, svwriter);
    }
    else
    {
      // is the data stored with the parent of a finalized cell 
      if (type == 1)
      {
        // check tetrahedra stored with parent of finalized cell
        if ((type = ss3d->prepareParent()) != -1)
        {
          CheckFinalizeAndWriteParent(dt, type, ss3d, svwriter);
        }
      }
      else if (type == 2)
      {
        // check tetrahedra stored with finalized cell
        while ((type = ss3d->prepareLeaf()) != -1)
        {
          CheckFinalizeAndWriteLeaf(dt, type, ss3d, svwriter);
        }
      }
      // then we loop over all new tetrahedra
      CheckFinalizeAndWriteNewTets(dt, ss3d, svwriter);
      dt->resetDelaunay3TetrahedronNewChunks();
    }
  }
#ifndef NDEBUG		
//	dt->audit();
#endif
}

void SPdelaunay3D::close()
{
  if (bb_min_f == 0 || bb_max_f == 0)
  {
    // this finalize whatever space is remaining
    ss3d->close();
    // write whatever is remaining
  	FinalizeAndWriteNewTets(dt, svwriter);
#ifndef NDEBUG
//    dt->audit();
#endif
  } 

#ifdef COLLECT_STATISTICS
  fprintf(stderr, "new (%d) : dead %d infinite %d (fin %d) with_cell %d checked %d linked %d\n", stat_finalize_new, stat_finalize_new_dead, stat_finalize_new_infinite,stat_finalize_new_infinite_finalized, stat_finalize_new_with_cell, stat_finalize_new_checked, stat_finalize_new_linked);
  fprintf(stderr, "cell (%d=%d+%d) : dead %d checked %d relinked %d\n", stat_finalize_cell, stat_finalize_leaf_cell, stat_finalize_parent_cell, stat_finalize_cell_dead, stat_finalize_cell_checked, stat_finalize_cell_relinked);
  fprintf(stderr, "first %d again %d = %d + %d \n", stat_finalize_first, stat_finalize_again_parent+stat_finalize_again_leaf, stat_finalize_again_parent, stat_finalize_again_leaf);
  fprintf(stderr, "max number of allocated vertices/tetrahedra: %d/%d (%d/%d)\n", dt->delaunay_vertex_buffer_maxsize, dt->tetrahedron_buffer_maxsize, dt->delaunay_vertex_buffer_size, dt->tetrahedron_buffer_used);
#endif

  // delete spatial finalizer

  delete ss3d; ss3d = 0;
  delete dt; dt = 0;
}

SPdelaunay3D::SPdelaunay3D()
{
  // init of SPwriter interface
  ncomments = 0;
  comments = 0;

  npoints = -1;
  p_count = -1;

  datatype = SP_FLOAT;

  bb_min_d = 0;
  bb_max_d = 0;

  bb_min_f = 0;
  bb_max_f = 0;

  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPdelaunay3D interface
  svwriter = 0;
  ss3d = 0;
  dt = 0;
}

void SPdelaunay3D::getActiveTetrahedraInit()
{
  active_tetrahedra_next = 0;
}

int SPdelaunay3D::getActiveTetrahedraNext(float* v0, float* v1, float* v2, float* v3)
{
  while (active_tetrahedra_next <= dt->tetrahedron_buffer_maxsize)
  {
    if (dt->tpool[active_tetrahedra_next].RAD2_DEAD == RAD2_DEAD_TRUE)
    {
      active_tetrahedra_next++;
      continue;
    }

    if (dt->tpool[active_tetrahedra_next].V[0] == dt->pinf || dt->tpool[active_tetrahedra_next].V[1] == dt->pinf)
    {
      active_tetrahedra_next++;
      continue;
    }

    Delaunay3Vertex* v;

    v = dt->tpool[active_tetrahedra_next].V[0];
    v0[0] = (float)v->x[0];
    v0[1] = (float)v->x[1];
    v0[2] = (float)v->x[2];

		v = dt->tpool[active_tetrahedra_next].V[1];
    v1[0] = (float)v->x[0];
    v1[1] = (float)v->x[1];
    v1[2] = (float)v->x[2];

		v = dt->tpool[active_tetrahedra_next].V[2];
    v2[0] = (float)v->x[0];
    v2[1] = (float)v->x[1];
    v2[2] = (float)v->x[2];
    
    v = dt->tpool[active_tetrahedra_next].V[3];
    v3[0] = (float)v->x[0];
    v3[1] = (float)v->x[1];
    v3[2] = (float)v->x[2];

    if (dt->tpool[active_tetrahedra_next].rad2 == -2.0f) // hack to visualize inverted tets (marked with rad == -2)
    {
      active_tetrahedra_next++;
      return -2;
    }
    return active_tetrahedra_next++;
  }
  return -1;
}

bool SPdelaunay3D::getActiveTetrahedraCoords(int idx, float* v0, float* v1, float* v2, float* v3)
{
  if (dt->tpool[idx].RAD2_DEAD == RAD2_DEAD_TRUE)
  {
    return false;
  }

  if (dt->tpool[idx].V[0] == dt->pinf || dt->tpool[idx].V[1] == dt->pinf)
  {
    return false;
  }

  v0[0] = (float)dt->tpool[idx].V[0]->x[0];
  v0[1] = (float)dt->tpool[idx].V[0]->x[1];
  v0[2] = (float)dt->tpool[idx].V[0]->x[2];
  v1[0] = (float)dt->tpool[idx].V[1]->x[0];
  v1[1] = (float)dt->tpool[idx].V[1]->x[1];
  v1[2] = (float)dt->tpool[idx].V[1]->x[2];
  v2[0] = (float)dt->tpool[idx].V[2]->x[0];
  v2[1] = (float)dt->tpool[idx].V[2]->x[1];
  v2[2] = (float)dt->tpool[idx].V[2]->x[2];
  v3[0] = (float)dt->tpool[idx].V[3]->x[0];
  v3[1] = (float)dt->tpool[idx].V[3]->x[1];
  v3[2] = (float)dt->tpool[idx].V[3]->x[2];

  return true;
}

bool SPdelaunay3D::getActiveTetrahedraSphere(int idx, float* cen_rad)
{
  if (dt->tpool[idx].RAD2_DEAD == RAD2_DEAD_TRUE)
  {
    return false;
  }

  if (dt->tpool[idx].V[0] == dt->pinf || dt->tpool[idx].V[1] == dt->pinf)
  {
    return false;
  }

  if (dt->tpool[idx].rad == -1.0f)
  {
    fprintf(stderr, "WARNING: circumsphere of tet %d not computed yet\n", idx);
    return false;
  }

  cen_rad[0] = dt->tpool[idx].cen[0];
  cen_rad[1] = dt->tpool[idx].cen[1];
  cen_rad[2] = dt->tpool[idx].cen[2];
  cen_rad[3] = dt->tpool[idx].rad;

  return true;
}

bool SPdelaunay3D::getActiveTetrahedraCellIdx(int idx, int* cell_idx)
{
  if (dt->tpool[idx].RAD2_DEAD == RAD2_DEAD_TRUE)
  {
    return false;
  }

  if (dt->tpool[idx].V[0] == dt->pinf || dt->tpool[idx].V[1] == dt->pinf)
  {
    return false;
  }

  if (dt->tpool[idx].cell == 0)
  {
    return false;
  }

  *cell_idx = dt->tpool[idx].cell->idx;

  return true;
}

void SPdelaunay3D::getInfiniteTetrahedraBaseTriangleInit()
{
  infinite_tetrahedra_next = 0;
}

int SPdelaunay3D::getInfiniteTetrahedraBaseTriangleNext(float* v0, float* v1, float* v2)
{
  while (infinite_tetrahedra_next < dt->tetrahedron_buffer_maxsize)
  {
    if (dt->tpool[infinite_tetrahedra_next].RAD2_DEAD == RAD2_DEAD_TRUE)
    {
      infinite_tetrahedra_next++;
      continue;
    }

    Delaunay3Vertex* v;

    v = dt->tpool[infinite_tetrahedra_next].V[0];
    if (v == dt->pinf)
    {
			v = dt->tpool[infinite_tetrahedra_next].V[1]; 
      v0[0] = (float)v->x[0];
      v0[1] = (float)v->x[1];
      v0[2] = (float)v->x[2];

      v = dt->tpool[infinite_tetrahedra_next].V[2]; 
      v1[0] = (float)v->x[0];
      v1[1] = (float)v->x[1];
      v1[2] = (float)v->x[2];

      v = dt->tpool[infinite_tetrahedra_next].V[3]; 
      v2[0] = (float)v->x[0];
      v2[1] = (float)v->x[1];
      v2[2] = (float)v->x[2];

      infinite_tetrahedra_next++;
      return true;
    }
		v = dt->tpool[infinite_tetrahedra_next].V[1];
    if (v == dt->pinf)
    {
     	v = dt->tpool[infinite_tetrahedra_next].V[0]; 
      v0[0] = (float)v->x[0];
      v0[1] = (float)v->x[1];
      v0[2] = (float)v->x[2];

    	v = dt->tpool[infinite_tetrahedra_next].V[2]; 
      v1[0] = (float)v->x[0];
      v1[1] = (float)v->x[1];
      v1[2] = (float)v->x[2];

    	v = dt->tpool[infinite_tetrahedra_next].V[3]; 
      v2[0] = (float)v->x[0];
      v2[1] = (float)v->x[1];
      v2[2] = (float)v->x[2];

      infinite_tetrahedra_next++;
      return 1;
    }
    infinite_tetrahedra_next++;
  }
  return 0;
}


