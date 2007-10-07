/*
===============================================================================

  FILE:  SPreadScattered.cpp
  
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
#include "spreadscattered.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "vec3fv.h"
#include "vec3iv.h"

#define PRINT_CONTROL_OUTPUT
//#undef PRINT_CONTROL_OUTPUT

// data structures used for efficient memory management

static int scatterbuffer_num = 0;
static int scatterbuffer_last = 0;
static int scatterbuffer_next = 0;
static int scatterbuffer_alloc = 0;
static float* scatterbuffer = 0;

static int finalizebuffer_num = 0;
static int finalizebuffer_last = 0;
static int finalizebuffer_next = 0;
static int finalizebuffer_alloc = 0;
static int* finalizebuffer = 0;

static int next_finalize = -1;
static int next_scatter = -1;
static int next_point = -1;

static int scatter_0_step = 0;
static int scatter_1_step = 0;
static int scatter_2_step = 0;

static int scatter_0_next = -1;
static int scatter_1_next = -1;
static int scatter_2_next = -1;

static int scatter_0_more = 0;
static int scatter_1_more = 0;
static int scatter_2_more = 0;

bool SPreadScattered::open(SPreader* spreader, int max_scatter, int fct_scatter)
{
  if (spreader == 0 )
  {
    return false;
  }
  this->spreader = spreader;
  this->max_scatter = max_scatter;
  this->fct_scatter = fct_scatter;

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"maximal point scatter is %d\n", max_scatter);
#endif
  
  npoints = spreader->npoints;

  p_count = 0;

  bb_min_f = spreader->bb_min_f;
  bb_max_f = spreader->bb_max_f;

  scatterbuffer_alloc = max_scatter;
  scatterbuffer = (float*)malloc(sizeof(float)*4*scatterbuffer_alloc);

  finalizebuffer_alloc = max_scatter/10;
  finalizebuffer = (int*)malloc(sizeof(int)*2*finalizebuffer_alloc);

  next_point = 0;
  next_scatter = 0;
  scatter_0_step = max_scatter/10;
  scatter_1_step = max_scatter/100;
  scatter_2_step = max_scatter/1000;
  scatter_0_more = 10;
  scatter_1_more = 10;
  scatter_2_more = 100;

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"filling buffer ... ");
#endif

  fill_buffer();

#ifdef PRINT_CONTROL_OUTPUT
  fprintf(stderr,"%d %d.\n", scatterbuffer_num, finalizebuffer_num);
#endif

  return true;
}

void SPreadScattered::close()
{
  if (npoints != -1) if (npoints != p_count)  fprintf(stderr,"WARNING: got npoints %d but p_count %d\n",npoints,p_count);

  p_count = -1;

  free(scatterbuffer);
  free(finalizebuffer);
}

void SPreadScattered::fill_buffer()
{
  SPevent event;

  while (scatterbuffer_num < scatterbuffer_alloc && next_finalize != next_point)
  {
    event = spreader->read_event();
    if (event == SP_POINT)
    {
       // copy point into buffer
      scatterbuffer[4*scatterbuffer_last+0] = 1; // the point is not output yet
      VecCopy3fv(&(scatterbuffer[4*scatterbuffer_last+1]), spreader->p_pos_f);
      // increment pointer
      scatterbuffer_last++;
      if (scatterbuffer_last == scatterbuffer_alloc) scatterbuffer_last = 0;
      // increment counter
      scatterbuffer_num++;
    }
    else if (event == SP_FINALIZED_CELL)
    {
      if (next_finalize == -1) next_finalize = spreader->p_count;
      // copy finalization event into buffer
      finalizebuffer[2*finalizebuffer_last+0] = spreader->p_count; // finalize as soon as input point count
      finalizebuffer[2*finalizebuffer_last+1] = spreader->final_idx;
      // increment pointer
      finalizebuffer_last++;
      if (finalizebuffer_last == finalizebuffer_alloc) finalizebuffer_last = 0;
      // increment counter
      finalizebuffer_num++;
      // make sure we do not exceed our memory resources
      assert(finalizebuffer_num <= finalizebuffer_alloc);
    }
    else
    {
      spreader->close();
      spreader = 0;
      return;
    }
  }
}

SPevent SPreadScattered::read_event()
{
  while (true)
  {
    // first check if we need to have a finalize event
    if (next_finalize == next_point)
    {
      // just making sure
      assert(finalizebuffer[2*finalizebuffer_next+0] == next_point);
      // copy finalization event
      final_idx = finalizebuffer[2*finalizebuffer_next+1];
      // increment pointer
      finalizebuffer_next++;
      if (finalizebuffer_next == finalizebuffer_alloc) finalizebuffer_next = 0;
      // decrement counter
      finalizebuffer_num--;
      // when to do the next finalize event
      if (finalizebuffer_num)
      {
        next_finalize = finalizebuffer[2*finalizebuffer_next+0];
      }
      else
      {
        next_finalize = -1;
      }
      // refill buffer
      if (spreader)
      {
        fill_buffer();
      }
      return SP_FINALIZED_CELL;
    }
    else if (scatterbuffer_num > 0)
    {
      // then check if we need to do some scattering
      if (next_scatter == next_point)
      {
        while (scatter_0_more)
        {
          scatter_0_more--;
          scatter_0_next += scatter_0_step;
          if (scatter_0_next >= scatterbuffer_alloc) scatter_0_next -= scatterbuffer_alloc;
          // does the point exist
          if (scatterbuffer[scatter_0_next*4])
          {
            scatterbuffer[scatter_0_next*4] = 0;
            VecCopy3fv(p_pos_f, &(scatterbuffer[scatter_0_next*4+1]));
            // increment point counter
            p_count++;
            return SP_POINT;
          }
        }

        while (scatter_1_more)
        {
          scatter_1_more--;
          scatter_1_next += scatter_1_step;
          if (scatter_1_next >= scatterbuffer_alloc) scatter_1_next -= scatterbuffer_alloc;
          if (scatterbuffer[scatter_1_next*4])
          {
            scatterbuffer[scatter_1_next*4] = 0;
            VecCopy3fv(p_pos_f, &(scatterbuffer[scatter_1_next*4+1]));
            // increment point counter
            p_count++;
            return SP_POINT;
          }
        }

        while (scatter_2_more)
        {
          scatter_2_more--;
          scatter_2_next += scatter_2_step;
          if (scatter_2_next >= scatterbuffer_alloc) scatter_2_next -= scatterbuffer_alloc;
          if (scatterbuffer[scatter_2_next*4])
          {
            scatterbuffer[scatter_2_next*4] = 0;
            VecCopy3fv(p_pos_f, &(scatterbuffer[scatter_2_next*4+1]));
            // increment point counter
            p_count++;
            return SP_POINT;
          }
        }

        // nothing more to scatter ... setup next scatter

        next_scatter += scatter_0_step;
        scatter_0_more = 1;
        scatter_1_more = 10;
        scatter_2_more = 100;
      }

      // get the next point

      if (scatterbuffer[scatterbuffer_next*4])
      {
        scatterbuffer[scatterbuffer_next*4] = 0;
        VecCopy3fv(p_pos_f, &(scatterbuffer[scatterbuffer_next*4+1]));

        // advance pointer
        scatterbuffer_next++;
        if (scatterbuffer_next == scatterbuffer_alloc) scatterbuffer_next = 0;

        // increment input point counter
        next_point++;
        // decrement number of points in buffer
        scatterbuffer_num--;

        // increment point counter
        p_count++;
        return SP_POINT;
      }
      else
      {
        // advance pointer
        scatterbuffer_next++;
        if (scatterbuffer_next == scatterbuffer_alloc) scatterbuffer_next = 0;

        // increment input point counter
        next_point++;
        // decrement number of points in buffer
        scatterbuffer_num--;

        // refill buffer
        if (spreader)
        {
          fill_buffer();
        }
      }
    }
    else
    {
      return SP_EOF;
    }
  }
}

SPreadScattered::SPreadScattered()
{
  // init of SPreader interface
  ncomments = 0;
  comments = 0;

  npoints = -1;

  p_count = -1;

  bb_min_f = 0;
  bb_max_f = 0;

  // init of SPreadScattered
  spreader = 0;
  max_scatter = -1;
  buffered_points = -1;
}

SPreadScattered::~SPreadScattered()
{
}
