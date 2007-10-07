/*
===============================================================================

  FILE:  SPconverter.cpp
  
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
#include "spconverter.h"

#include "vec3fv.h"

bool SPconverter::open(SMreader* smreader)
{
  this->smreader = smreader;
  this->svreader = 0;

  if (smreader->nverts != -1)
  {
    npoints = smreader->nverts;
  }
  if (smreader->bb_min_f)
  {
    if (bb_min_f == 0) bb_min_f = new float[3];
    VecCopy3fv(bb_min_f, smreader->bb_min_f);
  }
  if (smreader->bb_max_f)
  {
    if (bb_max_f == 0) bb_max_f = new float[3];
    VecCopy3fv(bb_max_f, smreader->bb_max_f);
  }

  p_count = 0;

  return true;
}

bool SPconverter::open(SVreader* svreader)
{
  this->svreader = svreader;
  this->smreader = 0;

  if (svreader->nverts != -1)
  {
    npoints = svreader->nverts;
  }
  if (svreader->v_pmin_f)
  {
    if (bb_min_f == 0) bb_min_f = new float[3];
    VecCopy3fv(bb_min_f, svreader->v_pmin_f);
  }
  if (svreader->v_pmax_f)
  {
    if (bb_max_f == 0) bb_max_f = new float[3];
    VecCopy3fv(bb_max_f, svreader->v_pmax_f);
  }

  p_count = 0;

  return true;
}

void SPconverter::close()
{
  if (smreader) smreader->close();
  else if (svreader) svreader->close();
  p_count = -1;
}

SPevent SPconverter::read_event()
{
  if (smreader)
  {
    SMevent event;
    while (event = smreader->read_element())
    {
      if (event == SM_VERTEX)
      {
        VecCopy3fv(p_pos_f, smreader->v_pos_f);
        p_count++;
        return SP_POINT;
      }
      else if (event == SM_ERROR)
      {
        return SP_ERROR;
      }
    }
  }
  else if (svreader)
  {
    SVevent event;
    while (event = svreader->read_element())
    {
      if (event == SV_VERTEX)
      {
        VecCopy3fv(p_pos_f, svreader->v_prop_f);
        p_count++;
        return SP_POINT;
      }
      else if (event == SV_ERROR)
      {
        return SP_ERROR;
      }
    }
  }
  return SP_EOF;
}

SPconverter::SPconverter()
{
  // init of SPreader interface
  ncomments = 0;
  comments = 0;

  datatype = SP_VOID;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPconverter interface
  smreader = 0;
  svreader = 0;
}

SPconverter::~SPconverter()
{
  if (bb_min_d) delete [] bb_min_d;
  if (bb_max_d) delete [] bb_max_d;
  if (bb_min_f) delete [] bb_min_f;
  if (bb_max_f) delete [] bb_max_f;
  if (bb_min_i) delete [] bb_min_i;
  if (bb_max_i) delete [] bb_max_i;
}
