/*
===============================================================================

  FILE:  SPreadScattered.h
  
  CONTENTS:
  
    Reads a Streaming Mesh with a "little-cache" aware greedy reordering of
    triangles that is subject to a constraint of maximal delay of triangles.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    08 September 2005 -- created after cancelling the streaming talk for SIAM
  
===============================================================================
*/
#ifndef SPREAD_SCATTERED_H
#define SPREAD_SCATTERED_H

#include "spreader.h"

class SPreadScattered : public SPreader
{
public:
  // spreader interface function implementations

  SPevent read_event();

  void close();

  // spreadscattered functions

  bool open(SPreader* spreader, int max_scatter=100000, int fct_scatter=50);

  SPreadScattered();
  ~SPreadScattered();

private:
  SPreader* spreader;
  int max_scatter;
  int fct_scatter;
  int buffered_points;

  void fill_buffer();
};

#endif
