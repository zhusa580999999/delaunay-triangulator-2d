/*
===============================================================================

  FILE:  SPreader_tiles.h
  
  CONTENTS:
  
    Reads points from a spreader and replicates them into an x*y tiling to
    create an even larger terrain.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    20 January 2006 -- created after the "initial" SIGGRAPH meeting in Soda
  
===============================================================================
*/
#ifndef SPREADER_TILES_H
#define SPREADER_TILES_H

#include "spreader.h"

#include <stdio.h>

class SPreader_tiles : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_spb functions

  bool open(char* file_name_raw, char* file_name_hdr, int x, int y);

  SPreader_tiles();
  ~SPreader_tiles();

private:
  FILE* file;
  char* file_name;
  double bb_range_d[2];
  int x_max, x_cur;
  int y_max, y_cur;
};

#endif
