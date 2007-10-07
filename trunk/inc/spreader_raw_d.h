/*
===============================================================================

  FILE:  SPreader_spb.h
  
  CONTENTS:
  
    Reads points from raw binary double-precision format.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 January 2006 -- created to have an efficient reader for raw double data
  
===============================================================================
*/
#ifndef SPREADER_RAW_D_H
#define SPREADER_RAW_D_H

#include "spreader.h"

#include <stdio.h>

class SPreader_raw_d : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_spb functions

  bool open(FILE* file, bool precompute_bounding_box=false);
  void compute_bounding_box();
  bool load_header(FILE* file);

  SPreader_raw_d();
  ~SPreader_raw_d();

private:
  FILE* file;
};

#endif
