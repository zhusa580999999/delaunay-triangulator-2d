/*
===============================================================================

  FILE:  SPreader_spb.h
  
  CONTENTS:
  
    Reads points from raw binary single-precision format.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    19 January 2006 -- created to have an efficient reader for float data
  
===============================================================================
*/
#ifndef SPREADER_RAW_H
#define SPREADER_RAW_H

#include "spreader.h"

#include <stdio.h>

class SPreader_raw : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_spb functions

  bool open(FILE* file, bool precompute_bounding_box=false);
  void compute_bounding_box();
  bool load_header(FILE* file);

  SPreader_raw();
  ~SPreader_raw();

private:
  FILE* file;
};

#endif
