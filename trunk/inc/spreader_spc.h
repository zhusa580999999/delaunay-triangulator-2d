/*
===============================================================================

  FILE:  SPreader_spc.h
  
  CONTENTS:
  
    Reads a point cloud from a compressed Streaming Point format (SPC).

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    16 February 2005 -- created after Stefan's point compressor was ditched
  
===============================================================================
*/
#ifndef SPREADER_SPC_H
#define SPREADER_SPC_H

#include "spreader.h"

#include <stdio.h>

class SPreader_spc : public SPreader
{
public:
  // additional variables

  int nbits;

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_spc functions

  bool open(FILE* file);

  SPreader_spc();
  ~SPreader_spc();

private:
  void read_header();
  int decompress_kdtree();
};

#endif
