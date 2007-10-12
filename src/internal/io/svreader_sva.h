/*
===============================================================================

  FILE:  SVreader_sva.h
  
  CONTENTS:
  
    Reads a tetrahedral volume mesh from a Streaming Mesh ASCII format (SMA).
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    20 February 2005 -- created after streaming points turned out too tricky
  
===============================================================================
*/
#ifndef SVREADER_SMA_H
#define SVREADER_SMA_H

#include "svreader.h"

#include <stdio.h>

class SVreader_sva : public SVreader
{
public:

  // smreader interface function implementations

  void close();

  SVevent read_element();
  SVevent read_event();

  // smreader_sma functions

  bool open(FILE* file);

  SVreader_sva();
  ~SVreader_sva();

private:
  FILE* file;
  int skipped_lines;
  char* line;
  int have_finalized, next_finalized;
  int* finalized_vertices;
};

#endif
