/*
===============================================================================

  FILE:  SVreader_svc.h
  
  CONTENTS:
  
    Reads a Streaming Volume Mesh from a compressed format.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    22 April 2005 -- completed version v0.0 which only uses a vertex cache
    12 April 2005 -- created after completing the sketch for SIGGRAPH 2005

===============================================================================
*/
#ifndef SVREADER_SVC_H
#define SVREADER_SVC_H

#include "svreader.h"

#include <stdio.h>

class SVreader_svc : public SVreader
{
public:

  // additional mesh variables
  int nbits;

  // svreader interface function implementations

  void close();

  SVevent read_element();
  SVevent read_event();

  // svreader_svc functions

  bool open(FILE* file);

  SVreader_svc();
  ~SVreader_svc();

private:
  int have_new, next_new, *new_vertices;
  int have_tetrahedron;
  int have_finalized, next_finalized, *finalized_vertices;

  void read_header();
  int decompress_tetrahedron();
};

#endif
