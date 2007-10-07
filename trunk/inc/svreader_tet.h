/*
===============================================================================

  FILE:  SVreader_tet.h
  
  CONTENTS:
  
    Reads a mesh from my TET format and presents it as a (poor) streaming mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    22 February 2005 -- created after washing and drying two loads of laundry 
  
===============================================================================
*/
#ifndef SVREADER_TET_H
#define SVREADER_TET_H

#include "svreader.h"

#include <stdio.h>

class SVreader_tet : public SVreader
{
public:

  // smreader interface function implementations

  void close();

  SVevent read_element();
  SVevent read_event();

  // smreader_sma functions

  bool open(FILE* fp, bool compute_boundingbox=true);

  SVreader_tet();
  ~SVreader_tet();

private:
  FILE* file;
  int skipped_lines;
  int minus;
  char* line;
};

#endif
