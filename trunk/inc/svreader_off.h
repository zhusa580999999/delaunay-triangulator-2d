/*
===============================================================================

  FILE:  SVreader_off.h
  
  CONTENTS:
  
    Reads a mesh from an OFF format and presents it as a (poor) streaming mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    08 April 2005 -- created after doing my three US tax returns for 2004
  
===============================================================================
*/
#ifndef SVREADER_OFF_H
#define SVREADER_OFF_H

#include "svreader.h"

#include <stdio.h>

class SVreader_off : public SVreader
{
public:

  // smreader interface function implementations

  void close();

  SVevent read_element();
  SVevent read_event();

  // smreader_sma functions

  bool open(FILE* fp, bool compute_boundingbox=false);

  SVreader_off();
  ~SVreader_off();

private:
  FILE* file;
  int skipped_lines;
  char* line;
};

#endif
