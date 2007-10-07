/*
===============================================================================

  FILE:  SPreader_node.h
  
  CONTENTS:
  
    Reads points from Jonathan's node format.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 January 2006 -- added ability to read double precision points
    12 August 2005 -- created after collecting three eggs in the garden
  
===============================================================================
*/
#ifndef SPREADER_NODE_H
#define SPREADER_NODE_H

#include "spreader.h"

#include <stdio.h>

class SPreader_node : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_node functions

  bool open(FILE* fp);
  void force_double();

  SPreader_node();
  ~SPreader_node();

private:
  FILE* file;
  char* line;
  int ncoords;
};

#endif
