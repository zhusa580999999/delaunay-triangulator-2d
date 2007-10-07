/*
===============================================================================

  FILE:  SPreader_spa.h
  
  CONTENTS:
  
    Reads points from a simple streaming point ASCII format (SPA).
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 March 2006 -- added finalizemethod to header and support to skip it
    12 January 2006 -- added ability for reading double-precision input
    11 February 2005 -- created after eating a Germknoedel in the Mensa
  
===============================================================================
*/
#ifndef SPREADER_SPA_H
#define SPREADER_SPA_H

#include "spreader.h"

#include <stdio.h>

class SPreader_spa : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // smreader_sma functions

  bool open(FILE* file, bool skip_finalize_header = true);
  void force_datatype(SPdatatype datatype);

  SPreader_spa();
  ~SPreader_spa();

private:
  FILE* file;
  char* line;
};

#endif
