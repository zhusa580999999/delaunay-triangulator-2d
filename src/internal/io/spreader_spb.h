/*
===============================================================================

  FILE:  SPreader_spb.h
  
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
    13 January 2006 -- extended to handle double-precision floats
    15 August 2005 -- created as Don and Carlos are leaving for their SA trip
  
===============================================================================
*/
#ifndef SPREADER_SPB_H
#define SPREADER_SPB_H

#include "spreader.h"

#include <stdio.h>

class SPreader_spb : public SPreader
{
public:

  // spreader interface function implementations

  void close();

  SPevent read_event();

  // spreader_spb functions

  bool open(FILE* fp, bool skip_finalize_header = true);

  SPreader_spb();
  ~SPreader_spb();

private:
  FILE* file;

  void read_header();
  void read_buffer();

  bool endian_swap;

  int element_size;
  int element_number;
  int element_counter;
  unsigned int element_descriptor;
  int* element_buffer;
};

#endif
