/*
===============================================================================

  FILE:  SVreader_svb.h
  
  CONTENTS:
   
    Reads a Streaming Volume Mesh from an efficient binary format. 
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2004  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    27 May 2005 -- created after my 450 Euro trip to the BNP in Forbach
  
===============================================================================
*/
#ifndef SVREADER_SVB_H
#define SVREADER_SVB_H

#include "svreader.h"

#include <stdio.h>

class SVreader_svb : public SVreader
{
public:

  // svreader interface function implementations

  void close();

  SVevent read_element();
  SVevent read_event();

  // svreader_svb functions

  bool open(FILE* fp);

  SVreader_svb();
  ~SVreader_svb();

private:
  FILE* file;
  int have_finalized, next_finalized;
  int* finalized_vertices;

  void read_header();
  void read_buffer();

  bool endian_swap;

  int element_number;
  int element_counter;
  unsigned int element_descriptor;
  int* element_buffer;
};

#endif
