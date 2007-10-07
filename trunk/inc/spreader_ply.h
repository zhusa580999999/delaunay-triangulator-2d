/*
===============================================================================

  FILE:  SPreader_ply.h
  
  CONTENTS:
  
    Reads the vertices of a mesh in Stanford's non-streaming PLY format and
    provides access to it in form of a Streaming Point Cloud of maximal width.

    The open function takes one additional flag. It allows to first compute
    the bounding box. Because this is implemented with fseek() the input stream
    must allow backwards seeking. Therefore gipped input streams or stdin input
    cannot make use of this option.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 July 2005 -- created after eating the corn & potatoes stew by carlos.
  
===============================================================================
*/
#ifndef SPREADER_PLY_H
#define SPREADER_PLY_H

#include "spreader.h"

#include <stdio.h>

class SPreader_ply : public SPreader
{
public:

  // smreader interface function implementations

  void close();

  SPevent read_event();

  // smreader_ply functions

  bool open(FILE* fp, bool compute_bounding_box=false, bool skip_points=false);

  SPreader_ply();
  ~SPreader_ply();

private:
  int velem;
};

#endif
