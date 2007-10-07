/*
===============================================================================

  FILE:  SPconverter.h
  
  CONTENTS:
  
    Reads Streaming Points from a Streaming Mesh. While this is a reduction of
    functionality it is a useful tool to have for experiments.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    03 October 2005 -- added svreader on the most windy day so far in Berkeley
    18 February 2005 -- created just before going to see a movie with oli
  
===============================================================================
*/
#ifndef SPCONVERTER_H
#define SPCONVERTER_H

#include "spreader.h"
#include "smreader.h"
#include "svreader.h"

class SPconverter : public SPreader
{
public:

  // smreader interface function implementations

  void close();

  SPevent read_event();

  // smreader_converter functions

  bool open(SMreader* smreader);
  bool open(SVreader* svreader);

  SPconverter();
  ~SPconverter();

private:
  SMreader* smreader;
  SVreader* svreader;
};

#endif
