/*
===============================================================================

  FILE:  SPreader_spb.cpp
  
  CONTENTS:
  
    see corresponding header file
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    see corresponding header file
  
===============================================================================
*/
#include "spreader_spb.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

bool SPreader_spb::open(FILE* file, bool skip_finalize_header)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPreader_spb\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to binary (untranslated) mode\n");
    }
  }
#endif

  this->file = file;

  read_header();
  read_buffer();

  p_count = 0;

  return true;
}

void SPreader_spb::close()
{
  // close of SPreader interface
  p_count = -1;

  // close of SPreader_spb
  file = 0;

  element_number = 0;
  element_counter = 0;
}

static int swap_endian_int(int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
}

SPevent SPreader_spb::read_event()
{
  if (element_counter < element_number)
  {
    if (element_descriptor & 1) // next element is a point
    {
//*
      if (datatype == SP_DOUBLE)
      {
        if (endian_swap) VecCopy3dv_swap_endian(p_pos_d, &(((double*)(element_buffer))[element_counter*3]));
        else VecCopy3dv(p_pos_d, &(((double*)(element_buffer))[element_counter*3]));
      }
      else
//*/
      {
        if (endian_swap) VecCopy3fv_swap_endian(p_pos_f, &(((float*)(element_buffer))[element_counter*3]));
        else VecCopy3fv(p_pos_f, &(((float*)(element_buffer))[element_counter*3]));
      }
      p_count++;
      element_counter++;
      if (element_counter == element_number)
      {
        read_buffer();
      }
      else
      {
        element_descriptor = element_descriptor >> 1;
      }
      return SP_POINT;
    }
    else // next element is a finalization event
    {
//*
      if (datatype == SP_DOUBLE)
      {
        if (endian_swap) final_idx = swap_endian_int(element_buffer[element_counter*6]);
        else final_idx = element_buffer[element_counter*6];
      }
      else
//*/
      {
        if (endian_swap) final_idx = swap_endian_int(element_buffer[element_counter*3]);
        else final_idx = element_buffer[element_counter*3];
      }
      element_counter++;
      if (element_counter == element_number)
      {
        read_buffer();
      }
      else
      {
        element_descriptor = element_descriptor >> 1;
      }
      return SP_FINALIZED_CELL;
    }
  }
  if (npoints == -1)
  {
    npoints = p_count;
  }
  else
  {
    if (p_count != npoints)
    {
      fprintf(stderr,"ERROR: wrong point count: p_count (%d) != npoints (%d)\n", p_count, npoints);
    }
  }
  return SP_EOF;
}

static unsigned int swap_endian_uint(unsigned int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
}

#define SPB_VERSION 7
#define SPB_LITTLE_ENDIAN 0
#define SPB_BIG_ENDIAN 1

void SPreader_spb::read_header()
{
  int version = fgetc(file);
  // read version
  if (version != SPB_VERSION)
  {
    fprintf(stderr,"ERROR: wrong reader (data is %d but reader is SPB %d)\n", version, SPB_VERSION);
    exit(0);
  }

  // read endianness
#if (defined(i386) || defined(WIN32))   // if little endian machine
  if (fgetc(file) == SPB_LITTLE_ENDIAN) endian_swap = false;
  else endian_swap = true;
#else                                   // else big endian machine
  if (fgetc(file) == SPB_BIG_ENDIAN) endian_swap = false;
  else endian_swap = true;
#endif

  int flag = fgetc(file);

  // which datatype
  datatype = (SPdatatype)(flag & 3);
  switch(datatype)
  {
  case SP_FLOAT:
    fprintf(stderr, "SPdatatype ... SP_FLOAT\n");
    element_size = sizeof(float);
    break;
  case SP_DOUBLE:
    fprintf(stderr, "SPdatatype ... SP_DOUBLE\n");
    element_size = sizeof(double);
    break;
  case SP_INT:
    fprintf(stderr, "SPdatatype ... SP_INT\n");
    element_size = sizeof(int);
    break;
  default:
    fprintf(stderr, "WARNING: unknown SPdatatype %d ... assuming float\n",datatype);
    datatype = SP_FLOAT;
    element_size = sizeof(float);
    break;
  }

  // which finalize method
  finalizemethod = (SPfinalizemethod)(flag >> 2);

  switch(finalizemethod)
  {
    case SP_QUAD_TREE:
//      fprintf(stderr, "INFO: point are finalized with SP_QUAD_TREE\n");
      break;
    case SP_OCT_TREE:
//      fprintf(stderr, "INFO: point are finalized with SP_OCT_TREE\n");
      break;
    case SP_CLARKSON_2D:
//      fprintf(stderr, "INFO: point are finalized with SP_CLARKSON_2D\n");
      break;
    case SP_CLARKSON_3D:
//      fprintf(stderr, "INFO: point are finalized with SP_CLARKSON_3D\n");
      break;
  default:
      fprintf(stderr, "WARNING: SPfinalizemethod %d (maybe legacy point set) ... \n",finalizemethod);
  }

//  if (datatype != SP_FLOAT) fprintf(stderr,"ERROR: wrong reader .. this is for SP_FLOAT\n");
//  if (datatype != SP_DOUBLE) fprintf(stderr,"ERROR: wrong reader .. this is optimized for SP_DOUBLE\n");

  // read comments
  int input;
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) ncomments = swap_endian_int(input);
  else ncomments = input;
  if (ncomments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      fread(&input, sizeof(int), 1, file);
      if (endian_swap) input = swap_endian_int(input);
      comments[i] = (char*)malloc(sizeof(char)*input);
      fread(comments[i], sizeof(char), input, file);
    }
  }
  // read npoints
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  if (input != -1) npoints = input;
  // read bounding box
  if (getc(file))
  {
    if (datatype == SP_FLOAT)
    {
      if (bb_min_f) delete [] bb_min_f;
      if (bb_max_f) delete [] bb_max_f;
      bb_min_f = new float[3];
      bb_max_f = new float[3];
      if (endian_swap)
      {
        float temp[3];
        fread(temp, sizeof(float), 3, file);
        VecCopy3fv_swap_endian(bb_min_f, temp);
        fread(temp, sizeof(float), 3, file);
        VecCopy3fv_swap_endian(bb_max_f, temp);
      }
      else
      {
        fread(bb_min_f, sizeof(float), 3, file);
        fread(bb_max_f, sizeof(float), 3, file);
      }
    }
    else if (datatype == SP_DOUBLE)
    {
      if (bb_min_d) delete [] bb_min_d;
      if (bb_max_d) delete [] bb_max_d;
      bb_min_d = new double[3];
      bb_max_d = new double[3];
      if (endian_swap)
      {
        double temp[3];
        fread(temp, sizeof(double), 3, file);
        VecCopy3dv_swap_endian(bb_min_d, temp);
        fread(temp, sizeof(double), 3, file);
        VecCopy3dv_swap_endian(bb_max_d, temp);
      }
      else
      {
        fread(bb_min_d, sizeof(double), 3, file);
        fread(bb_max_d, sizeof(double), 3, file);
      }
      if (bb_min_f) delete [] bb_min_f;
      if (bb_max_f) delete [] bb_max_f;
      bb_min_f = new float[3];
      bb_max_f = new float[3];
      VecCopy3fv(bb_min_f, bb_min_d);
      VecCopy3fv(bb_max_f, bb_max_d);
    }
    else
    {
      if (bb_min_i) delete [] bb_min_i;
      if (bb_max_i) delete [] bb_max_i;
      bb_min_i = new int[3];
      bb_max_i = new int[3];
      if (endian_swap)
      {
        int temp[3];
        fread(temp, sizeof(int), 3, file);
        VecCopy3iv_swap_endian(bb_min_i, temp);
        fread(temp, sizeof(int), 3, file);
        VecCopy3iv_swap_endian(bb_max_i, temp);
      }
      else
      {
        fread(bb_min_i, sizeof(int), 3, file);
        fread(bb_max_i, sizeof(int), 3, file);
      }
    }
  }
  // allocate buffer
  element_buffer = (int*)malloc(element_size*3*32);
}

void SPreader_spb::read_buffer()
{
  fread(&element_descriptor, sizeof(int), 1, file);
  if (endian_swap) element_descriptor = swap_endian_uint(element_descriptor);
  element_number = fread(element_buffer, element_size, 32*3, file) / 3;
  element_counter = 0;
}

SPreader_spb::SPreader_spb()
{
  // init of SPreader interface
  ncomments = 0;
  comments = 0;

  npoints = -1;
  p_count = -1;

  datatype = SP_VOID;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPreader_spb
  file = 0;
  element_size = -1;
  element_buffer = 0;
}

SPreader_spb::~SPreader_spb()
{
  // clean-up for SPreader interface
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
  }

  if (bb_min_d) delete [] bb_min_d;
  if (bb_max_d) delete [] bb_max_d;
  if (bb_min_f) delete [] bb_min_f;
  if (bb_max_f) delete [] bb_max_f;
  if (bb_min_i) delete [] bb_min_i;
  if (bb_max_i) delete [] bb_max_i;

  // clean-up for SPreader_spb interface
  if (element_buffer) free(element_buffer);
}
