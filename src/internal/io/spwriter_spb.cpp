/*
===============================================================================

  FILE:  SPwriter_spb.cpp
  
  CONTENTS:

    see corresponding header file

  PROGRAMMERS:

    martin isenburg@cs.unc.edu

  COPYRIGHT:

    copyright (C) 2005  martin isenburg@cs.unc.edu

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  CHANGE HISTORY:

    see corresponding header file

===============================================================================
*/
#include "spwriter_spb.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

bool SPwriter_spb::open(FILE* file)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPwriter_spb\n");
    return false;
  }

#ifdef _WIN32
  if (file == stdout)
  {
    if(_setmode( _fileno( stdout ), _O_BINARY ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdout to binary (untranslated) mode\n");
    }
  }
#endif

  this->file = file;
  ncomments = 0;
  p_count = 0;
  datatype = SP_FLOAT; // default data type

  element_number = 0;
  element_descriptor = 0;

  return true;
}

void SPwriter_spb::add_comment(const char* comment)
{
  if (comments == 0)
  {
    ncomments = 0;
    comments = (char**)malloc(sizeof(char*)*10);
    comments[9] = (char*)-1;
  }
  else if (comments[ncomments] == (char*)-1)
  {
    comments = (char**)realloc(comments,sizeof(char*)*ncomments*2);
    comments[ncomments*2-1] = (char*)-1;
  }
  comments[ncomments] = strdup(comment);
  ncomments++;
}

void SPwriter_spb::set_datatype(SPdatatype datatype)
{
  this->datatype = datatype;
}

void SPwriter_spb::set_finalizemethod(SPfinalizemethod finalizemethod)
{
  this->finalizemethod = finalizemethod;
}

void SPwriter_spb::set_npoints(int npoints)
{
  this->npoints = npoints;
}

void SPwriter_spb::set_boundingbox(const double* bb_min_d, const double* bb_max_d)
{
  if (this->bb_min_d == 0) this->bb_min_d = new double[3];
  if (this->bb_max_d == 0) this->bb_max_d = new double[3];
  VecCopy3dv(this->bb_min_d, bb_min_d);
  VecCopy3dv(this->bb_max_d, bb_max_d);
  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];
  VecCopy3fv(this->bb_min_f, bb_min_d);
  VecCopy3fv(this->bb_max_f, bb_max_d);
}

void SPwriter_spb::set_boundingbox(const float* bb_min_f, const float* bb_max_f)
{
  if (this->bb_min_f == 0) this->bb_min_f = new float[3];
  if (this->bb_max_f == 0) this->bb_max_f = new float[3];
  VecCopy3fv(this->bb_min_f, bb_min_f);
  VecCopy3fv(this->bb_max_f, bb_max_f);
}

void SPwriter_spb::set_boundingbox(const int* bb_min_i, const int* bb_max_i)
{
  if (this->bb_min_i == 0) this->bb_min_i = new int[3];
  if (this->bb_max_i == 0) this->bb_max_i = new int[3];
  VecCopy3iv(this->bb_min_i, bb_min_i);
  VecCopy3iv(this->bb_max_i, bb_max_i);
}

void SPwriter_spb::set_endianness(bool big_endian)
{
#if (defined(i386) || defined(WIN32))   // if little endian machine
  endian_swap = big_endian;
#else                                   // else big endian machine
  endian_swap = !big_endian;
#endif
}

#define SPB_VERSION 7
#define SPB_LITTLE_ENDIAN 0
#define SPB_BIG_ENDIAN 1

static int swap_endian_int(int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
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

void SPwriter_spb::write_header()
{
  // write version
  fputc(SPB_VERSION, file);

  // write endianness
#if (defined(i386) || defined(WIN32))   // if little endian machine
  if (endian_swap) fputc(SPB_BIG_ENDIAN, file);
  else fputc(SPB_LITTLE_ENDIAN, file);
#else                                    // else big endian machine
  if (endian_swap) fputc(SPB_LITTLE_ENDIAN, file);
  else fputc(SPB_BIG_ENDIAN, file);
#endif

  char flag;

  // which datatype
  switch (datatype)
  {
  case SP_FLOAT:
    flag = SP_FLOAT;
    element_size = sizeof(float);
    break;
  case SP_DOUBLE:
    flag = SP_DOUBLE;
    element_size = sizeof(double);
    break;
  case SP_INT:
    flag = SP_INT;
    element_size = sizeof(int);
    break;
  default:
    fprintf(stderr, "WARNING: unknown SPdatatype %d ... assuming float\n",datatype);
    datatype = SP_FLOAT;
    flag = SP_FLOAT;
    element_size = sizeof(float);
    break;
  }

  // which finalizemethod
  flag = flag | (finalizemethod << 2);

  putc(flag, file);

  // write comments
  int output;
  if (endian_swap) output = swap_endian_int(ncomments);
  else output = ncomments;
  fwrite(&output, sizeof(int), 1, file);
  if (ncomments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      if (endian_swap) output = swap_endian_int(strlen(comments[i]));
      else output = strlen(comments[i]);
      fwrite(&output, sizeof(int), 1, file);
      fwrite(comments[i], sizeof(char), output, file);
    }
  }

  // write npoints
  if (endian_swap) output = swap_endian_int(npoints);
  else output = npoints;
  fwrite(&output, sizeof(int), 1, file);

  // write bounding box
  if (datatype == SP_FLOAT && bb_min_f && bb_max_f)
  {
    fputc(1, file);
    if (endian_swap)
    {
      float temp[3];
      VecCopy3fv_swap_endian(temp, bb_min_f);
      fwrite(temp, sizeof(float), 3, file);
      VecCopy3fv_swap_endian(temp, bb_max_f);
      fwrite(temp, sizeof(float), 3, file);
    }
    else
    {
      fwrite(bb_min_f, sizeof(float), 3, file);
      fwrite(bb_max_f, sizeof(float), 3, file);
    }
  }
  else if (datatype == SP_DOUBLE && bb_min_d && bb_max_d)
  {
    fputc(1, file);
    if (endian_swap)
    {
      double temp[3];
      VecCopy3dv_swap_endian(temp, bb_min_d);
      fwrite(temp, sizeof(double), 3, file);
      VecCopy3dv_swap_endian(temp, bb_max_d);
      fwrite(temp, sizeof(double), 3, file);
    }
    else
    {
      fwrite(bb_min_d, sizeof(double), 3, file);
      fwrite(bb_max_d, sizeof(double), 3, file);
    }
  }
  else if (datatype == SP_INT && bb_min_i && bb_max_i)
  {
    fputc(1, file);
    if (endian_swap)
    {
      int temp[3];
      VecCopy3iv_swap_endian(temp, bb_min_i);
      fwrite(temp, sizeof(int), 3, file);
      VecCopy3iv_swap_endian(temp, bb_max_i);
      fwrite(temp, sizeof(int), 3, file);
    }
    else
    {
      fwrite(bb_min_i, sizeof(float), 3, file);
      fwrite(bb_max_i, sizeof(float), 3, file);
    }
  }
  else
  {
    fputc(0, file);
  }

  // allocate buffer
  element_buffer = (int*)malloc(element_size*3*32);
}

void SPwriter_spb::write_point(const double* p_pos_d)
{
  if (endian_swap) VecCopy3dv_swap_endian(&(((double*)element_buffer)[element_number*3]), p_pos_d);
  else VecCopy3dv(&(((double*)element_buffer)[element_number*3]), p_pos_d);
  element_descriptor = 0x80000000 | (element_descriptor >> 1);
  element_number++;

  if (element_number == 32) write_buffer();

  p_count++;
}

void SPwriter_spb::write_point(const float* p_pos_f)
{
  if (endian_swap) VecCopy3fv_swap_endian(&(((float*)element_buffer)[element_number*3]), p_pos_f);
  else VecCopy3fv(&(((float*)element_buffer)[element_number*3]), p_pos_f);
  element_descriptor = 0x80000000 | (element_descriptor >> 1);
  element_number++;

  if (element_number == 32) write_buffer();

  p_count++;
}

void SPwriter_spb::write_point(const int* p_pos_i)
{
  if (endian_swap) VecCopy3iv_swap_endian(&(((int*)element_buffer)[element_number*3]), p_pos_i);
  else VecCopy3iv(&(((int*)element_buffer)[element_number*3]), p_pos_i);
  element_descriptor = 0x80000000 | (element_descriptor >> 1);
  element_number++;

  if (element_number == 32) write_buffer();

  p_count++;
}

void SPwriter_spb::write_finalize_cell(int idx)
{
  if (datatype == SP_DOUBLE)
  {
    if (endian_swap)
    {
      element_buffer[element_number*6] = swap_endian_int(idx);
    }
    else
    {
      element_buffer[element_number*6] = idx;
    }
    element_buffer[element_number*6+1] = 0;
    element_buffer[element_number*6+2] = 0;
    element_buffer[element_number*6+3] = 0;
    element_buffer[element_number*6+4] = 0;
    element_buffer[element_number*6+5] = 0;
  }
  else
  {
    if (endian_swap)
    {
      element_buffer[element_number*3] = swap_endian_int(idx);
    }
    else
    {
      element_buffer[element_number*3] = idx;
    }
    element_buffer[element_number*3+1] = 0;
    element_buffer[element_number*3+2] = 0;
  }
  element_descriptor = (element_descriptor >> 1);
  element_number++;
  if (element_number == 32) write_buffer();
}

void SPwriter_spb::close()
{
  write_buffer_remaining();

  file = 0;

  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
    ncomments = 0;
    comments = 0;
  }

  if (npoints != -1) if (npoints != p_count)  fprintf(stderr,"WARNING: set npoints %d but p_count %d\n",npoints,p_count);

  p_count = -1;
}

void SPwriter_spb::write_buffer()
{
  if (endian_swap) element_descriptor = swap_endian_uint(element_descriptor);
  fwrite(&element_descriptor, sizeof(unsigned int), 1, file);
  element_descriptor = 0;
  fwrite(element_buffer, element_size, 32*3, file);
  element_number = 0;
}

void SPwriter_spb::write_buffer_remaining()
{
  element_descriptor = element_descriptor >> (32 - element_number);
  if (endian_swap) element_descriptor = swap_endian_uint(element_descriptor);
  fwrite(&element_descriptor, sizeof(unsigned int), 1, file);
  element_descriptor = 0;
  fwrite(element_buffer, element_size, element_number*3, file);
  element_number = 0;
}

SPwriter_spb::SPwriter_spb()
{
  // init of SPwriter interface
  ncomments = 0;
  comments = 0;

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

  // init of SPwriter_spb interface
  file = 0;
  element_size = -1;
  endian_swap = false;
  element_buffer = 0;
}

SPwriter_spb::~SPwriter_spb()
{
  // clean-up for SPwriter interface
  if (p_count != -1)
  {
    close(); // user must have forgotten to close the mesh
  }
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

  // clean-up for SPwriter_spb interface
  if (element_buffer) free(element_buffer);
}
