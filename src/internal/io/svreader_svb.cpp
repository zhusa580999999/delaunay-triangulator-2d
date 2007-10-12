/*
===============================================================================

  FILE:  SVreader_svb.cpp
  
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
#include "svreader_svb.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "vecnfv.h"
#include "vecniv.h"

#define SV_VERSION 4 // this is the first SVB implementation

bool SVreader_svb::open(FILE* file)
{
  if (file == 0)
  {
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

  int input = fgetc(file);
  // read version
  if (input != SV_VERSION)
  {
    fprintf(stderr,"ERROR: wrong SVreader (need %d but this is SVreader_svb %d)\n",input,SV_VERSION);
    exit(0);
  }

  read_header();
  read_buffer();

  if (element_descriptor & 1)
  {
    post_order = false;
  }
  else
  {
    post_order = true;
  }

  v_count = 0;
  c_count = 0;

  if (v_pnum_f) v_prop_f = new float[v_pnum_f];
  if (v_pnum_i) v_prop_i = new int[v_pnum_i];
  if (c_pnum_f) c_prop_f = new float[c_pnum_f];
  if (c_pnum_i) c_prop_i = new int[c_pnum_i];

  // currently only TET meshes are supported
  c_type = SV_TET;
  c_idx = new int[4];
  c_final = new bool[4];
  finalized_vertices = new int[4]; 
  return true;
}

void SVreader_svb::close()
{
  // close of SVreader interface
  v_count = -1;
  c_count = -1;

  // close of SVreader_svb
  file = 0;
  have_finalized = 0; next_finalized = 0;

  element_number = 0;
  element_counter = 0;
}

SVevent SVreader_svb::read_element()
{
  if (element_counter < element_number)
  {
    have_finalized = next_finalized = 0;
    if (element_descriptor & 1) // next element is a vertex
    {
      if (endian_swap)
      {
        if (v_pnum_f)
        {
          VecCopyNfv_swap_endian(v_prop_f, (float*)(&element_buffer[element_counter*4]), v_pnum_f);
          if (v_pnum_i)
          {
            VecCopyNiv_swap_endian(v_prop_i, (int*)(&element_buffer[element_counter*4+v_pnum_f]), v_pnum_i);
          }
        }
        else if (v_pnum_i)
        {
          VecCopyNiv_swap_endian(v_prop_i, (int*)(&element_buffer[element_counter*4]), v_pnum_i);
        }
        else
        {
          // strange
        }
      }
      else
      {
        if (v_pnum_f)
        {
          VecCopyNfv(v_prop_f, (float*)(&element_buffer[element_counter*4]), v_pnum_f);
          if (v_pnum_i)
          {
            VecCopyNiv(v_prop_i, (int*)(&element_buffer[element_counter*4+v_pnum_f]), v_pnum_i);
          }
        }
        else if (v_pnum_i)
        {
          VecCopyNiv(v_prop_i, (int*)(&element_buffer[element_counter*4]), v_pnum_i);
        }
        else
        {
          // strange
        }
      }
//      v_idx = v_count;
      if (post_order) {finalized_vertices[have_finalized] = v_count; have_finalized++;}
      v_count++;
      element_counter++;
      if (element_counter == element_number)
      {
        read_buffer();
      }
      else
      {
        element_descriptor = element_descriptor >> 1;
      }
      return SV_VERTEX;
    }
    else // next element is a tetrahedron
    {
      if (endian_swap) VecCopyNiv_swap_endian(c_idx, (int*)(&element_buffer[element_counter*4]), 4);
      else VecCopyNiv(c_idx, (int*)(&element_buffer[element_counter*4]),4);
      c_count++;
      for (int i = 0; i < 4; i++)
      {
        if (c_idx[i] < 0)
        {
          c_idx[i] = v_count+c_idx[i];
          c_final[i] = true;
          finalized_vertices[have_finalized] = c_idx[i];
          have_finalized++;
        }
        else
        {
          c_idx[i] = c_idx[i]-1;
          c_final[i] = false;
        }
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
      return SV_TETRAHEDRON;
    }
  }

  if (nverts != -1 && v_count != nverts)
  {
    fprintf(stderr,"WARNING: wrong vertex count: v_count (%d) != nverts (%d)\n", v_count, nverts);
  }
  nverts = v_count;
  if (ncells != -1 && c_count != ncells)
  {
    fprintf(stderr,"WARNING: wrong face count: c_count (%d) != ncells (%d)\n", c_count, ncells);
  }
  ncells = c_count;
  return SV_EOF;
}

SVevent SVreader_svb::read_event()
{
  if (have_finalized)
  {
    final_idx = finalized_vertices[next_finalized];
    have_finalized--; next_finalized++;
    return SV_FINALIZED;
  }
  else
  {
    return read_element();
  }
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

static unsigned int swap_endian_uint(unsigned int input)
{
  int output;
  ((char*)&output)[0] = ((char*)&input)[3];
  ((char*)&output)[1] = ((char*)&input)[2];
  ((char*)&output)[2] = ((char*)&input)[1];
  ((char*)&output)[3] = ((char*)&input)[0];
  return output;
}

#define SV_LITTLE_ENDIAN 0
#define SV_BIG_ENDIAN 1

void SVreader_svb::read_header()
{
  int input;
  // read endianness
#if (defined(i386) || defined(WIN32))   // if little endian machine
  if (fgetc(file) == SV_LITTLE_ENDIAN) endian_swap = false;
  else endian_swap = true;
#else                                   // else big endian machine
  if (fgetc(file) == SV_BIG_ENDIAN) endian_swap = false;
  else endian_swap = true;
#endif
  // read compression flags (not used yet)
  fgetc(file);
  fgetc(file);
  // read comments
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
  // read nverts
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  if (input != -1) nverts = input;
  // read ncells
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  if (input != -1) ncells = input;
  // read how many floating-point attributes (per vertex)
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  v_pnum_f = input;
  // read how many integer attributes (per vertex)
  fread(&input, sizeof(int), 1, file);
  if (endian_swap) input = swap_endian_int(input);
  v_pnum_i = input;
  // read bounding box of floating-point attributes (per vertex)
  if (v_pnum_f)
  {
    if (getc(file))
    {
      if (v_pmin_f) delete [] v_pmin_f;
      if (v_pmax_f) delete [] v_pmax_f;
      v_pmin_f = new float[v_pnum_f];
      v_pmax_f = new float[v_pnum_f];
      if (endian_swap)
      {
        float temp[4];
        fread(temp, sizeof(float), v_pnum_f, file);
        VecCopyNfv_swap_endian(v_pmin_f, temp, v_pnum_f);
        fread(temp, sizeof(float), v_pnum_f, file);
        VecCopyNfv_swap_endian(v_pmax_f, temp, v_pnum_f);
      }
      else
      {
        fread(v_pmin_f, sizeof(float), v_pnum_f, file);
        fread(v_pmax_f, sizeof(float), v_pnum_f, file);
      }
    }
  }
  // read bounding box of integer attributes (per vertex)
  if (v_pnum_i)
  {
    if (getc(file))
    {
      if (v_pmin_i) delete [] v_pmin_i;
      if (v_pmax_i) delete [] v_pmax_i;
      v_pmin_i = new int[v_pnum_i];
      v_pmax_i = new int[v_pnum_i];
      if (endian_swap)
      {
        int temp[4];
        fread(temp, sizeof(int), v_pnum_i, file);
        VecCopyNiv_swap_endian(v_pmin_i, temp, v_pnum_i);
        fread(temp, sizeof(int), v_pnum_i, file);
        VecCopyNiv_swap_endian(v_pmax_i, temp, v_pnum_i);
      }
      else
      {
        fread(v_pmin_i, sizeof(int), v_pnum_i, file);
        fread(v_pmax_i, sizeof(int), v_pnum_i, file);
      }
    }
  }
}

void SVreader_svb::read_buffer()
{
  fread(&element_descriptor, sizeof(int), 1, file);
  if (endian_swap) element_descriptor = swap_endian_uint(element_descriptor);
  element_number = fread(element_buffer, sizeof(int), 32*4, file) / 4;
  element_counter = 0;
}

SVreader_svb::SVreader_svb()
{
  // init of SVreader interface
  ncomments = 0;
  comments = 0;

  ncells = -1;
  nverts = -1;

  c_count = -1;
  v_count = -1;

  v_pnum_f = 0;
  v_pmin_f = 0;
  v_pmax_f = 0;

  v_pnum_i = 0;
  v_pmin_i = 0;
  v_pmax_i = 0;

  c_pnum_f = 0;
  c_pmin_f = 0;
  c_pmax_f = 0;

  c_pnum_i = 0;
  c_pmin_i = 0;
  c_pmax_i = 0;

  c_idx = 0;
  c_final = 0;

  post_order = false;

  // init of SVreader_svb
  file = 0;
  have_finalized = 0; next_finalized = 0;

  element_buffer = (int*)malloc(sizeof(int)*4*32);
  element_number = 0;
  element_counter = 0;
}

SVreader_svb::~SVreader_svb()
{
  // clean-up for SVreader interface
  if (v_count != -1)
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
  if (v_pmin_f) delete [] v_pmin_f;
  if (v_pmax_f) delete [] v_pmax_f;
  if (v_pmin_i) delete [] v_pmin_i;
  if (v_pmax_i) delete [] v_pmax_i;
  if (c_pmin_f) delete [] c_pmin_f;
  if (c_pmax_f) delete [] c_pmax_f;
  if (c_pmin_i) delete [] c_pmin_i;
  if (c_pmax_i) delete [] c_pmax_i;

  // clean-up for SVwriter_svb interface
  free (element_buffer);
}
