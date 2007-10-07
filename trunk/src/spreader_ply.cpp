/*
===============================================================================

  FILE:  SPreader_ply.cpp
  
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
#include "spreader_ply.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "vec3fv.h"
#include "vec3iv.h"

#include "ply.h"

/* vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;
} Vertex;

static PlyProperty vertex_props[] = { /* list of property information for a vertex */
  {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
};

static PlyFile* in_ply;

bool SPreader_ply::open(FILE* file, bool compute_bounding_box, bool skip_points)
{
  int i;
  int elem_count;
  char *elem_name;
  
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPreader_ply\n");
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

  in_ply = read_ply(file);

  if (in_ply == 0)
  {
    fprintf(stderr, "FATAL ERROR: input PLY file is corrupt\n");
    return false;
  }

  npoints = 0;

  for (i = 0; i < in_ply->num_elem_types; i++)
  {
    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply(in_ply, i, &elem_count);

		if (equal_strings ("vertex", elem_name))
		{
  		npoints = elem_count;
      velem = i;

      if (compute_bounding_box)
      {
        fprintf(stderr, "computing bounding box ...\n");

        setup_property_ply (in_ply, &vertex_props[0]);
			  setup_property_ply (in_ply, &vertex_props[1]);
			  setup_property_ply (in_ply, &vertex_props[2]);

        // alloc bounding box
        if (bb_min_f == 0) bb_min_f = new float[3];
        if (bb_max_f == 0) bb_max_f = new float[3];
        // get current file pointer position
        long here = ftell(in_ply->fp);
        // read first element
  		  get_element_ply (in_ply, (void *)p_pos_f);
        // init bounding box
        VecCopy3fv(bb_min_f, p_pos_f);
        VecCopy3fv(bb_max_f, p_pos_f);
        // process remaining elements
        for (elem_count = 1; elem_count < npoints; elem_count++)
        {
          // read an element
  		    get_element_ply (in_ply, (void *)p_pos_f);
          // update bounding box
          VecUpdateMinMax3fv(bb_min_f, bb_max_f, p_pos_f);
        }
        fprintf(stderr, "bb_min_f %g %g %g\n", bb_min_f[0], bb_min_f[1], bb_min_f[2]);
        fprintf(stderr, "bb_max_f %g %g %g\n", bb_max_f[0], bb_max_f[1], bb_max_f[2]);

        // should we move file pointer back to beginning of point block
        if (!skip_points)
        {
          // this better be a seakable stream ....
          fseek(in_ply->fp, here, SEEK_SET);
        }
      }
      p_count = 0;
      break;
    }
  }
  
  if (npoints == 0)
  {
    fprintf(stderr, "WARNING: no vertices in PLY file\n");
  }

  for (i = 0; i < 3; i++)
  {
    p_pos_f[i] = 0.0f;
  }
  return true;
}

void SPreader_ply::close()
{
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

  npoints = -1;

  p_count = -1;

  close_ply (in_ply);
  free_ply (in_ply);
}

SPevent SPreader_ply::read_event()
{
  if (p_count < npoints)
  {
    if (p_count == 0)
    {
      int elem_count;

      /* prepare to read the vertex elements */
      setup_element_read_ply(in_ply, velem, &elem_count);

      /* set up for getting vertex elements */
			setup_property_ply (in_ply, &vertex_props[0]);
			setup_property_ply (in_ply, &vertex_props[1]);
			setup_property_ply (in_ply, &vertex_props[2]);
    }
		get_element_ply(in_ply, (void *)p_pos_f);
    p_count++;
    return SP_POINT;
  }
  else
  {
    return SP_EOF;
  }
}

SPreader_ply::SPreader_ply()
{
  // init of SMreader interface
  ncomments = 0;
  comments = 0;

  datatype = SP_FLOAT;
  finalizemethod = SP_NONE;

  npoints = -1;
  p_count = -1;

  bb_min_d = 0;
  bb_max_d = 0;
  bb_min_f = 0;
  bb_max_f = 0;
  bb_min_i = 0;
  bb_max_i = 0;

  // init of SPreader_ply
  velem = -1;
}

SPreader_ply::~SPreader_ply()
{
  if (comments)
  {
    for (int i = 0; i < ncomments; i++)
    {
      free(comments[i]);
    }
    free(comments);
  }

  if (bb_min_f) delete [] bb_min_f;
  if (bb_max_f) delete [] bb_max_f;
}
