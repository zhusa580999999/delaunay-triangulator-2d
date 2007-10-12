/*
===============================================================================

  FILE:  SVreader_sva.cpp
  
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
#include "svreader_sva.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "vecnfv.h"
#include "vecniv.h"

bool SVreader_sva::open(FILE* file)
{
  if (file == 0)
  {
    return false;
  }

#ifdef _WIN32
  if (file == stdin)
  {
    if(_setmode( _fileno( stdin ), _O_TEXT ) == -1 )
    {
      fprintf(stderr, "ERROR: cannot set stdin to text (translated) mode\n");
    }
  }
#endif

  this->file = file;

  skipped_lines = 0;
  line = (char*)malloc(sizeof(char)*256);
  if (fgets(line, sizeof(char) * 256, file) == 0)
  {
    free(line);
    line = 0;
    return false;
  }

  // look for header information
  char dummy[256];

  while (line && line[0] != 'v' && line[0] != 'c' && line[0] != 'x')
  {
    if (strstr(line, "svtype"))
    {
      char svtype[256];
      sscanf(&(line[1]), "%s %s", dummy, svtype);
      if (strcmp(svtype, "tet") == 0)
      {
        type = SV_TET;
      }
      else
      {
        fprintf(stderr,"FATAL ERROR: input type '%s' not supported\n",svtype);
        exit(0);
      }
    }
    else if (strstr(line, "v_pnum_f"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &v_pnum_f);
    }
    else if (strstr(line, "v_pmin_f"))
    {
      if (v_pmin_f) delete v_pmin_f;
      v_pmin_f = new float[v_pnum_f];
      if (v_pnum_f == 3)
      {
        sscanf(&(line[1]), "%s %f %f %f", dummy, &(v_pmin_f[0]), &(v_pmin_f[1]), &(v_pmin_f[2]));
      }
      else if (v_pnum_f == 4)
      {
        sscanf(&(line[1]), "%s %f %f %f %f", dummy, &(v_pmin_f[0]), &(v_pmin_f[1]), &(v_pmin_f[2]), &(v_pmin_f[3]));
      }
      else if (v_pnum_f == 5)
      {
        sscanf(&(line[1]), "%s %f %f %f %f %f", dummy, &(v_pmin_f[0]), &(v_pmin_f[1]), &(v_pmin_f[2]), &(v_pmin_f[3]), &(v_pmin_f[4]));
      }
    }
    else if (strstr(line, "v_pmax_f"))
    {
      if (v_pmax_f) delete v_pmax_f;
      v_pmax_f = new float[v_pnum_f];
      if (v_pnum_f == 3)
      {
        sscanf(&(line[1]), "%s %f %f %f", dummy, &(v_pmax_f[0]), &(v_pmax_f[1]), &(v_pmax_f[2]));
      }
      else if (v_pnum_f == 4)
      {
        sscanf(&(line[1]), "%s %f %f %f %f", dummy, &(v_pmax_f[0]), &(v_pmax_f[1]), &(v_pmax_f[2]), &(v_pmax_f[3]));
      }
      else if (v_pnum_f == 5)
      {
        sscanf(&(line[1]), "%s %f %f %f %f %f", dummy, &(v_pmax_f[0]), &(v_pmax_f[1]), &(v_pmax_f[2]), &(v_pmax_f[3]), &(v_pmax_f[4]));
      }
    }
    else if (strstr(line, "v_pnum_i"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &v_pnum_i);
    }
    else if (strstr(line, "v_pmin_i"))
    {
      if (v_pmin_i) delete v_pmin_i;
      v_pmin_i = new int[v_pnum_i];
      if (v_pnum_i == 1)
      {
        sscanf(&(line[1]), "%s %d", dummy, &(v_pmin_i[0]));
      }
      else if (v_pnum_i == 2)
      {
        sscanf(&(line[1]), "%s %d %d", dummy, &(v_pmin_i[0]), &(v_pmin_i[1]));
      }
      else if (v_pnum_i == 3)
      {
        sscanf(&(line[1]), "%s %d %d %d", dummy, &(v_pmin_i[0]), &(v_pmin_i[1]), &(v_pmin_i[2]));
      }
      else if (v_pnum_i == 4)
      {
        sscanf(&(line[1]), "%s %d %d %d %d", dummy, &(v_pmin_i[0]), &(v_pmin_i[1]), &(v_pmin_i[2]), &(v_pmin_i[3]));
      }
      else if (v_pnum_i == 5)
      {
        sscanf(&(line[1]), "%s %d %d %d %d %d", dummy, &(v_pmin_i[0]), &(v_pmin_i[1]), &(v_pmin_i[2]), &(v_pmin_i[3]), &(v_pmin_i[4]));
      }
    }
    else if (strstr(line, "v_pmax_i"))
    {
      if (v_pmax_i) delete v_pmax_i;
      v_pmax_i = new int[v_pnum_i];
      if (v_pnum_i == 1)
      {
        sscanf(&(line[1]), "%s %d", dummy, &(v_pmax_i[0]));
      }
      else if (v_pnum_i == 2)
      {
        sscanf(&(line[1]), "%s %d %d", dummy, &(v_pmax_i[0]), &(v_pmax_i[1]));
      }
      else if (v_pnum_i == 3)
      {
        sscanf(&(line[1]), "%s %d %d %d", dummy, &(v_pmax_i[0]), &(v_pmax_i[1]), &(v_pmax_i[2]));
      }
      else if (v_pnum_i == 4)
      {
        sscanf(&(line[1]), "%s %d %d %d %d", dummy, &(v_pmax_i[0]), &(v_pmax_i[1]), &(v_pmax_i[2]), &(v_pmax_i[3]));
      }
      else if (v_pnum_i == 5)
      {
        sscanf(&(line[1]), "%s %d %d %d %d %d", dummy, &(v_pmax_i[0]), &(v_pmax_i[1]), &(v_pmax_i[2]), &(v_pmax_i[3]), &(v_pmax_i[4]));
      }
    }
    else if (strstr(line, "c_pnum_f"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &c_pnum_f);
    }
    else if (strstr(line, "c_pmin_f"))
    {
      if (c_pmin_f) delete c_pmin_f;
      c_pmin_f = new float[c_pnum_f];
      if (c_pnum_f == 1)
      {
        sscanf(&(line[1]), "%s %f", dummy, &(c_pmin_f[0]));
      }
      else if (c_pnum_f == 2)
      {
        sscanf(&(line[1]), "%s %f %f", dummy, &(c_pmin_f[0]), &(c_pmin_f[1]));
      }
      else if (c_pnum_f == 3)
      {
        sscanf(&(line[1]), "%s %f %f %f", dummy, &(c_pmin_f[0]), &(c_pmin_f[1]), &(c_pmin_f[2]));
      }
      else if (c_pnum_f == 4)
      {
        sscanf(&(line[1]), "%s %f %f %f %f", dummy, &(c_pmin_f[0]), &(c_pmin_f[1]), &(c_pmin_f[2]), &(c_pmin_f[3]));
      }
      else if (c_pnum_f == 5)
      {
        sscanf(&(line[1]), "%s %f %f %f %f %f", dummy, &(c_pmin_f[0]), &(c_pmin_f[1]), &(c_pmin_f[2]), &(c_pmin_f[3]), &(c_pmin_f[4]));
      }
    }
    else if (strstr(line, "c_pmax_f"))
    {
      if (c_pmax_f) delete c_pmax_f;
      c_pmax_f = new float[c_pnum_f];
      if (c_pnum_f == 1)
      {
        sscanf(&(line[1]), "%s %f", dummy, &(c_pmax_f[0]));
      }
      else if (c_pnum_f == 2)
      {
        sscanf(&(line[1]), "%s %f %f", dummy, &(c_pmax_f[0]), &(c_pmax_f[1]));
      }
      else if (c_pnum_f == 3)
      {
        sscanf(&(line[1]), "%s %f %f %f", dummy, &(c_pmax_f[0]), &(c_pmax_f[1]), &(c_pmax_f[2]));
      }
      else if (c_pnum_f == 4)
      {
        sscanf(&(line[1]), "%s %f %f %f %f", dummy, &(c_pmax_f[0]), &(c_pmax_f[1]), &(c_pmax_f[2]), &(c_pmax_f[3]));
      }
      else if (c_pnum_f == 5)
      {
        sscanf(&(line[1]), "%s %f %f %f %f %f", dummy, &(c_pmax_f[0]), &(c_pmax_f[1]), &(c_pmax_f[2]), &(c_pmax_f[3]), &(c_pmax_f[4]));
      }
    }
    else if (strstr(line, "c_pnum_i"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &c_pnum_i);
    }
    else if (strstr(line, "c_pmin_i"))
    {
      if (c_pmin_i) delete c_pmin_i;
      c_pmin_i = new int[c_pnum_i];
      if (c_pnum_i == 1)
      {
        sscanf(&(line[1]), "%s %d", dummy, &(c_pmin_i[0]));
      }
      else if (c_pnum_i == 2)
      {
        sscanf(&(line[1]), "%s %d %d", dummy, &(c_pmin_i[0]), &(c_pmin_i[1]));
      }
      else if (c_pnum_i == 3)
      {
        sscanf(&(line[1]), "%s %d %d %d", dummy, &(c_pmin_i[0]), &(c_pmin_i[1]), &(c_pmin_i[2]));
      }
      else if (c_pnum_i == 4)
      {
        sscanf(&(line[1]), "%s %d %d %d %d", dummy, &(c_pmin_i[0]), &(c_pmin_i[1]), &(c_pmin_i[2]), &(c_pmin_i[3]));
      }
      else if (c_pnum_i == 5)
      {
        sscanf(&(line[1]), "%s %d %d %d %d %d", dummy, &(c_pmin_i[0]), &(c_pmin_i[1]), &(c_pmin_i[2]), &(c_pmin_i[3]), &(c_pmin_i[4]));
      }
    }
    else if (strstr(line, "c_pmax_i"))
    {
      if (c_pmax_i) delete c_pmax_i;
      c_pmax_i = new int[c_pnum_i];
      if (c_pnum_i == 1)
      {
        sscanf(&(line[1]), "%s %d", dummy, &(c_pmax_i[0]));
      }
      else if (c_pnum_i == 3)
      {
        sscanf(&(line[1]), "%s %d %d", dummy, &(c_pmax_i[0]), &(c_pmax_i[1]));
      }
      else if (c_pnum_i == 3)
      {
        sscanf(&(line[1]), "%s %d %d %d", dummy, &(c_pmax_i[0]), &(c_pmax_i[1]), &(c_pmax_i[2]));
      }
      else if (c_pnum_i == 4)
      {
        sscanf(&(line[1]), "%s %d %d %d %d", dummy, &(c_pmax_i[0]), &(c_pmax_i[1]), &(c_pmax_i[2]), &(c_pmax_i[3]));
      }
      else if (c_pnum_i == 5)
      {
        sscanf(&(line[1]), "%s %d %d %d %d %d", dummy, &(c_pmax_i[0]), &(c_pmax_i[1]), &(c_pmax_i[2]), &(c_pmax_i[3]), &(c_pmax_i[4]));
      }
    }
    else if (strstr(line, "bb_min")) // for claudio's *.sma variation
    {
      v_pnum_f = 3;
      if (v_pmin_f) delete v_pmin_f;
      v_pmin_f = new float[v_pnum_f];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(v_pmin_f[0]), &(v_pmin_f[1]), &(v_pmin_f[2]));
    }
    else if (strstr(line, "bb_max")) // for claudio's *.sma variation
    {
      v_pnum_f = 3;
      if (v_pmax_f) delete v_pmax_f;
      v_pmax_f = new float[v_pnum_f];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(v_pmax_f[0]), &(v_pmax_f[1]), &(v_pmax_f[2]));
    }
    else if (strstr(line, "nverts"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &nverts);
    }
    else if (strstr(line, "ncells"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &ncells);
    }
    else if (line[0] == '#')
    {
      if (sscanf(&(line[1]), "%s", dummy))
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
        comments[ncomments] = strdup(dummy);
        ncomments++;
      }
    }
    else if (line[0] == 'p')
    {
      fprintf(stderr,"WARNING: input seems to be a point cloud\n");
    }
    else if (line[0] == 'f')
    {
      fprintf(stderr,"WARNING: input seems to be a surface mesh\n");
    }
    else
    {
      int j = (int)strlen(line)-1;
      if (j > 0)
      {
        for (int i = 0; i < j; i++)
        {
          if (!isprint(line[i]))
          {
            fprintf(stderr,"WARNING: input file seems to be binary\n");
            break;
          }
        }
        line[(10<j?10:j)] = '\0';
        if (skipped_lines < 5)
        {
          fprintf(stderr,"WARNING: skipping line '%s...'.\n",line);
        }
        else if (skipped_lines == 5)
        {
          fprintf(stderr,"WARNING: skipping more than 5 lines. being quiet.\n");
        }
        skipped_lines++;
      }
    }
    if (fgets(line, sizeof(char) * 256, file) == 0)
    {
      free(line);
      line = 0;
      return false;
    }
  }
  if (line[0] == 'v')
  {
    post_order = false;
  }
  else if (line[0] == 'c')
  {
    post_order = true;
  }
  else
  {
    fprintf(stderr,"WARNING: cannot determine order. assuming pre-order SM.\n");
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

void SVreader_sva::close()
{
  // close of SMreader interface
  v_count = -1;
  c_count = -1;

  delete [] c_idx; c_idx = 0;
  delete [] c_final; c_final = 0;

  // close of SVreader_sva
  if (skipped_lines) fprintf(stderr,"WARNING: skipped %d lines.\n",skipped_lines);
  file = 0;
  skipped_lines = 0;
  if (line)
  {
    free(line);
    line = 0;
  }
  have_finalized = 0; next_finalized = 0;
  delete [] finalized_vertices; finalized_vertices = 0;
}

SVevent SVreader_sva::read_element()
{
  have_finalized = next_finalized = 0;
  while (line)
  {
    if ((line[0] == 'v') && (line[1] == ' '))
    {
      if (v_pnum_f == 3)
      {
        if (v_pnum_i == 0)
        {
          sscanf(&(line[1]), "%f %f %f", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]));
        }
        else if (v_pnum_i == 1)
        {
          sscanf(&(line[1]), "%f %f %f %d", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]), &(v_prop_i[0]));
        }
        else if (v_pnum_i == 2)
        {
          sscanf(&(line[1]), "%f %f %f %d %d", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]), &(v_prop_i[0]), &(v_prop_i[1]));
        }
      }
      else if (v_pnum_f == 4)
      {
        if (v_pnum_i == 0)
        {
          sscanf(&(line[1]), "%f %f %f %f", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]), &(v_prop_f[3]));
        }
        else if (v_pnum_i == 1)
        {
          sscanf(&(line[1]), "%f %f %f %f %d", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]), &(v_prop_f[3]), &(v_prop_i[0]));
        }
        else if (v_pnum_i == 2)
        {
          sscanf(&(line[1]), "%f %f %f %f %d %d", &(v_prop_f[0]), &(v_prop_f[1]), &(v_prop_f[2]), &(v_prop_f[3]), &(v_prop_i[0]), &(v_prop_i[1]));
        }
      }
      if (post_order) {finalized_vertices[have_finalized] = v_count; have_finalized++;}
      v_count++;
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SV_VERTEX;
    }
    else if ((line[0] == 'c') && (line[1] == ' '))
    {
      if (c_pnum_f == 0)
      {
        if (c_pnum_i == 0)
        {
          sscanf(&(line[1]), "%d %d %d %d", &(c_idx[0]), &(c_idx[1]), &(c_idx[2]), &(c_idx[3]));
        }
        else if (c_pnum_i == 1)
        {
          sscanf(&(line[1]), "%d %d %d %d %d", &(c_idx[0]), &(c_idx[1]), &(c_idx[2]), &(c_idx[3]), &(c_prop_i[0]));
        }
        else if (c_pnum_i == 2)
        {
          sscanf(&(line[1]), "%d %d %d %d %d %d", &(c_idx[0]), &(c_idx[1]), &(c_idx[2]), &(c_idx[3]), &(c_prop_i[0]), &(c_prop_i[1]));
        }
      }
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
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SV_TETRAHEDRON;
    }
    else if ((line[0] == 'x') && (line[1] == ' '))
    {
      sscanf(&(line[1]), "%d", &(final_idx));
      if (final_idx < 0)
      {
        final_idx = v_count+final_idx;
      }
      else
      {
        final_idx = final_idx - 1;
      }
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SV_FINALIZED;
    }
    else if (line[0] == '#')
    {
      // comments in the body are silently ignored
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
    else
    {
      int j = (int)strlen(line)-1;
      if (j)
      {
        for (int i = 0; i < j; i++)
        {
          if (!isprint(line[i]))
          {
            fprintf(stderr,"WARNING: input file seems to be binary\n");
            break;
          }
        }
        line[(10<j?10:j)] = '\0';
        if (skipped_lines < 5)
        {
          fprintf(stderr,"WARNING: skipping line '%s...'.\n",line);
        }
        else if (skipped_lines == 5)
        {
          fprintf(stderr,"WARNING: skipping more than 5 lines. being quiet.\n");
        }
        skipped_lines++;
      }
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
  }

  if (nverts == -1)
  {
    nverts = v_count;
  }
  else
  {
    if (v_count != nverts)
    {
      fprintf(stderr,"WARNING: wrong vertex count: v_count (%d) != nverts (%d)\n", v_count, nverts);
    }
  }
  if (ncells == -1)
  {
    ncells = c_count;
  }
  else
  {
    if (c_count != ncells)
    {
      fprintf(stderr,"WARNING: wrong cell count: c_count (%d) != ncells (%d)\n", c_count, ncells);
    }
  }
  return SV_EOF;
}

SVevent SVreader_sva::read_event()
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

SVreader_sva::SVreader_sva()
{
  // init of SMreader interface
  ncomments = 0;
  comments = 0;

  nverts = -1;
  ncells = -1;

  v_count = -1;
  c_count = -1;

  v_pnum_f = 0;
  v_pmin_f = 0;
  v_pmax_f = 0;

  v_pnum_i = 0;
  v_pmin_i = 0;
  v_pmax_i = 0;

  v_prop_f = 0;
  v_prop_i = 0;

  c_pnum_f = 0;
  c_pmin_f = 0;
  c_pmax_f = 0;

  c_pnum_i = 0;
  c_pmin_i = 0;
  c_pmax_i = 0;

  c_prop_f = 0;
  c_prop_i = 0;

  c_idx = 0;
  c_final = 0;

  post_order = false;

  // init of SVreader_sva
  file = 0;
  line = 0;
  have_finalized = 0; next_finalized = 0;
  finalized_vertices = 0;
}

SVreader_sva::~SVreader_sva()
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
}
