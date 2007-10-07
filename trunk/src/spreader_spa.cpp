/*
===============================================================================

  FILE:  SPreader_spa.cpp
  
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
#include "spreader_spa.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "vec3dv.h"
#include "vec3fv.h"
#include "vec3iv.h"

bool SPreader_spa::open(FILE* file, bool skip_finalize_header)
{
  if (file == 0)
  {
    fprintf(stderr, "ERROR: zero file pointer not supported by SPreader_spa\n");
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

  line = (char*)malloc(sizeof(char)*256);
  if (fgets(line, sizeof(char) * 256, file) == 0)
  {
    free(line);
    return false;
  }

  // look for header information
  char dummy[256];

  while (line && (line[0] != 'v') && (line[0] != 'p') && (line[0] != 'x'))
  {
    if (strstr(line, "npoints"))
    {
      sscanf(&(line[1]), "%s %d", dummy, &npoints);
    }
    else if (strstr(line, "bb_min"))
    {
      if (bb_min_f == 0) bb_min_f = new float[3];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(bb_min_f[0]), &(bb_min_f[1]), &(bb_min_f[2]));
    }
    else if (strstr(line, "bb_max"))
    {
      if (bb_max_f == 0) bb_max_f = new float[3];
      sscanf(&(line[1]), "%s %f %f %f", dummy, &(bb_max_f[0]), &(bb_max_f[1]), &(bb_max_f[2]));
    }
    else if (strstr(line, "bb_min_i"))
    {
      if (bb_min_f == 0) bb_min_f = new float[3];
      sscanf(&(line[1]), "%s %d %d %d", dummy, &(bb_min_i[0]), &(bb_min_i[1]), &(bb_min_i[2]));
    }
    else if (strstr(line, "bb_max_i"))
    {
      if (bb_max_f == 0) bb_max_f = new float[3];
      sscanf(&(line[1]), "%s %d %d %d", dummy, &(bb_max_i[0]), &(bb_max_i[1]), &(bb_max_i[2]));
    }
    else if (strstr(line, "datatype"))
    {
      if (strstr(line, "SP_FLOAT"))
      {
        datatype = SP_FLOAT;
      }
      else if (strstr(line, "SP_DOUBLE"))
      {
        datatype = SP_DOUBLE;
      }
      else if (strstr(line, "SP_INT"))
      {
        datatype = SP_INT;
      }
      else
      {
        fprintf(stderr,"WARNING: unknown datatype specified in input line:\n");
        fprintf(stderr,line);
      }
    }
    else if (strstr(line, "finalizemethod"))
    {
      if (strstr(line, "SP_QUAD_TREE"))
      {
        finalizemethod = SP_QUAD_TREE;
      }
      else if (strstr(line, "SP_OCT_TREE"))
      {
        finalizemethod = SP_OCT_TREE;
      }
      else if (strstr(line, "SP_CLARKSON_2D"))
      {
        finalizemethod = SP_CLARKSON_2D;
      }
      else if (strstr(line, "SP_CLARKSON_3D"))
      {
        finalizemethod = SP_CLARKSON_3D;
      }
      else if (strstr(line, "SP_NONE"))
      {
        finalizemethod = SP_NONE;
      }
      else
      {
        fprintf(stderr,"WARNING: unknown finalizemethod specified in input line:\n");
        fprintf(stderr,line);
      }
    }
    else if (line[0] == '#')
    {
      if (sscanf(&(line[1]), "%s", dummy) == 1)
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
    else
    {
      float help;
      if ((strlen(line) >= 5) && (sscanf(&(line[0]), "%f %f %f", &(help), &(help), &(help)) == 3))
      {
        break;
      }
      if ((strlen(line) >= 5) && (sscanf(&(line[0]), "%f,%f,%f", &(help), &(help), &(help)) == 3))
      {
        break;
      }
    }

    if (fgets(line, sizeof(char) * 256, file) == 0)
    {
      free(line);
      return false;
    }
  }

  p_count = 0;

  return true;
}

void SPreader_spa::force_datatype(SPdatatype datatype)
{
  this->datatype = datatype;
}

void SPreader_spa::close()
{
  p_count = -1;
}

SPevent SPreader_spa::read_event()
{
  while (line)
  {
    if (((line[0] == 'p') || (line[0] == 'v')) && (line[1] == ' '))
    {
      if (datatype == SP_FLOAT)
      {
        sscanf(&(line[1]), "%f %f %f", &(p_pos_f[0]), &(p_pos_f[1]), &(p_pos_f[2]));
      }
      else if (datatype == SP_DOUBLE)
      {
        sscanf(&(line[1]), "%lf %lf %lf", &(p_pos_d[0]), &(p_pos_d[1]), &(p_pos_d[2]));
      }
      else if (datatype == SP_INT)
      {
        sscanf(&(line[1]), "%d %d %d", &(p_pos_i[0]), &(p_pos_i[1]), &(p_pos_i[2]));
      }
      p_count++;
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SP_POINT;
    }
    else if ((line[0] == 'x') && (line[1] == ' '))
    {
      if ((line[2] == 'c') && (line[3] == 'e')  && (line[4] == 'l') && (line[5] == 'l'))
      {
        sscanf(&(line[7]), "%d", &final_idx);
        if (fgets(line, sizeof(char) * 256, file) == 0)
        {
          free(line);
          line = 0;
        }
        return SP_FINALIZED_CELL;
      }
      else
      {
        // skipping an unknown finalization event
        if (fgets(line, sizeof(char) * 256, file) == 0)
        {
          free(line);
          line = 0;
        }
      }
    }
    else if (line[0] == 'f')
    {
      // skipping a face
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
    else if (line[0] == 'c')
    {
      // skipping a cell
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
    else if (line[0] == '#')
    {
      // skipping a comment
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
    }
    else if (strlen(line) >= 5 && sscanf(&(line[0]), "%f %f %f", &(p_pos_f[0]), &(p_pos_f[1]), &(p_pos_f[2])) == 3)
    {
      p_count++;
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SP_POINT;
    }
    else if (strlen(line) >= 5 && sscanf(&(line[0]), "%f,%f,%f", &(p_pos_f[0]), &(p_pos_f[1]), &(p_pos_f[2])) == 3)
    {
      p_count++;
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
      return SP_POINT;
    }
    else
    {
      for (int i = 0; i < (int)strlen(line); i++)
      {
        if (!isprint(line[i]))
        {
          fprintf(stderr,"FATAL ERROR: input file is binary\n");
        }
      }
      line[10] = '\0';
      fprintf(stderr,"WARNING: skipping line '%s...'",line);
      if (fgets(line, sizeof(char) * 256, file) == 0)
      {
        free(line);
        line = 0;
      }
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

SPreader_spa::SPreader_spa()
{
  // init of SPreader interface
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

  // init of SPreader_spa
  file = 0;
  line = 0;
}

SPreader_spa::~SPreader_spa()
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
  if (bb_min_i) delete [] bb_min_i;
  if (bb_max_i) delete [] bb_max_i;
}
