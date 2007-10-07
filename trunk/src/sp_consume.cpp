/*
===============================================================================

  FILE:  sp_consume.cpp
  
  CONTENTS:
  
    A program to consume a sample streaming point set that is finalized with
    the CLARKSON_2D method and read the header in a separate file in form
    of a streaming mesh.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2006  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 March 2006 -- created as an example for mario and jack
  
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spreader_spa.h"
#include "spreader_spb.h"

#include "smreader_sma.h"
#include "smreader_smb.h"

#include "vec3fv.h"
#include "vec3iv.h"

static void usage()
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"sp_consume -i points.spa -hi header.spa\n");
  fprintf(stderr,"sp_consume -i points.spb -hi header.spb\n");
  fprintf(stderr,"sp_consume -ispb -hi header.spb\n");
  exit(1);
}

int main(int argc, char *argv[])
{
  int i;
  bool ispa = false;
  bool ispb = false;
  char* file_name_in = 0;
  char* file_name_header_in = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-ispa") == 0)
    {
      ispa = true;
    }
    else if (strcmp(argv[i],"-ispb") == 0)
    {
      ispb = true;
    }
    else if (strcmp(argv[i],"-i") == 0)
    {
      i++;
      file_name_in = argv[i];
    }
    else if (strcmp(argv[i],"-hi") == 0)
    {
      i++;
      file_name_header_in = argv[i];
    }
    else
    {
      usage();
    }
  }

  // open input point file 
  SPreader* spreader;
  FILE* file_in;

  if (file_name_in || ispa || ispb)
  {
    if (file_name_in)
    {
      if (strstr(file_name_in, ".spa") || ispa)
      {
        file_in = fopen(file_name_in, "r");
      }
      else if (strstr(file_name_in, ".spb") || ispb)
      {
        file_in = fopen(file_name_in, "rb");
      }
      else
      {
        if (file_name_in)
        {
          fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name_in);
        }
        else
        {
          fprintf(stderr,"ERROR: no ouput format specified\n");
        }
        exit(0);
      }
      if (file_in == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_in);
        exit(0);
      }
    }
    else
    {
      file_in = stdin;
    }

    if ((file_name_in && strstr(file_name_in, ".spa")) || ispa)
    {
      SPreader_spa* spreader_spa = new SPreader_spa();
      spreader_spa->open(file_in);
      spreader = spreader_spa;
    }
    else if ((file_name_in && strstr(file_name_in, ".spb")) || ispb)
    {
      SPreader_spb* spreader_spb = new SPreader_spb();
      spreader_spb->open(file_in);
      spreader = spreader_spb;
    }
  }

  // make sure we opened a point cloud with the right type of finalization

  if (spreader->finalizemethod != SP_CLARKSON_2D)
  {
    fprintf(stderr,"ERROR: the finalizemethod of the point cloud is not SP_CLARKSON_2D\n");
    exit(0);
  }

  // open output header file ... important ... do not open any earlier to make sure header file is written

  SMreader* smreader;
  FILE* file_header_in;

  if (file_name_header_in)
  {
    if (strstr(file_name_header_in, ".sma"))
    {
      file_header_in = fopen(file_name_header_in, "r");
    }
    else if (strstr(file_name_header_in, ".smb"))
    {
      file_header_in = fopen(file_name_header_in, "rb");
    }
    else
    {
      fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name_header_in);
      exit(0);
    }
    if (file_header_in == 0)
    {
      fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_header_in);
      exit(0);
    }
  }
  else
  {
    fprintf(stderr,"ERROR: no file name given for storing header of finalizemethod \n");
    exit(0);
  }

  if (strstr(file_name_header_in, ".sma"))
  {
    SMreader_sma* smreader_sma = new SMreader_sma();
    smreader_sma->open(file_header_in);
    smreader = smreader_sma;
  }
  else if (strstr(file_name_header_in, ".smb"))
  {
    SMreader_smb* smreader_smb = new SMreader_smb();
    smreader_smb->open(file_header_in);
    smreader = smreader_smb;
  }

  // simple example of reading in a streaming point cloud with SP_CLARKSON_2D finalization

  // *first* we read the header
  
  float* clarkson2d_mesh_verts = 0;
  int clarkson2d_mesh_nverts = smreader->nverts;

  if (clarkson2d_mesh_nverts > 0)
  {
    clarkson2d_mesh_verts = (float*)malloc(clarkson2d_mesh_nverts*sizeof(float)*3);
  }
  else 
  {
    fprintf(stderr,"ERROR: the clarkson finalization mesh does not have any vertices\n");
    exit(0);
  }

  int* clarkson2d_mesh_faces = 0;
  int clarkson2d_mesh_nfaces = smreader->nfaces;

  if (clarkson2d_mesh_nfaces > 0)
  {
    clarkson2d_mesh_faces = (int*)malloc(clarkson2d_mesh_nfaces*sizeof(int)*3);
  }
  else 
  {
    fprintf(stderr,"ERROR: the clarkson finalization mesh does not have any faces\n");
    exit(0);
  }

  SMevent event;

  int clarkson2d_mesh_v_count = 0;
  int clarkson2d_mesh_f_count = 0;

  while (event = smreader->read_element())
  {
    if (event == SM_TRIANGLE)
    {
      VecCopy3iv(&(clarkson2d_mesh_faces[clarkson2d_mesh_f_count*3]), smreader->t_idx);
      clarkson2d_mesh_f_count++;
    }
    else if (event == SM_VERTEX)
    {
      VecCopy3fv(&(clarkson2d_mesh_verts[clarkson2d_mesh_v_count*3]), smreader->v_pos_f);
      clarkson2d_mesh_v_count++;
    }
  }

  smreader->close();
  fclose(file_header_in);
  delete smreader;

  if (clarkson2d_mesh_v_count != clarkson2d_mesh_nverts)
  {
    fprintf(stderr,"WARNING: clarkson2d_mesh_v_count %d is different from clarkson2d_mesh_nverts %d\n", clarkson2d_mesh_v_count, clarkson2d_mesh_nverts);
    exit(0);
  }

  if (clarkson2d_mesh_f_count != clarkson2d_mesh_nfaces)
  {
    fprintf(stderr,"WARNING: clarkson2d_mesh_f_count %d is different from clarkson2d_mesh_nfaces %d\n", clarkson2d_mesh_f_count, clarkson2d_mesh_nfaces);
    exit(0);
  }

  fprintf(stderr,"INFO: read %d vertices and %d triangles from clarkson finalization header\n", clarkson2d_mesh_v_count, clarkson2d_mesh_f_count);

  // *then* we read the points

  int data_npoints = spreader->npoints;

  if (data_npoints > 0)
  {
    fprintf(stderr,"INFO: streaming in %d points ...\n", data_npoints);
  }
  else 
  {
    fprintf(stderr,"ERROR: the point stream does not have any vertices\n");
    exit(0);
  }

  SPevent spevent;

  int data_p_count = 0;
  int data_fin_count = 0;

  while (spevent = spreader->read_event())
  {
    if (spevent == SP_POINT)
    {
      data_p_count++;
    }
    else if (spevent == SP_FINALIZED_CELL)
    {
      data_fin_count++;
    }
  }
  spreader->close();
  if (file_in != stdin) fclose(file_in);
  delete spreader;

  if (data_p_count != data_npoints)
  {
    fprintf(stderr,"WARNING: clarkson2d_mesh_f_count %d is different from clarkson2d_mesh_nfaces %d\n", clarkson2d_mesh_f_count, clarkson2d_mesh_nfaces);
    exit(0);
  }

  fprintf(stderr,"INFO: read %d points and %d finalization events from the point stream\n", data_p_count, data_fin_count);


  return 1;
}
