/*
===============================================================================

  FILE:  sp_create.cpp
  
  CONTENTS:
  
    A program to create a sample streaming point set that is finalized with
    the CLARKSON_2D method and stores the header in a separate file in form
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

#include "smwriter_sma.h"
#include "smwriter_smb.h"

#include "spwriter_spa.h"
#include "spwriter_spb.h"

#include "vec3iv.h"
#include "vec3fv.h"

static void usage()
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"sp_create -o points.spa -ho header.spa\n");
  fprintf(stderr,"sp_create -o points.spb -ho header.spb\n");
  fprintf(stderr,"sp_create -ospb -ho header.spb\n");
  exit(1);
}

int clarkson2d_mesh_nverts = 12;

float clarkson2d_mesh_verts[12][3] = {
  {0.0f, 0.0f, -99999.9f},  // infinite vertex
  {89.4063f, 114.5f, 197.2f},
  {90.0625f, 46.5f, 196.93f},
  {95.6563f, 23.0f, 195.99f},
  {75.1875f, 21.5f, 197.87f},
  {56.5313f, 23.5f, 198.07f},
  {41.1563f, 25.5f, 195.31f},
  {68.9688f, 146.5f, 202.12f},
  {9.6875f, 153.0f, 194.5f},
  {60.6875f, 81.0f, 195.63f},
  {75.75f, 136.5f, 195.88f},
  {68.1875f, 81.5f, 196.55f}
};

int clarkson2d_mesh_nfaces = 20;

int clarkson2d_mesh_faces[20][3] = {
  {0, 1, 2},
  {3, 0, 2},
  {0, 3, 4},
  {0, 4, 5},
  {6, 0, 5},
  {3, 2, 4},
  {7, 8, 9},
  {10, 7, 9},
  {1, 10, 9},
  {1, 11, 2},
  {11, 1, 9},
  {11, 9, 2},
  {2, 9, 5},
  {4, 2, 5},
  {9, 6, 5},
  {8, 6, 9},
  {6, 8, 0},
  {8, 7, 0},
  {7, 10, 0},
  {10, 1, 0}
};

int data_npoints = 40;

float data_points[40][3] = {
  {11.5313f, 41.5f, 193.73f},
  {1.75f, 112.5f, 192.63f},
  {19.7813f, 91.0f, 193.98f},
  {5.84375f, 128.0f, 195.79f},
  {9.6875f, 153.0f, 194.5f},
  {40.9688f, 17.5f, 195.33f},
  {41.1563f, 25.5f, 195.31f},
  {45.0f, 25.5f, 195.44f},
  {36.25f, 115.5f, 194.58f},
  {116.406f, 7.0f, 195.77f},
  {89.4063f, 114.5f, 197.2f},
  {90.0625f, 46.5f, 196.93f},
  {95.6563f, 23.0f, 195.99f},
  {75.1875f, 21.5f, 197.87f},
  {56.5313f, 23.5f, 198.07f},
  {41.1563f, 25.5f, 195.31f},
  {68.9688f, 146.5f, 202.12f},
  {9.6875f, 153.0f, 194.5f},
  {60.6875f, 81.0f, 195.63f},
  {75.75f, 136.5f, 195.88f},
  {55.6875f, 4.1f, 199.83f},
  {56.5313f, 23.5f, 198.07f},
  {68.6875f, 0.5f, 196.12f},
  {50.5f, 86.5f, 195.3f},
  {55.625f, 63.3f, 196.08f},
  {60.6875f, 81.2f, 195.63f},
  {68.1875f, 81.5f, 196.55f},
  {62.6563f, 136.1f, 195.84f},
  {70.6563f, 9.5f, 196.05f},
  {75.1875f, 21.5f, 197.87f},
  {83.1563f, 40.8f, 196.43f},
  {90.1f, 42.5f, 196.65f},
  {90.0625f, 46.5f, 196.93f},
  {89.4063f, 114.5f, 197.2f},
  {75.75f, 136.5f, 195.88f},
  {68.9688f, 146.5f, 202.12f},
  {95.6563f, 23.23f, 195.99f},
  {116.406f, 7.1f, 195.77f},
  {112.813f, 50.5f, 197.98f},
  {104.563f, 93.4f, 198.06f}
};

int main(int argc, char *argv[])
{
  int i;
  bool ospa = false;
  bool ospb = false;
  char* file_name_out = 0;
  char* file_name_header_out = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-ospa") == 0)
    {
      ospa = true;
    }
    else if (strcmp(argv[i],"-ospb") == 0)
    {
      ospb = true;
    }
    else if (strcmp(argv[i],"-o") == 0)
    {
      i++;
      file_name_out = argv[i];
    }
    else if (strcmp(argv[i],"-ho") == 0)
    {
      i++;
      file_name_header_out = argv[i];
    }
    else
    {
      usage();
    }
  }

  // open output header file 
  SMwriter* smwriter;
  FILE* file_header_out;

  if (file_name_header_out)
  {
    if (strstr(file_name_header_out, ".sma"))
    {
      file_header_out = fopen(file_name_header_out, "w");
    }
    else if (strstr(file_name_header_out, ".smb"))
    {
      file_header_out = fopen(file_name_header_out, "wb");
    }
    else
    {
      fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name_header_out);
      exit(0);
    }
    if (file_header_out == 0)
    {
      fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_header_out);
      exit(0);
    }
  }
  else
  {
    fprintf(stderr,"ERROR: no file name given for storing header of finalizemethod \n");
    exit(0);
  }

  if (strstr(file_name_header_out, ".sma"))
  {
    SMwriter_sma* smwriter_sma = new SMwriter_sma();
    smwriter_sma->open(file_header_out);
    smwriter = smwriter_sma;
  }
  else if (strstr(file_name_header_out, ".smb"))
  {
    SMwriter_smb* smwriter_smb = new SMwriter_smb();
    smwriter_smb->open(file_header_out);
    smwriter = smwriter_smb;
  }

  // open output point file 

  SPwriter* spwriter;
  FILE* file_out;

  if (file_name_out || ospa || ospb)
  {
    if (file_name_out)
    {
      if (strstr(file_name_out, ".spa") || ospa)
      {
        file_out = fopen(file_name_out, "w");
      }
      else if (strstr(file_name_out, ".spb") || ospb)
      {
        file_out = fopen(file_name_out, "wb");
      }
      else
      {
        if (file_name_out)
        {
          fprintf(stderr,"ERROR: cannot guess desired output from file name '%s'\n",file_name_out);
        }
        else
        {
          fprintf(stderr,"ERROR: no ouput format specified\n");
        }
        exit(0);
      }
      if (file_out == 0)
      {
        fprintf(stderr,"ERROR: cannot open '%s' for write\n", file_name_out);
        exit(0);
      }
    }
    else
    {
      file_out = stdout;
    }

    if ((file_name_out && strstr(file_name_out, ".spa")) || ospa)
    {
      SPwriter_spa* spwriter_spa = new SPwriter_spa();
      spwriter_spa->open(file_out);
      spwriter = spwriter_spa;
    }
    else if ((file_name_out && strstr(file_name_out, ".spb")) || ospb)
    {
      SPwriter_spb* spwriter_spb = new SPwriter_spb();
      spwriter_spb->open(file_out);
      spwriter = spwriter_spb;
    }
  }

  // simple example of outputting a streaming point cloud with SP_CLARKSON_2D finalization

  float bb_min_f[3];
  float bb_max_f[3];

  bb_min_f[0] = 0;
  bb_min_f[1] = 0;
  bb_min_f[2] = 192.25;

  bb_max_f[0] =  122.438f;
  bb_max_f[1] =  158.5f;
  bb_max_f[2] = 215.71f;

  // *first* we write the header in form of a separate mesh
  // assume we have 12 vertices and 20 triangles in the coarse
  // clarkson finalization mesh. the first vertex is the infinite
  // vertex ... by giving it a high negative z coordinate we can
  // nicely render it for verification purposes with sm_viewer

  // all triangles referening the first vertex are infinite triangles

  smwriter->set_nverts(clarkson2d_mesh_nverts);
  smwriter->set_nfaces(clarkson2d_mesh_nfaces);
  smwriter->set_boundingbox(bb_min_f, bb_max_f);

  for (i = 0; i < clarkson2d_mesh_nverts; i++)
  {
    smwriter->write_vertex(clarkson2d_mesh_verts[i]);
  }

  for (i = 0; i < clarkson2d_mesh_nfaces; i++)
  {
    smwriter->write_triangle(clarkson2d_mesh_faces[i]);
  }

  smwriter->close();
  fclose(file_header_out);
  delete smwriter;

  // *then* we write the finalized point stream

  spwriter->set_npoints(data_npoints);
  spwriter->set_boundingbox(bb_min_f, bb_max_f);
  spwriter->set_datatype(SP_FLOAT);
  spwriter->set_finalizemethod(SP_CLARKSON_2D);

  spwriter->write_header();

  // write all the points

  for (i = 0; i < data_npoints; i++)
  {
    spwriter->write_point(data_points[i]);
  }

  // finalize all the cells of the coarse clarkson finalization mesh 

  for (i = 0; i < clarkson2d_mesh_nfaces; i++)
  {
    spwriter->write_finalize_cell(i);
  }

  // usually the point and cell finalization should be intermixed, but
  // since i do not have computed "real" finalization in this example
  // i have to finalize all cells all at the end.

  spwriter->close();
  if (file_out != stdout) fclose(file_out);
  delete spwriter;

  return 1;
}
