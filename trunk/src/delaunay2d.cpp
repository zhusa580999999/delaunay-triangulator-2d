/*
===============================================================================

  FILE:  delaunay2d.cpp
  
  CONTENTS:
  
    a little module that takes 2d points as input and produces a streaming
    delaunay mesh as output

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
    yuanxin liuy@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    15 January 2006 -- also process double-precision input
    10 September 2005 -- added a '-scatter' command line option
    15 August 2005 -- added a binary point format for efficiency
    13 August 2005 -- supporting also Jonathan's unstreaming ASCII format
    09 August 2005 -- made the initial quadtree hack better and faster
    07 August 2005 -- adapted from Leo's code to support *real* streaming
    22 July 2005 -- last update by Yuanxin Liu
  
===============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "spreader_spa.h"
#include "spreader_spb.h"
#include "spreader_node.h"

#include "spreadscattered.h"

#include "spdelaunay2d.h"

#include "smwriter_sma.h"
#include "smwriter_smb.h"
#include "smwriter_nil.h"

#include "vec3fv.h"

#ifdef _WIN32
extern "C" FILE* fopenGzipped(const char* filename, const char* mode);
#endif

void usage()
{
	printf("usage: \n\n");
	printf("spdelaunay2d -i pointcloud.spa\n");
	printf("spdelaunay2d -ispb -osma > mesh.sma\n");
	printf("spdelaunay2d -i pointcloud.spb -o mesh.smb\n");
	printf("spdelaunay2d -i pointcloud.spa -o mesh.sma\n");
	printf("spdelaunay2d -ispb -o mesh.smb < pointcloud.spb\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  bool ispa = false;
  bool ispb = false;
  bool osma = false;
  bool osmb = false;
  bool osmc = false;
	char* file_name_in = 0;
	char* file_name_out = 0;
  int scatter = 0;
  clock_t timer_start, timer_finish;

  if (argc == 1)
  {
    usage();
  }

  for (int i = 1; i < argc; i++)
  {
    if ( strcmp("-i", argv[i]) == 0)
    {
      i++;
			file_name_in = argv[i];
		}
    else if ( strcmp("-o", argv[i]) == 0)
    {
      i++;
			file_name_out = argv[i];
		}
    else if (strcmp(argv[i],"-ispa") == 0)
    {
      ispa = true;
    }
    else if (strcmp(argv[i],"-ispb") == 0)
    {
      ispb = true;
    }
    else if (strcmp(argv[i],"-osma") == 0)
    {
      osma = true;
    }
    else if (strcmp(argv[i],"-osmb") == 0)
    {
      osmb = true;
    }
    else if (strcmp(argv[i],"-osmc") == 0)
    {
      osmc = true;
    }
    else if (strcmp(argv[i],"-scatter") == 0)
    {
      i++;
      scatter = atoi(argv[i]);;
    }
    else
    {
      usage();
    }
	}

  // open input file
  FILE* infile = 0;
  SPreader* spreader = 0;

  if ((file_name_in && strstr(file_name_in, ".spa")) || ispa)
  {
    if (file_name_in)
    {
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        infile = fopenGzipped(file_name_in,"r");
#else
		    fprintf(stderr,"ERROR: cannot open gzipped file '%s' for read.\n", file_name_in);
#endif
      }
      else
      {
        infile = fopen(file_name_in,"r");
      }
	    if (!infile)
      {
		    fprintf(stderr,"ERROR: cannot open file '%s' for read.\n", file_name_in);
		    return 0;
	    }
    }
    else
    {
      infile = stdin;
    }
    SPreader_spa* spreader_spa = new SPreader_spa();
    spreader_spa->open(infile);
    spreader = spreader_spa;
  }
  else if (file_name_in && strstr(file_name_in, ".node"))
  {
    if (strstr(file_name_in, ".gz"))
    {
#ifdef _WIN32
      infile = fopenGzipped(file_name_in,"r");
#else
		  fprintf(stderr,"ERROR: cannot open gzipped file '%s' for read.\n", file_name_in);
#endif
    }
    else
    {
      infile = fopen(file_name_in,"r");
    }
	  if (!infile)
    {
		  fprintf(stderr,"ERROR: cannot open file '%s' for read.\n", file_name_in);
		  return 0;
	  }
    SPreader_node* spreader_node = new SPreader_node();
    spreader_node->open(infile);
    // compute bounding box
    spreader_node->bb_min_f = new float[3];
    spreader_node->bb_max_f = new float[3];
    // read first point
    spreader_node->read_event();
    VecCopy3fv(spreader_node->bb_min_f, spreader_node->p_pos_f);
    VecCopy3fv(spreader_node->bb_max_f, spreader_node->p_pos_f);
    // read remaining points
    while(spreader_node->read_event())
    {
      VecUpdateMinMax3fv(spreader_node->bb_min_f, spreader_node->bb_max_f, spreader_node->p_pos_f);
    }
    spreader_node->close();
    fclose(infile);
    if (strstr(file_name_in, ".gz"))
    {
#ifdef _WIN32
      infile = fopenGzipped(file_name_in,"r");
#endif
    }
    else
    {
      infile = fopen(file_name_in,"r");
    }
	  if (!infile)
    {
		  fprintf(stderr,"ERROR: cannot open file '%s' for second read.\n", file_name_in);
		  return 0;
	  }
    spreader_node->open(infile);
    spreader = spreader_node;
  }
  else if ((file_name_in && strstr(file_name_in, ".spb")) || ispb)
  {
    if (file_name_in)
    {
      if (strstr(file_name_in, ".gz"))
      {
#ifdef _WIN32
        infile = fopenGzipped(file_name_in,"rb");
#else
		    fprintf(stderr,"ERROR: cannot open gzipped file '%s' for read.\n", file_name_in);
#endif
      }
      else
      {
        infile = fopen(file_name_in,"rb");
      }
	    if (!infile)
      {
		    fprintf(stderr,"ERROR: cannot open file '%s' for read.\n", file_name_in);
		    return 0;
	    }
    }
    else
    {
      infile = stdin;
    }
    SPreader_spb* spreader_spb = new SPreader_spb();
    spreader_spb->open(infile);
    spreader = spreader_spb;
  }

  // maybe read scattered

  if (scatter)
  {
    SPreadScattered* spreadscattered = new SPreadScattered();
    spreadscattered->open(spreader, scatter);
    spreader = spreadscattered;
  }

  // open output file

  FILE* outfile = 0;
  SMwriter* smwriter = 0;

  if ((file_name_out && strstr(file_name_out, ".sma")) || osma)
  {
    if (file_name_out)
    {
      outfile = fopen(file_name_out,"w");
	    if (!outfile)
      {
		    fprintf(stderr,"ERROR: cannot open file '%s' for write.\n", file_name_out);
		    return 0;
	    }
    }
    else
    {
      outfile = stdout;
    }
    SMwriter_sma* smwriter_sma = new SMwriter_sma(); 
	  smwriter_sma->open(outfile);
    smwriter = smwriter_sma;
  }
  else if ((file_name_out && strstr(file_name_out, ".smb")) || osmb)
  {
    if (file_name_out)
    {
      outfile = fopen(file_name_out,"wb");
	    if (!outfile)
      {
		    fprintf(stderr,"ERROR: cannot open file '%s' for write.\n", file_name_out);
		    return 0;
	    }
    }
    else
    {
      outfile = stdout;
    }
    SMwriter_smb* smwriter_smb = new SMwriter_smb(); 
	  smwriter_smb->open(outfile);
    smwriter = smwriter_smb;
  }
  else
  {
		fprintf(stderr,"producing no output.\n");
    smwriter = new SMwriter_nil(); 
  }

  // create delaunay triangulator

  SPdelaunay2D* spdelaunay2d = new SPdelaunay2D();

  spdelaunay2d->open(smwriter);

  // start clock
  timer_start = clock();

  SPevent event;

  if (spreader->datatype == SP_FLOAT)
  {
    if (spreader->bb_min_f && spreader->bb_max_f)
    {
      spdelaunay2d->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
    }
    else
    {
  		fprintf(stderr,"ERROR: single-precicion input point cloud has no bounding box.\n");
      exit(0);
    }
    while (event = spreader->read_event())
    {
      if (event == SP_POINT)
      {
        spdelaunay2d->write_point(spreader->p_pos_f);
      }
      else if (event == SP_FINALIZED_CELL)
      {
        spdelaunay2d->write_finalize_cell(spreader->final_idx);
      }
      else
      {
        fprintf(stderr,"WARNING: unknown event %d at p_count\n",event, spreader->p_count);
      }
    }
  }
  else // should be SP_DOUBLE
  {
    if (spreader->bb_min_d && spreader->bb_max_d)
    {
      spdelaunay2d->set_boundingbox(spreader->bb_min_d, spreader->bb_max_d);
    }
    else
    {
  		fprintf(stderr,"ERROR: double-precicion input point cloud has no bounding box.\n");
      exit(0);
    }
    while (event = spreader->read_event())
    {
      if (event == SP_POINT)
      {
        spdelaunay2d->write_point(spreader->p_pos_d);
      }
      else if (event == SP_FINALIZED_CELL)
      {
        spdelaunay2d->write_finalize_cell(spreader->final_idx);
      }
      else
      {
        fprintf(stderr,"WARNING: unknown event %d at p_count\n",event, spreader->p_count);
      }
    }
  }

  spdelaunay2d->close();

  // stop clock
  timer_finish = clock();
  fprintf(stderr,"needed %f seconds\n",  (double)(timer_finish-timer_start)/CLOCKS_PER_SEC);	

  // delete spreader
  fprintf(stderr,"read %d points\n", spreader->p_count);
  spreader->close();
  delete spreader;

  if (smwriter)
  {
    fprintf(stderr,"produced mesh with %d vertices and %d triangles\n", smwriter->v_count, smwriter->f_count);
    smwriter->close();
    delete smwriter;
  }

  return 1;
}

