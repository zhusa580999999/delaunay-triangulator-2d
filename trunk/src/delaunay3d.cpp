/*
===============================================================================

  FILE:  delaunay3d.cpp
  
  CONTENTS:
  
    a little module that takes 3d points as input and produces a streaming
    delaunay-tetrahedralized volume mesh as output

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
    yuanxin liuy@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 September 2005 -- joining leos and my code on the boardgame evening 
  
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

#include "spdelaunay3d.h"

#include "svwriter_nil.h"
#include "svwriter_sva.h"
#include "svwriter_svb.h"

#include "vec3fv.h"

#ifdef _WIN32
extern "C" FILE* fopenGzipped(const char* filename, const char* mode);
#endif

void usage()
{
	printf("usage: \n\n");
	printf("spdelaunay3d -i points.spa -o mesh.sva\n");
	printf("spdelaunay3d -i points.spb -o mesh.svb\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  bool ispa = false;
  bool ispb = false;
  bool osva = false;
  bool osvb = false;
  bool osvc = false;
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
    else if (strcmp(argv[i],"-osva") == 0)
    {
      osva = true;
    }
    else if (strcmp(argv[i],"-osvb") == 0)
    {
      osvb = true;
    }
    else if (strcmp(argv[i],"-osvc") == 0)
    {
      osvc = true;
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
  SVwriter* svwriter = 0;

  if ((file_name_out && strstr(file_name_out, ".sva")) || osva)
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
    SVwriter_sva* svwriter_sva = new SVwriter_sva(); 
	  svwriter_sva->open(outfile);
    svwriter = svwriter_sva;
  }
  else if ((file_name_out && strstr(file_name_out, ".svb")) || osvb)
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
    SVwriter_svb* svwriter_svb = new SVwriter_svb(); 
	  svwriter_svb->open(outfile);
    svwriter = svwriter_svb;
  }
  else
  {
    svwriter = new SVwriter_nil(); 
		fprintf(stderr,"producing no output.\n");
  }

  // create delaunay triangulator

  SPdelaunay3D* spdelaunay3d = new SPdelaunay3D();

  spdelaunay3d->open(svwriter);

	if ((spreader->bb_min_f != 0) && (spreader->bb_max_f != 0))
  {
    spdelaunay3d->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
	}
  else
  {
		fprintf(stderr,"WARNING: input point cloud has no bounding box.\n");
    exit(1);
  }

  // start clock
  timer_start = clock();

  SPevent event;
  while (event = spreader->read_event())
  {
    if (event == SP_POINT)
    {
      spdelaunay3d->write_point(spreader->p_pos_f);
    }
    else if (event == SP_FINALIZED_CELL)
    {
//      fprintf(stderr,"finalizing ... %d\n", spreader->final_idx);	
      if (spreader->final_idx >= 0) spdelaunay3d->write_finalize_cell(spreader->final_idx);
    }
    else
    {
      fprintf(stderr,"WARNING: unknown event %d at p_count\n",event, spreader->p_count);
    }
  }

  spdelaunay3d->close();

  // stop clock
  timer_finish = clock();
  fprintf(stderr,"needed %f seconds\n",  (double)(timer_finish-timer_start)/CLOCKS_PER_SEC);	

  // delete spreader
  fprintf(stderr,"read %d points\n", spreader->p_count);
  spreader->close();
  delete spreader;

  if (svwriter)
  {
    fprintf(stderr,"produced mesh with %d vertices and %d tetrahedra\n", svwriter->v_count, svwriter->c_count);
    svwriter->close();
    delete svwriter;
  }

  return 1;
}
