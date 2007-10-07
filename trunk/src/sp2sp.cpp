/*
===============================================================================

  FILE:  sp2sp.cpp
  
  CONTENTS:
  
    This (experimental) program converts streaming points to and from the
    specified format.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    11 February 2005 -- created despite this anhaltenden muskelkaters
  
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spreader_spa.h"
#include "spreader_spb.h"
#include "spreader_spc.h"
#include "spreader_node.h"
#include "spreader_ply.h"
#include "spreader_raw_d.h"
#include "spreader_raw.h"

#include "spwriter_spa.h"
#include "spwriter_spb.h"
#include "spwriter_spc.h"
#include "spwriter_raw_d.h"
#include "spwriter_raw.h"

#include "smreader_sma.h"
#include "smreader_smb.h"
#include "smreader_smc.h"

#include "vec3dv.h"
#include "vec3fv.h"

#ifdef _WIN32
extern "C" FILE* fopenGzipped(const char* filename, const char* mode);
extern "C" int gettime_in_msec();
extern "C" int gettime_in_sec();
extern "C" void settime();
#endif

static void usage()
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"sp2sp -i mesh.ply -o points.spb \n");
  fprintf(stderr,"sp2sp -i points.txt.gz -o points.spb\n");
  fprintf(stderr,"sp2sp -i points.raw -ho points.hdr\n");
  fprintf(stderr,"sp2sp -i mesh.smc -ho points.spb\n");
  fprintf(stderr,"sp2sp -i mesh.ply.gz -o points.raw\n");
  fprintf(stderr,"sp2sp -h\n");
  exit(1);
}

static bool saveHeader(FILE* file, int npoints, SPdatatype datatype, const void* bb_min, const void* bb_max)
{
  fwrite(&npoints,sizeof(int),1,file);
  fwrite(&datatype,sizeof(SPdatatype),1,file);
  if (datatype == SP_DOUBLE)
  {
    fwrite(bb_min,sizeof(double),3,file);
    fwrite(bb_max,sizeof(double),3,file);
  }
  else if (datatype == SP_FLOAT)
  {
    fwrite(bb_min,sizeof(float),3,file);
    fwrite(bb_max,sizeof(float),3,file);
  }
  else if (datatype == SP_INT)
  {
    fwrite(bb_min,sizeof(int),3,file);
    fwrite(bb_max,sizeof(int),3,file);
  }
  else
  {
    return false;
  }
  return true;
}

static bool loadHeader(FILE* file, SPreader* spreader)
{
  fread(&(spreader->npoints),sizeof(int),1,file);
  fread(&(spreader->datatype),sizeof(SPdatatype),1,file);
  if (spreader->datatype == SP_DOUBLE)
  {
    if (spreader->bb_min_d == 0) spreader->bb_min_d = new double[3];
    if (spreader->bb_max_d == 0) spreader->bb_max_d = new double[3];
    fread(spreader->bb_min_d,sizeof(double),3,file);
    fread(spreader->bb_max_d,sizeof(double),3,file);
  }
  else if (spreader->datatype == SP_FLOAT)
  {
    if (spreader->bb_min_f == 0) spreader->bb_min_f = new float[3];
    if (spreader->bb_max_f == 0) spreader->bb_max_f = new float[3];
    fread(spreader->bb_min_f,sizeof(float),3,file);
    fread(spreader->bb_max_f,sizeof(float),3,file);
  }
  else if (spreader->datatype == SP_INT)
  {
    if (spreader->bb_min_i == 0) spreader->bb_min_i = new int[3];
    if (spreader->bb_max_i == 0) spreader->bb_max_i = new int[3];
    fread(spreader->bb_min_i,sizeof(int),3,file);
    fread(spreader->bb_max_i,sizeof(int),3,file);
  }
  else
  {
    return false;
  }
  return true;
}

int main(int argc, char *argv[])
{
  int i;
  bool dry = false;
  bool isma = false;
  bool ismb = false;
  bool ismc = false;
  bool ispa = false;
  bool ispb = false;
  bool ispc = false;
  bool inode = false;
  bool ospa = false;
  bool ospb = false;
  bool ospc = false;
  bool clamp2 = false;
  bool clamp3 = false;
  int subsample = 0;
  float clamp_min_f[3];
  float clamp_max_f[3];
  bool pertube = false;
  bool unfinalize = false; 
  bool force_double = false;
  bool omit_bb = false;
  bool store_ = false;
  int bits = 16;
  char* file_name_in = 0;
  char* file_name_out = 0;
  char* file_name_header_out = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-dry") == 0)
    {
      dry = true;
    }
    else if (strcmp(argv[i],"-isma") == 0)
    {
      isma = true;
    }
    else if (strcmp(argv[i],"-ismb") == 0)
    {
      ismb = true;
    }
    else if (strcmp(argv[i],"-ismc") == 0)
    {
      ismc = true;
    }
    else if (strcmp(argv[i],"-ispa") == 0)
    {
      ispa = true;
    }
    else if (strcmp(argv[i],"-ispb") == 0)
    {
      ispb = true;
    }
    else if (strcmp(argv[i],"-ispc") == 0)
    {
      ispc = true;
    }
    else if (strcmp(argv[i],"-ospa") == 0)
    {
      ospa = true;
    }
    else if (strcmp(argv[i],"-ospb") == 0)
    {
      ospb = true;
    }
    else if (strcmp(argv[i],"-ospc") == 0)
    {
      ospc = true;
    }
    else if (strcmp(argv[i],"-pertube") == 0)
    {
      pertube = true;
    }
    else if (strcmp(argv[i],"-unfinalize") == 0)
    {
      unfinalize = true;
    }
    else if (strcmp(argv[i],"-double") == 0)
    {
      force_double = true;
    }
    else if (strcmp(argv[i],"-omit_bb") == 0)
    {
      omit_bb = true;
    }
    else if (strcmp(argv[i],"-clamp2") == 0)
    {
      if (i + 4 >= argc)
      {
        usage();
      }
      clamp2 = true;
      i++;
      clamp_min_f[0] = (float)atof(argv[i]);
      i++;
      clamp_min_f[1] = (float)atof(argv[i]);
      i++;
      clamp_max_f[0] = (float)atof(argv[i]);
      i++;
      clamp_max_f[1] = (float)atof(argv[i]);
      clamp_min_f[2] = 0;
      clamp_max_f[2] = 0;
    }
    else if (strcmp(argv[i],"-clamp3") == 0)
    {
      if (i + 6 >= argc)
      {
        usage();
      }
      clamp3 = true;
      i++;
      clamp_min_f[0] = (float)atof(argv[i]);
      i++;
      clamp_min_f[1] = (float)atof(argv[i]);
      i++;
      clamp_min_f[2] = (float)atof(argv[i]);
      i++;
      clamp_max_f[0] = (float)atof(argv[i]);
      i++;
      clamp_max_f[1] = (float)atof(argv[i]);
      i++;
      clamp_max_f[2] = (float)atof(argv[i]);
    }
    else if (strcmp(argv[i],"-subsample") == 0 || strcmp(argv[i],"-sub") == 0)
    {
      i++;
      subsample = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-b") == 0 || strcmp(argv[i],"-bits") == 0)
    {
      i++;
      bits = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-i") == 0)
    {
      i++;
      file_name_in = argv[i];
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

  SMreader* smreader = 0;
  SPreader* spreader = 0;
  FILE* file_in;
  
  if (file_name_in)
  {
    if (strstr(file_name_in, ".gz"))
    {
#ifdef _WIN32
      if (strstr(file_name_in, ".sma.gz") || strstr(file_name_in, ".spa.gz") || strstr(file_name_in, ".node.gz") || ispa || isma || inode)
      {
        file_in = fopenGzipped(file_name_in, "r");
      }
      else
      {
        file_in = fopenGzipped(file_name_in, "rb");
      }
#else
      fprintf(stderr,"ERROR: cannot open gzipped file '%s'\n",file_name_in);
      exit(1);
#endif
    }
    else
    {
      if (strstr(file_name_in, ".sma") || strstr(file_name_in, ".spa")  || strstr(file_name_in, ".node") || ispa || isma || inode)
      {
        file_in = fopen(file_name_in, "r");
      }
      else
      {
        file_in = fopen(file_name_in, "rb");
      }
    }
  }
  else
  {
    if (ispa || ispb || ispc || isma || ismb || ismc || inode)
    {
      file_in = stdin;
    }
    else
    {
      fprintf(stderr,"ERROR: no input file format specified\n");
      exit(0);
    }
  }

  if (file_in == 0)
  {
    fprintf(stderr,"ERROR: cannot open '%s' for read\n", file_name_in);
    exit(0);
  }
  
  if (file_name_in)
  {
    if (strstr(file_name_in, ".sma") || isma)
    {
      SMreader_sma* smreader_sma = new SMreader_sma();
      smreader_sma->open(file_in);

      if (smreader_sma->bb_min_f == 0 || smreader_sma->bb_max_f == 0)
      {
        SMevent event;

        fprintf(stderr,"need additional pass to compute bounding box\n");
        
        smreader_sma->bb_min_f = new float[3];
        smreader_sma->bb_max_f = new float[3];

        while ((event = smreader_sma->read_element()) > SM_EOF)
        {
          if (event == SM_VERTEX)
          {
            VecCopy3fv(smreader_sma->bb_min_f, smreader_sma->v_pos_f);
            VecCopy3fv(smreader_sma->bb_max_f, smreader_sma->v_pos_f);
            break;
          }
        }
        while ((event = smreader_sma->read_element()) > SM_EOF)
        {
          if (event == SM_VERTEX)
          {
            VecUpdateMinMax3fv(smreader_sma->bb_min_f, smreader_sma->bb_max_f, smreader_sma->v_pos_f);
          }
        }
        smreader_sma->close();
        fclose(file_in);
        fprintf(stderr, "bb_min_f %g %g %g\n", smreader_sma->bb_min_f[0], smreader_sma->bb_min_f[1], smreader_sma->bb_min_f[2]);
        fprintf(stderr, "bb_max_f %g %g %g\n", smreader_sma->bb_max_f[0], smreader_sma->bb_max_f[1], smreader_sma->bb_max_f[2]);
        if (strstr(file_name_in, ".gz"))
        {
#ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "r");
#endif
        }
        else
        {
          file_in = fopen(file_name_in, "r");
        }
        smreader_sma->open(file_in);
      }
      smreader = smreader_sma;
    }
    else if (strstr(file_name_in, ".spa") || ispa)
    {
      SPreader_spa* spreader_spa = new SPreader_spa();
      spreader_spa->open(file_in);

      if (spreader_spa->bb_min_f == 0 || spreader_spa->bb_max_f == 0)
      {
        SPevent event;

        fprintf(stderr,"need additional pass to compute bounding box\n");
        
        spreader_spa->bb_min_f = new float[3];
        spreader_spa->bb_max_f = new float[3];

        while ((event = spreader_spa->read_event()) > SP_EOF)
        {
          if (event == SP_POINT)
          {
            VecCopy3fv(spreader_spa->bb_min_f, spreader_spa->p_pos_f);
            VecCopy3fv(spreader_spa->bb_max_f, spreader_spa->p_pos_f);
            break;
          }
        }
        while ((event = spreader_spa->read_event()) > SP_EOF)
        {
          if (event == SP_POINT)
          {
            VecUpdateMinMax3fv(spreader_spa->bb_min_f, spreader_spa->bb_max_f, spreader_spa->p_pos_f);
          }
        }
        spreader_spa->close();
        fclose(file_in);
        fprintf(stderr, "bb_min_f %g %g %g\n", spreader_spa->bb_min_f[0], spreader_spa->bb_min_f[1], spreader_spa->bb_min_f[2]);
        fprintf(stderr, "bb_max_f %g %g %g\n", spreader_spa->bb_max_f[0], spreader_spa->bb_max_f[1], spreader_spa->bb_max_f[2]);
        if (strstr(file_name_in, ".gz"))
        {
#ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "r");
#endif
        }
        else
        {
          file_in = fopen(file_name_in, "r");
        }
        spreader_spa->open(file_in);
      }
      spreader = spreader_spa;
    }
    else if (strstr(file_name_in, ".smb") || ismb)
    {
      SMreader_smb* smreader_smb = new SMreader_smb();
      smreader_smb->open(file_in);

      if (smreader_smb->bb_min_f == 0 || smreader_smb->bb_max_f == 0)
      {
        SMevent event;

        fprintf(stderr,"need additional pass to compute bounding box\n");
        
        smreader_smb->bb_min_f = new float[3];
        smreader_smb->bb_max_f = new float[3];

        while ((event = smreader_smb->read_element()) > SM_EOF)
        {
          if (event == SM_VERTEX)
          {
            VecCopy3fv(smreader_smb->bb_min_f, smreader_smb->v_pos_f);
            VecCopy3fv(smreader_smb->bb_max_f, smreader_smb->v_pos_f);
            break;
          }
        }
        while ((event = smreader_smb->read_element()) > SM_EOF)
        {
          if (event == SM_VERTEX)
          {
            VecUpdateMinMax3fv(smreader_smb->bb_min_f, smreader_smb->bb_max_f, smreader_smb->v_pos_f);
          }
        }
        smreader_smb->close();
        fclose(file_in);
        fprintf(stderr, "bb_min_f %g %g %g\n", smreader_smb->bb_min_f[0], smreader_smb->bb_min_f[1], smreader_smb->bb_min_f[2]);
        fprintf(stderr, "bb_max_f %g %g %g\n", smreader_smb->bb_max_f[0], smreader_smb->bb_max_f[1], smreader_smb->bb_max_f[2]);
        if (strstr(file_name_in, ".gz"))
        {
#ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "rb");
#endif
        }
        else
        {
          file_in = fopen(file_name_in, "rb");
        }
        smreader_smb->open(file_in);
      }
      smreader = smreader_smb;
    }
    else if (strstr(file_name_in, ".spb") || ispb)
    {
      SPreader_spb* spreader_spb = new SPreader_spb();
      spreader_spb->open(file_in);

      if (spreader_spb->bb_min_f == 0 || spreader_spb->bb_max_f == 0)
      {
        SPevent event;

        fprintf(stderr,"need additional pass to compute bounding box\n");
        
        spreader_spb->bb_min_f = new float[3];
        spreader_spb->bb_max_f = new float[3];

        while ((event = spreader_spb->read_event()) > SP_EOF)
        {
          if (event == SP_POINT)
          {
            VecCopy3fv(spreader_spb->bb_min_f, spreader_spb->p_pos_f);
            VecCopy3fv(spreader_spb->bb_max_f, spreader_spb->p_pos_f);
            break;
          }
        }
        while ((event = spreader_spb->read_event()) > SP_EOF)
        {
          if (event == SP_POINT)
          {
            VecUpdateMinMax3fv(spreader_spb->bb_min_f, spreader_spb->bb_max_f, spreader_spb->p_pos_f);
          }
        }
        spreader_spb->close();
        fclose(file_in);
        fprintf(stderr, "bb_min_f %g %g %g\n", spreader_spb->bb_min_f[0], spreader_spb->bb_min_f[1], spreader_spb->bb_min_f[2]);
        fprintf(stderr, "bb_max_f %g %g %g\n", spreader_spb->bb_max_f[0], spreader_spb->bb_max_f[1], spreader_spb->bb_max_f[2]);
        if (strstr(file_name_in, ".gz"))
        {
#ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "rb");
#endif
        }
        else
        {
          file_in = fopen(file_name_in, "rb");
        }
        spreader_spb->open(file_in);
      }
      spreader = spreader_spb;
    }
    else if (strstr(file_name_in, ".node"))
    {
      SPevent event;
      SPreader_node* spreader_node = new SPreader_node();
      spreader_node->open(file_in);

      if (force_double)
      {      
        spreader_node->force_double();
      }

      if (!omit_bb)
      {
        fprintf(stderr,"need additional pass to compute bounding box\n");

        if (force_double)
        {      
          spreader_node->bb_min_d = new double[3];
          spreader_node->bb_max_d = new double[3];

          while ((event = spreader_node->read_event()) > SP_EOF)
          {
            if (event == SP_POINT)
            {
              VecCopy3dv(spreader_node->bb_min_d, spreader_node->p_pos_d);
              VecCopy3dv(spreader_node->bb_max_d, spreader_node->p_pos_d);
              break;
            }
          }
          while ((event = spreader_node->read_event()) > SP_EOF)
          {
            if (event == SP_POINT)
            {
              VecUpdateMinMax3dv(spreader_node->bb_min_d, spreader_node->bb_max_d, spreader_node->p_pos_d);
            }
          }
          fprintf(stderr, "bb_min_d %f %f %f\n", spreader_node->bb_min_d[0], spreader_node->bb_min_d[1], spreader_node->bb_min_d[2]);
          fprintf(stderr, "bb_max_d %f %f %f\n", spreader_node->bb_max_d[0], spreader_node->bb_max_d[1], spreader_node->bb_max_d[2]);
        }
        else
        {
          spreader_node->bb_min_f = new float[3];
          spreader_node->bb_max_f = new float[3];

          while ((event = spreader_node->read_event()) > SP_EOF)
          {
            if (event == SP_POINT)
            {
              VecCopy3fv(spreader_node->bb_min_f, spreader_node->p_pos_f);
              VecCopy3fv(spreader_node->bb_max_f, spreader_node->p_pos_f);
              break;
            }
          }
          while ((event = spreader_node->read_event()) > SP_EOF)
          {
            if (event == SP_POINT)
            {
              VecUpdateMinMax3fv(spreader_node->bb_min_f, spreader_node->bb_max_f, spreader_node->p_pos_f);
            }
          }
          fprintf(stderr, "bb_min_f %g %g %g\n", spreader_node->bb_min_f[0], spreader_node->bb_min_f[1], spreader_node->bb_min_f[2]);
          fprintf(stderr, "bb_max_f %g %g %g\n", spreader_node->bb_max_f[0], spreader_node->bb_max_f[1], spreader_node->bb_max_f[2]);
        }
        spreader_node->close();
        fclose(file_in);

        if (strstr(file_name_in, ".gz"))
        {
  #ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "r");
  #endif
        }
        else
        {
          file_in = fopen(file_name_in, "r");
        }
        spreader_node->open(file_in);
        if (force_double)
        {
          spreader_node->force_double();
        }
      }
      spreader = spreader_node;
    }
    else if (strstr(file_name_in, ".raw_d"))
    {
      SPreader_raw_d* spreader_raw_d = new SPreader_raw_d();
      spreader_raw_d->open(file_in);

      if (!omit_bb)
      {
        fprintf(stderr,"need additional pass to compute bounding box\n");
        spreader_raw_d->compute_bounding_box();

        fprintf(stderr, "bb_min_d %f %f %f\n", spreader_raw_d->bb_min_d[0], spreader_raw_d->bb_min_d[1], spreader_raw_d->bb_min_d[2]);
        fprintf(stderr, "bb_max_d %f %f %f\n", spreader_raw_d->bb_max_d[0], spreader_raw_d->bb_max_d[1], spreader_raw_d->bb_max_d[2]);

        spreader_raw_d->close();
        fclose(file_in);

        if (strstr(file_name_in, ".gz"))
        {
  #ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "rb");
  #endif
        }
        else
        {
          file_in = fopen(file_name_in, "rb");
        }
        spreader_raw_d->open(file_in);
      }
      spreader = spreader_raw_d;
    }
    else if (strstr(file_name_in, ".raw"))
    {
      SPreader_raw* spreader_raw = new SPreader_raw();
      spreader_raw->open(file_in);

      if (!omit_bb)
      {
        fprintf(stderr,"need additional pass to compute bounding box\n");
        spreader_raw->compute_bounding_box();

        fprintf(stderr, "bb_min_f[0] = %f; bb_min_f[1] = %f; bb_min_f[2] = %f;\n", spreader_raw->bb_min_f[0], spreader_raw->bb_min_f[1], spreader_raw->bb_min_f[2]);
        fprintf(stderr, "bb_max_f[0] = %f; bb_max_f[1] = %f; bb_max_f[2] = %f;\n", spreader_raw->bb_max_f[0], spreader_raw->bb_max_f[1], spreader_raw->bb_max_f[2]);

        spreader_raw->close();
        fclose(file_in);

        if (strstr(file_name_in, ".gz"))
        {
  #ifdef _WIN32
          file_in = fopenGzipped(file_name_in, "rb");
  #endif
        }
        else
        {
          file_in = fopen(file_name_in, "rb");
        }
        spreader_raw->open(file_in);
      }
      spreader = spreader_raw;
    }
    else
    {
      fprintf(stderr,"ERROR: cannot determine which reader to use for '%s'\n", file_name_in);
      exit(0);
    }
  }
  else
  {
    if (isma)
    {
      SMreader_sma* smreader_sma = new SMreader_sma();
      smreader_sma->open(file_in);
      smreader = smreader_sma;
    }
    else if (ismb)
    {
      SMreader_smb* smreader_smb = new SMreader_smb();
      smreader_smb->open(file_in);
      smreader = smreader_smb;
    }
    else
    {
      fprintf(stderr,"ERROR: cannot determine which reader to use\n");
      exit(0);
    }
  }

  SPwriter* spwriter;
  FILE* file_out;

  if (file_name_out || ospa || ospb || ospc)
  {
    if (dry)
    {
      fprintf(stderr,"dry write pass. no output.\n");
      file_out = 0;
    }
    else
    {
      if (file_name_out)
      {
        if (strstr(file_name_out, ".spa") || ospa)
        {
          file_out = fopen(file_name_out, "w");
        }
        else if (strstr(file_name_out, ".spb") || strstr(file_name_out, ".spc") || strstr(file_name_out, ".raw") || ospb || ospc)
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
    }

    if (file_name_out)
    {
      if (strstr(file_name_out, ".spa") || ospa)
      {
        if (file_out)
        {
          SPwriter_spa* spwriter_spa = new SPwriter_spa();
          spwriter_spa->open(file_out);
          spwriter = spwriter_spa;
        }
        else
        {
          spwriter = 0;
        }
      }
      else if (strstr(file_name_out, ".spb") || ospb)
      {
        SPwriter_spb* spwriter_spb = new SPwriter_spb();
        spwriter_spb->open(file_out);
        spwriter = spwriter_spb;
      }
      else if (strstr(file_name_out, ".raw"))
      {
        SPwriter_raw* spwriter_raw = new SPwriter_raw();
        spwriter_raw->open(file_out);
        spwriter = spwriter_raw;
      }
      else
      {
        fprintf(stderr,"ERROR: cannot determine which writer to use for '%s'\n", file_name_out);
        exit(0);
      }
    }
    else
    {
      if (ospa)
      {
        if (file_out)
        {
          SPwriter_spa* spwriter_spa = new SPwriter_spa();
          spwriter_spa->open(file_out);
          spwriter = spwriter_spa;
        }
        else
        {
          spwriter = 0;
        }
      }
      else if (ospb)
      {
        if (file_out)
        {
          SPwriter_spb* spwriter_spb = new SPwriter_spb();
          spwriter_spb->open(file_out);
          spwriter = spwriter_spb;
        }
        else
        {
          spwriter = 0;
        }
      }
      else
      {
        fprintf(stderr,"ERROR: cannot determine which writer to use\n");
        exit(0);
      }
    }
  }
  else
  {
    fprintf(stderr,"read pass only.\n");
    spwriter = 0;
  }

  int event;

#ifdef _WIN32
  settime();
#endif

  if (smreader)
  {
    if (spwriter)
    {
      if (smreader->nverts != -1) spwriter->set_npoints(smreader->nverts);

      if (pertube)
      {
        float pmax[3];
        float pvec[3];
        VecSubtract3fv(pmax, smreader->bb_max_f, smreader->bb_min_f);
        VecSelfScalarDiv3fv(pmax, (float)(1 << 18));
        
        if (smreader->bb_min_f || smreader->bb_max_f) fprintf(stderr,"WARNING: bounding box info not written (because of pertubation).\n");

        while (event = smreader->read_element())
        {
          switch (event)
          {
          case SM_VERTEX:
            pvec[0] = smreader->v_pos_f[0]+pmax[0]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
            pvec[1] = smreader->v_pos_f[1]+pmax[1]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
            pvec[2] = smreader->v_pos_f[2]+pmax[2]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
            spwriter->write_point(pvec);
            break;
          default:
            break;
          }
        }
      }
      else
      {
        if (smreader->bb_min_f || smreader->bb_max_f) spwriter->set_boundingbox(smreader->bb_min_f, smreader->bb_max_f);
  
        spwriter->write_header();

        while (event = smreader->read_element())
        {
          switch (event)
          {
          case SM_VERTEX:
            spwriter->write_point(smreader->v_pos_f);
            break;
          case SM_TRIANGLE:
          case SM_FINALIZED:
          default:
            break;
          }
        }
      }

      fprintf(stderr,"v_count %d p_count %d\n",smreader->v_count,spwriter->p_count);

      spwriter->close();
      if (file_out && file_name_out) fclose(file_out);
      delete spwriter;
    }
    else
    {
      float bb_min_f[3];
      float bb_max_f[3];

      while (event = smreader->read_element())
      {
        switch (event)
        {
        case SM_VERTEX:
          if (smreader->v_count == 1)
          {
            VecCopy3fv(bb_min_f, smreader->v_pos_f);
            VecCopy3fv(bb_max_f, smreader->v_pos_f);
          }
          else
          {
            VecUpdateMinMax3fv(bb_min_f, bb_max_f, smreader->v_pos_f);
          }
          break;
        case SM_TRIANGLE:
        case SM_FINALIZED:
        default:
          break;
        }
      };

      fprintf(stderr,"nverts %d\n",smreader->v_count);
      fprintf(stderr,"bb_min_f %f %f %f\n",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
      fprintf(stderr,"bb_max_f %f %f %f\n",bb_max_f[0],bb_max_f[1],bb_max_f[2]);
    }

#ifdef _WIN32
    fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif

    smreader->close();
    if (file_in && file_name_in) fclose(file_in);
    delete smreader;
  }
  else if (spreader)
  {
    if (spwriter)
    {
	  spwriter->set_finalizemethod(spreader->finalizemethod);
      if (spreader->npoints != -1) spwriter->set_npoints(spreader->npoints);

      if (pertube)
      {
        float pmax[3];
        float pvec[3];

        if (spreader->bb_min_f || spreader->bb_max_f) fprintf(stderr,"WARNING: bounding box info not written (because of pertubation).\n");
        fprintf(stderr,"WARNING: any existing finalization info will be omitted (because of pertubation).\n");

        if (clamp2 || clamp3)
        {
          VecSubtract3fv(pmax, clamp_max_f, clamp_min_f);
          VecSelfScalarDiv3fv(pmax, (float)(1 << 18));
          while (event = spreader->read_event())
          {
            switch (event)
            {
            case SP_POINT:
              VecClamp3fv(pvec, clamp_min_f, clamp_max_f, spreader->p_pos_f);
              pvec[0] = pvec[0]+pmax[0]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              pvec[1] = pvec[1]+pmax[1]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              pvec[2] = pvec[2]+pmax[2]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              spwriter->write_point(pvec);
              break;
            default:
              break;
            }
          }
        }
        else
        {
          VecSubtract3fv(pmax, spreader->bb_max_f, spreader->bb_min_f);
          VecSelfScalarDiv3fv(pmax, (float)(1 << 18));
          while (event = spreader->read_event())
          {
            switch (event)
            {
            case SP_POINT:
              pvec[0] = spreader->p_pos_f[0]+pmax[0]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              pvec[1] = spreader->p_pos_f[1]+pmax[1]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              pvec[2] = spreader->p_pos_f[2]+pmax[2]*(((float)(rand()-RAND_MAX))/((float)RAND_MAX));
              spwriter->write_point(pvec);
              break;
            default:
              break;
            }
          }
        }
      }
      else
      {

        if (spreader->datatype == SP_FLOAT)
        {
          if (clamp2 || clamp3)
          {
            spwriter->set_boundingbox(clamp_min_f, clamp_max_f);
  
            spwriter->write_header();

            fprintf(stderr,"writing clamped single-precision points ...\n");

            float pvec[3];
            while (event = spreader->read_event())
            {
              switch (event)
              {
              case SP_POINT:
                VecClamp3fv(pvec, clamp_min_f, clamp_max_f, spreader->p_pos_f);
                spwriter->write_point(pvec);
                break;
              default:
                break;
              }
            }
          }
          else if (subsample)
          {
            if (spreader->bb_min_f && spreader->bb_max_f) spwriter->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
  
            spwriter->write_header();

            fprintf(stderr,"writing subsampled %d single-precision points ...\n",subsample);

            while (event = spreader->read_event())
            {
              switch (event)
              {
              case SP_POINT:
                if (rand()%subsample == 0)
                {
                  spwriter->write_point(spreader->p_pos_f);
                }
                break;
              case SP_FINALIZED_CELL:
                if (!unfinalize) spwriter->write_finalize_cell(spreader->final_idx);
                break;
              default:
                break;
              }
            }
          }
          else
          {
            if (spreader->bb_min_f && spreader->bb_max_f) spwriter->set_boundingbox(spreader->bb_min_f, spreader->bb_max_f);
  
            spwriter->write_header();

            fprintf(stderr,"writing single-precision points ...\n");

            while (event = spreader->read_event())
            {
              switch (event)
              {
              case SP_POINT:
                spwriter->write_point(spreader->p_pos_f);
                break;
              case SP_FINALIZED_CELL:
                if (!unfinalize) spwriter->write_finalize_cell(spreader->final_idx);
                break;
              default:
                break;
              }
            }
          }
        }
        else if (spreader->datatype == SP_DOUBLE)
        {
          if (spreader->bb_min_d && spreader->bb_max_d) spwriter->set_boundingbox(spreader->bb_min_d, spreader->bb_max_d);
          spwriter->write_header();

          if (subsample)
          {  
            fprintf(stderr,"writing subsampled %d double-precision points ...\n",subsample);

            while (event = spreader->read_event())
            {
              switch (event)
              {
              case SP_POINT:
                if (rand()%subsample == 0)
                {
                  spwriter->write_point(spreader->p_pos_d);
                }
                break;
              case SP_FINALIZED_CELL:
                if (!unfinalize) spwriter->write_finalize_cell(spreader->final_idx);
                break;
              default:
                break;
              }
            }
          }
          else
          {
            fprintf(stderr,"writing double-precision points ...\n");
            while (event = spreader->read_event())
            {
              switch (event)
              {
              case SP_POINT:
                spwriter->write_point(spreader->p_pos_d);
                break;
              case SP_FINALIZED_CELL:
                if (!unfinalize) spwriter->write_finalize_cell(spreader->final_idx);
                break;
              default:
                break;
              }
            }
          }
        }
      }

      fprintf(stderr,"p_count %d p_count %d\n",spreader->p_count,spwriter->p_count);

      spwriter->close();
      if (file_out && file_name_out) fclose(file_out);
      delete spwriter;
    }
    else
    {
      if (spreader->datatype == SP_FLOAT)
      {
        int p_count = 0;
        float bb_min_f[3];
        float bb_max_f[3];

        while (event = spreader->read_event())
        {
          switch (event)
          {
          case SP_POINT:
            if (p_count == 0)
            {
              VecCopy3fv(bb_min_f, spreader->p_pos_f);
              VecCopy3fv(bb_max_f, spreader->p_pos_f);
            }
            else
            {
              VecUpdateMinMax3fv(bb_min_f, bb_max_f, spreader->p_pos_f);
            }
            p_count++;
            break;
          default:
            break;
          }
        };

        fprintf(stderr,"npoints %d\n",spreader->p_count);
        fprintf(stderr,"bb_min_f %g %g %g\n",bb_min_f[0],bb_min_f[1],bb_min_f[2]);
        fprintf(stderr,"bb_max_f %g %g %g\n",bb_max_f[0],bb_max_f[1],bb_max_f[2]);

        if (file_name_header_out)
        {
          FILE* file_header = fopen(file_name_header_out, "wb");
          saveHeader(file_header,p_count,SP_FLOAT,bb_min_f,bb_max_f);
          fclose(file_header);
        }      
      }
      else if (spreader->datatype == SP_DOUBLE)
      {
        int p_count = 0;
        double bb_min_d[3];
        double bb_max_d[3];

        while (event = spreader->read_event())
        {
          switch (event)
          {
          case SP_POINT:
            if (p_count == 0)
            {
              VecCopy3dv(bb_min_d, spreader->p_pos_d);
              VecCopy3dv(bb_max_d, spreader->p_pos_d);
            }
            else
            {
              VecUpdateMinMax3dv(bb_min_d, bb_max_d, spreader->p_pos_d);
            }
            p_count++;
            break;
          default:
            break;
          }
        };

        fprintf(stderr,"npoints %d\n",spreader->p_count);
        fprintf(stderr,"bb_min_d %f %f %f\n",bb_min_d[0],bb_min_d[1],bb_min_d[2]);
        fprintf(stderr,"bb_max_d %f %f %f\n",bb_max_d[0],bb_max_d[1],bb_max_d[2]);

        if (file_name_header_out)
        {
          FILE* file_header = fopen(file_name_header_out, "wb");
          saveHeader(file_header,p_count,SP_DOUBLE,bb_min_d,bb_max_d);
          fclose(file_header);
        }
      }
    }

    spreader->close();
    if (file_in && file_name_in) fclose(file_in);
    delete spreader;
  }

#ifdef _WIN32
  fprintf(stderr,"needed %6.3f seconds\n",0.001f*gettime_in_msec());
#endif

  return 1;
}
