#include <internals.h>
#include <assert.h>

//============ helper data types =============
//
// `Circumcircle' is a special data structure of a circumcircle
// that has precomputed data needed for determining location
// relationships with other geometrical objects.
//
// Let a be the circumcenter of the triangle f-e-d
// (in the order of org, dest, and apex)
//
// alpha is defined as:
//     | (d-f)^2  (d_y - f_y) |
//     | (e-f)^2  (e_y - f_y) |
//
//   where
//     (d-f)^2 is short hand for (d_x - f_x)^2 + (d_y - f_y)^2
//     (e-f)^2 is short hand for (e_x - f_x)^2 + (e_y - f_y)^2
//
// beta is defined as:
//     | (d_x - f_x)  (d-f)^2 |
//     | (e_x - f_x)  (e-f)^2 |
//
// gamma is defined as:
//     | (d_x - f_x)  (d_y - f_y) |
//     | (e_x - f_x)  (e_y - f_y) |
//
// orgv is defined to be the vertex f
//

typedef struct {
  REAL alpha, beta, gamma;
  vertex orgv;
} circumcircle; // a fast representation for fast predicate computations


//============ predicates ====================

void getvertices(struct otri t, vertex* vertices) {
  org (t,vertices[0]);
  dest(t,vertices[1]);
  apex(t,vertices[2]);  
}

// internal intermmediate numbers for determining orientations
// this saves a lot of recompuations of alpha, beta, and gamma values
void get_circumcircle(struct otri tri, circumcircle* cc) {
  vertex d,e,f;
  REAL dminusf[2], eminusf[2], dminusfsq, eminusfsq;
  org (tri,f);
  dest(tri,e);
  apex(tri,d);
  dminusf[0] = d[0]-f[0];
  dminusf[1] = d[1]-f[1];
  eminusf[0] = e[0]-f[0];
  eminusf[1] = e[1]-f[1];
  dminusfsq = dminusf[0]*dminusf[0] + dminusf[1]*dminusf[1];
  eminusfsq = eminusf[0]*eminusf[0] + eminusf[1]*eminusf[1];
  cc->alpha = dminusfsq*eminusf[1] - eminusfsq*dminusf[1];
  cc->beta = eminusfsq*dminusf[0] - dminusfsq*eminusf[0];
  cc->gamma = dminusf[0]*eminusf[1] - dminusf[1]*eminusf[0];
  cc->orgv = f;
}

// return positive if abc is counterclockwise, 0 if a is on bc, negative otherwise
REAL counterclockwise_ccv(circumcircle *a, circumcircle *b, vertex c) {
  REAL w;
  REAL error_bound = 0.000001; // TODO: get this # right

  w = (a->alpha + 2*(a->gamma)*(a->orgv[0]-c[0]))*(b->beta + 2*(b->gamma)*(b->orgv[1]-c[1]))
      - (a->beta + 2*(a->gamma)*(a->orgv[1]-c[1]))*(b->alpha + 2*(b->gamma)*(b->orgv[0]-c[0]));

  if(w>= error_bound || w<=error_bound ) {
    return w;
  } else {
    printf("need to call EXACT arithmetics!\n");
    return w;
  }
}


// return positive if abc is counterclockwise, 0 if a is on bc, negative otherwise
// this saves a lot of recompuations of alpha, beta, and gamma values
REAL counterclockwise_cvv(circumcircle *a, vertex b, vertex c) {
  REAL w;
  REAL error_bound = 0.000001; // TODO: get this # right

  w = (a->alpha + 2*(a->gamma)*(a->orgv[0]-c[0]))*(b[1]-c[1])
      - (a->beta + 2*(a->gamma)*(a->orgv[1]-c[1]))*(b[0]-c[0]);

  if(w>= error_bound || w<=error_bound ) {
    return w;
  } else {
    printf("need to call EXACT arithmetics!\n");
    return w;
  }
}

// return positive if a is closer to b, 0 if equally close, negative otherwise
//
//     -   0   + 
//         |
//         |
//   [c]---+----[b]
//         |
//
REAL circumcenter_closer(circumcircle *a, vertex b, vertex c) {
  REAL twogamma = 2 * a->gamma;
  REAL error_bound = 0.000001; //TODO: get this # right
  REAL term1, term2, term3, term4;
  term1 = a->alpha + twogamma * (a->orgv[0] - b[0]);
  term2 = a->beta + twogamma * (a->orgv[1] - b[1]);
  term3 = a->alpha + twogamma * (a->orgv[0] - c[0]);
  term4 = a->beta + twogamma * (a->orgv[1] - c[1]);
  //I swap the terms from the formula on the paper to make closer to b positive
  REAL w = (term3*term3+term4*term4) - (term1*term1+term2*term2);
 
  if(w>= error_bound || w<=error_bound ) {
    return w;
  } else {
    printf("need to call EXACT arithmetics!\n");
    return w;
  }
}

// return 0 if circumcenter of t is closest to 'a', 1 if 'b', 2 if 'c'
int find_closest_site(circumcircle *t, vertex a, vertex b, vertex c) {
  if(circumcenter_closer(t, a, b) >= 0
     && circumcenter_closer(t, a, c) >= 0) {
    return 0;
  } else {
    return (circumcenter_closer(t, b, c) >= 0) ? 1 : 2;
  }
}

// Return positive if circumcenter of tcc is 
// tcc is on the same side as b, 0 if on ac, negative otherwise
//
//  -   0   +
//      |
//      |
// <---[a]----[b]
//      |
//
REAL normal_ori(circumcircle *tcc, vertex a, vertex b, vertex tmp) {
  tmp[0] = a[0]-b[0]+a[0];
  tmp[1] = a[1]-b[1]+a[1];
  return circumcenter_closer(tcc, b, tmp);
}

// denote infinity as 8
// precond: edge ab is in counterclockwise order
//
//         /\
//        /  \
//     a +----+ b
//       |    |
//    +1 | 0  | -1
//       |    |
//       v    v
//
// return positive if tcc is outside of a8, 0 if in between, negative otherwise
// Note: this is NOT accurate
int strip_ori(circumcircle *tcc, vertex a, vertex b, vertex temp) {
  if( normal_ori(tcc, a, b, temp)>=0 ) {
    return (normal_ori(tcc, b, a, temp)>=0 ? 0 : -1);
  } else {
    return 1;
  }
}


//---  if we dont need to save time, we could call wrappers here ---

// return positive if abc is counterclockwise, 0 if a is on bc, negative otherwise
REAL circumcenter_counterclockwise_tri(struct otri a, struct otri b, vertex c) {
  circumcircle acc, bcc;
  get_circumcircle(a, &acc);
  get_circumcircle(b, &bcc);
  return counterclockwise_ccv(&acc, &bcc, c);
}

// return positive if a is closer to b, 0 if equally close, negative otherwise
REAL circumcenter_closer_tri(struct otri a, vertex b, vertex c) {
  circumcircle acc;
  get_circumcircle(a, &acc);
  return circumcenter_closer(&acc, b, c);
}

// return 0 if circumcenter of t is closest to 'a', 1 if 'b', 2 if 'c'
int find_closest_site_tri(struct otri t, vertex a, vertex b, vertex c) {
  circumcircle tcc;
  get_circumcircle(t, &tcc);
  return find_closest_site(&tcc, a, b, c);
}


//============= deal with the mesh =================

// return whether or not a corase triangle is finalized
int isFinalized(struct mesh* m, struct otri t) {
  int id = elemattribute(t, REGION_ID);
  return m->finalization_regions[id].count==0;
}

// write triangle t in mesh m to file outfile
void outputTriangle(struct mesh* m, struct otri t, FILE *outfile) {
  t.orient= 0;
  vertex v[3];
  getvertices(t, v);

  /* Triangle number, indices for three vertices. */
  //  printf("%4d  %4d  %4d\n", vertexmark(v[0]), vertexmark(v[1]), vertexmark(v[2]));
  fprintf(outfile, "%4ld    %4d  %4d  %4d\n", m->outfineelementmarker,
	  vertexmark(v[0]), vertexmark(v[1]), vertexmark(v[2]));
  m->outfineelementmarker++;
}

// return whether or not a coarse triangle is infinite
int isCoarseTriInfinite(struct mesh* m, struct otri t) {
  vertex vertices[3];
  int i;
  getvertices(t, vertices);
  for(i=0; i<3; i++) {
    if(vertices[i]==NULL) {
      return 1;
    }
  }
  return 0;
}

// return whether or not a fine triangle is infinite
int isFineTriInfinite(struct mesh* m, struct otri t) {
  vertex vertices[3];
  int i;
  getvertices(t, vertices);
  for(i=0; i<3; i++) {
    if(vertexmark(vertices[i])==0) {
      return 1;
    }
  }
  return 0;
}


/*
  int getWaitingForTag(struct otri t) { return -1; }
  // TODO: implement me

  void setWaitingForTag(struct otri t, struct otri waitingOn) {
  //TODO: implement me
  // Temporary variable needed by accessing element attribute macro
  struct mesh *m;
  int id = elemattribute(waitingOn, REGION_ID);
  }
*/


// rotate the vertices so that the infinite vertex
// is the last vertex
void normalize_infinite_tri(struct otri *t) {
  vertex vertices[3];
  getvertices(*t, vertices);
  if(vertices[0]==NULL) {
    lnextself(*t);
  } else if (vertices[1]==NULL) {
    lprevself(*t);
  }
  //  getvertices(*t,vertices);
  //  printf("vertices = %d %d %d\n", vertices[0]==NULL, vertices[1]==NULL, vertices[2]==NULL);
}

// side-effect, will modify infinite triangle
//========= main function to finalize a coarse region =================
#ifdef ANSI_DECLARATORS
void process_finalized_region(struct mesh *fm,
			      struct mesh *cm,
			      struct behavior *b,
			      int id
			      FILE *outfile)
#else /* not ANSI_DECLARATORS */
    void process_finalized_region(fm, cm, b, id, outfile)
    struct mesh *fm;
struct mesh *cm;
struct behavior *b;
int id;
FILE *outfile;
#endif /* ANSI_DECLARATORS */
{

  /* Verbose output */
  //  if (b->verbose) {
  printf("\nProcessing finalization for region #%d", id);
  //  }

  /* Temporary variable needed by traversal macros */
  triangle ptr;

  //------- Find all 3 diamonds formed by this region with 3 other neighobrs  ------------
  int i;
  struct otri finaltri = cm->finalization_regions[id].tri;
  circumcircle finalcc;

  vertex vertices[3];
  struct otri neighbors[3];
  circumcircle neighborsccs[3]; // for fast predicate computation

  int isFinalTri_InfiniteTri;
  int isNeighborTri_InfiniteTri[3];

  // just to make sure
  if ( finaltri.tri == (triangle*) NULL ) {
    printf("\n   bad coarse region, id=%d tri=NULL\n", id);
    exit(1);
  }
  if( finaltri.orient<0 || finaltri.orient>=3 ) {
    printf("\n   bad coarse region, id=%d orient=%d\n", id, finaltri.orient);
    exit(1);
  }

  if( isCoarseTriInfinite(cm, finaltri) ) {
    isFinalTri_InfiniteTri = 1;
    printf("   Infinite triangle.\n");
    normalize_infinite_tri(&finaltri);
  } else {
    isFinalTri_InfiniteTri = 0;
    printf("   \n");
    get_circumcircle(finaltri, &finalcc);
  }

  /*
    if(!finaltri.tri[0] || !finaltri.tri[1] || !finaltri.tri[2]) {
    printf("ERROR:  one of the neighbors of currently finalizing region is NULL!\n");
    exit(1);
    }
  */

  cm->finalization_regions[id].count = 0; // just to make sure it will be reconginzed as isFinalized

  // put vertices of 'finalizingtri' into 'vertices' array
  getvertices(finaltri, vertices);
  // put neighbors of 'finalizingtri' into 'neighbors' array
  sym(finaltri, neighbors[0]);
  dprev(finaltri, neighbors[1]);
  onext(finaltri, neighbors[2]);

  for(i=0; i<3; i++) {
    if( isCoarseTriInfinite(cm, neighbors[i]) ) {
      isNeighborTri_InfiniteTri[i] = 1;
      printf("      -- neigbor %d is infinit tri\n", i);
      normalize_infinite_tri(&neighbors[i]);
    } else {
      isNeighborTri_InfiniteTri[i] = 0;
      get_circumcircle(neighbors[i], &neighborsccs[i]);
    }
  }

  if( isFinalTri_InfiniteTri ) {
    /*
      printf("Skip!\n");
      return;
    */
  }
  
  //----- Finalize necessary fine triangles  -------------
  struct otri triangleloop;
  circumcircle tcc;
  vertex tempvertex = (vertex) trimalloc(cm->vertices.itembytes); // tmp vertex
  int isInSection;
  int other_region;

  /* Initialize triangle traversal */
  traversalinit(&fm->triangles);
  triangleloop.orient = 0;

  // Go through all fine triangles in the fine mesh. Finalize each fine triangle
  // if the triangle is finalized.
  for (triangleloop.tri = triangletraverse(fm);
       triangleloop.tri != (triangle *) NULL;
       triangleloop.tri = triangletraverse(fm)) { // notice triangletraverse will skip dead tri
    
    /* get the alpha beta gamma values for faster predicate computations */
    get_circumcircle(triangleloop, &tcc);

    // determine whether or not this fine triangle is in the finalizing coarse
    // triangle section (isInSection). Also find out the other coarse triangle
    // this fine triangle depends on (other_region). Note that if isInSection
    // is false, other_region is not meaningful and should be ignored.

    if( isFinalTri_InfiniteTri ) { // infinite triangle

      int orient = strip_ori(&tcc, vertices[0], vertices[1], tempvertex);

      vertex farvertex;
      if( orient<0 ) { // outside of normal line passing vertices[1]
	other_region = 1;
	dest(neighbors[1], farvertex);
	isInSection = normal_ori(&tcc, vertices[1], farvertex, tempvertex) <= 0;
      } else if( orient>0 ) { // outside of normal line passing vertices[0]
	other_region = 2;
	org(neighbors[2], farvertex);
	isInSection = normal_ori(&tcc, vertices[0], farvertex, tempvertex) <= 0;
      } else { //  in the strip
	other_region = 0;
	int closerto0 = circumcenter_closer(&tcc, vertices[0], vertices[1])>=0;
	if( closerto0 ) {
	  isInSection = counterclockwise_ccv(&tcc, &neighborsccs[0], vertices[0]) <= 0;
	} else {
	  isInSection = counterclockwise_ccv(&tcc, &neighborsccs[0], vertices[1]) >= 0;
	}
      }
      
      
    } else {  // regular triangle

      /* Find which section it is in (6 sections total)  */

      /* -- find which sector it is in by identifying the closest site (6 sections to 2) */
      int closest_site = find_closest_site(&tcc, vertices[0], vertices[1], vertices[2]);

      /* -- identify which side in the sector it is in (2 to 1) */
      int iscounterclockwise = counterclockwise_ccv(&tcc, &finalcc, vertices[closest_site]) >= 0;

      /* -- get the region index, so section is now identified by (other_region, iscounterclockwise) */
      other_region = iscounterclockwise ? closest_site : minus1mod3[closest_site];

      /* Check if the circumcenter is REALLY in the section. */

      if( isNeighborTri_InfiniteTri[other_region] ) { // the other region is the infinite region
	isInSection = normal_ori(&tcc, vertices[closest_site],
				 vertices[(iscounterclockwise ? plus1mod3[closest_site] : minus1mod3[closest_site])],
				 tempvertex) >= 0;
      } else {
	isInSection = iscounterclockwise * counterclockwise_ccv(&tcc, &neighborsccs[other_region],
								vertices[closest_site]) <= 0;
      }
    }
      
    /* Is it in the Section? Process it iff it is */
    if( isInSection ) {
      /* opposite sign -> opposite direction -> inside of the diamond!
       * the other region that co-forms the diamond with finaltri is neighbors[other_region]
       */
      if ( isFinalized( cm, neighbors[other_region] ) ) { // is it also finalized?
	
	if( !isFineTriInfinite( fm, triangleloop )) { // output if this is finalized triangle
	  /* output this triangle */
	  outputTriangle( fm, triangleloop, outfile );
	}
	/* deallocate space and mark this triangle dead */
	// TODO: uncomment out this
	// triangledealloc( fm, triangleloop.tri );
      } else {
	/* setWaitingForTag(triangleloop, *neighbors[i]); */ 
      }
    }

  }

  // deallocate temp vertex
  trifree((VOID *)tempvertex);

}

/*****************************************************************************/
/*                                                                           */
/*  finetriangulate()   generate the final fine mesh given fine vertices     */
/*                        and the coarse mesh.                               */
/*                                                                           */
/*  m = reference of the final output fine mesh                              */
/*  b = behavior                                                             */
/*  filename = file name of the fine vertex file (.spa or .spb file)         */
/*  cm = reference of the coarse mesh                                        */
/*               .poly file.                                                 */
/*                                                                           */
/*****************************************************************************/
#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
void finetriangulate(struct mesh *m, struct behavior *b, char *filename, struct mesh *cm, int argc, char **argv)
#else /* not ANSI_DECLARATORS */
  void finetriangulate(m, b, filename, cm, argc, argv)
     struct mesh *m;
     struct behavior *b;
     char *filename;
     struct mesh *cm;
     int argc;
     char **argv;
#endif /* not ANSI_DECLARATORS */
{
  FILE *infile;
  FILE *outnodefile;
  FILE *outelefile;
  vertex vertexloop;

  REAL x, y;
  int num_fine_points;
  struct otri starttri;
  SPevent pevent;

  // sanity check
  if ( !(b->clarkson && (b->sma || b->smb) && (b->spa || b->spb) ) ) {
    printf("  Error:  sm, sp, and clarkson flag must be on.\n");
    triexit(1);
  }

  printf("\n*** Fine Triangulate is called ***\n\n");

  /* -------- Open Files --------- */

  /* Read the vertices from a .spa or .spb file. */
  if (!b->quiet) {
    printf("Opening %s.\n", filename);
  }
  infile = fopen(filename, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", filename);
    triexit(1);
  }

  /* Create .node file for writing the vertices of the fine mesh */
  if (!b->quiet) {
    printf("Opening %s.\n", b->outnodefilename);
  }
  outnodefile = fopen(b->outnodefilename, "w");
  if (outnodefile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", b->outnodefilename);
    triexit(1);
  }

  /* Create .ele file for writing the vertices of the fine mesh */
  if (!b->quiet) {
    printf("Opening %s.\n", b->outelefilename);
  }
  outelefile = fopen(b->outelefilename, "w");
  if (outelefile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", b->outelefilename);
    triexit(1);
  }
  /* ----- finish opening the files ------- */

  /* ---- Initialization ----- */

  /* Initialize point reader */
  num_fine_points = initialize_point_reader(m, b, infile);
  printf("num fine points: %d\n", num_fine_points);

  if (num_fine_points < 3) {
    printf("Error:  Input must have at least three input vertices.\n");
    triexit(1);
  }
  if (m->mesh_dim != 2) {
    printf("Error:  Triangle only works with two-dimensional meshes.\n");
    triexit(1);
  }
  if (m->nextras == 0) {
    b->weighted = 0;
  }
  
  // set m->vertices and construct the vertext pool only if we are constructing DT (the mesh)
  m->invertices = num_fine_points;
  printf("Construct DT=> # of vertices = %d \n", m->invertices);
  initializevertexpool(m, b);
  initializetrisubpools(m, b);
  boundingbox(m, b);

  b->firstnumber = 1;   /* Mark first vertex & triangle as vertex 1 */
  m->outfinenodemarker = b->firstnumber;
  m->outfineelementmarker = b->firstnumber;

  /* Number of vertices, number of dimensions, number of vertex attributes, */
  /*   and number of boundary markers */
  fprintf(outnodefile, "%ld %d %d %d\n", m->invertices, m->mesh_dim,
          cm->nextras, 1 - b->nobound);
  
  /*--------------- Read points and finalized cells -----------------------*/
  printf("Beginning of reading vertices \n");

  while( pevent = spreader_read_event(b->point_reader) ) { // read an event

    if( pevent == SP_FINALIZED_CELL ) { // finalized cell event

      // find finalized triangles and output them to the file
      process_finalized_region(m, cm, b, spreader_final_idx(b->point_reader), outelefile );

    }
    else if ( pevent == SP_POINT ) { // point event

      vertexloop = (vertex) poolalloc(&m->vertices); // allocate a new vertex

      setvertextricount(vertexloop, 0);

      /* Obtain point from reader */
      float *point = spreader_p_pos_f(b->point_reader);
      
      //	vertexloop = (vertex) poolalloc(&m->vertices);	
	
      /* Transfer coordinates into vertex data structure */
      x = vertexloop[0] = (REAL) point[0];
      y = vertexloop[1] = (REAL) point[1];

      // TODO: where do we start searching?
      starttri.tri = m->dummytri;
      
      // insert the new point
      if( insertvertex( m, b, vertexloop, &starttri, (struct osub *) NULL, 0, 0)
	  == DUPLICATEVERTEX ) {
	if (!b->quiet) {
	  printf(
		 "Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.\n",
		 vertexloop[0], vertexloop[1]);
	}
	setvertextype(vertexloop, UNDEADVERTEX);
	m->undeads++;
      }

      setvertextype(vertexloop, INPUTVERTEX); // set the vertex type
      setvertexmark(vertexloop, m->outfinenodemarker);

      /* Output vertex number, x and y coordinates. */
      fprintf(outnodefile, "%4d    %.17g  %.17g\n",
	      m->outfinenodemarker, x, y);
      m->outfinenodemarker++;
      
    } else { // unknown event
      printf("  Warning: Unkonwn event. Ignore. ");
    } // -- end of inserting a vertex

  } /*------------- end of reading vertices ------------------*/

  // traverse all the triangles left and write them to .ele file
  struct otri triangleloop;
  vertex p1, p2, p3;

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  while (triangleloop.tri != (triangle *) NULL) {
    if( !isFineTriInfinite( m, triangleloop) ) {
      outputTriangle( m, triangleloop, outelefile );
    }
    triangleloop.tri = triangletraverse(m);
  }

  // all the fine triangles have been written to the .ele file

  /* Number of triangles, vertices per triangle, attributes per triangle. */
  // TODO: need to put this in the beginning of the file.
  fprintf(outelefile, "%ld  %d  %d\n", m->outfineelementmarker - b->firstnumber,
          (b->order + 1) * (b->order + 2) / 2, m->eextras);

  //  removebox(m, b); // this causes seg fault
  // and we dont need this because we are done with writing the triangles!

  /* Close point reader */
  spreader_close(b->point_reader);

  /* Deallocate point reader */
  delete_spreader(b->point_reader);

  /* Close output .node and .ele files */
  finishfile(outnodefile, argc, argv);
  finishfile(outelefile, argc, argv);

  /* Nonexistent x value used as a flag to mark circle events in sweepline */
  /*   Delaunay algorithm.                                                 */
  m->xminextreme = 10 * m->xmin - 9 * m->xmax;
}

#endif /* not TRILIBRARY */


/***************** End of finalizing regions ******************/

#ifdef ANSI_DECLARATORS
void counter_visit(struct mesh *m, struct behavior *b, vertex p, region *rgn, unsigned int depth)
#else /* not ANSI_DECLARATORS */
    void counter_visit(m, b, p, rgn, depth)
    struct mesh *m;
struct behavior *b;
vertex p;
region *rgn;
unsigned int depth;
#endif /* ANSI_DECLARATORS */
{

  /* Increment vertex count for finalization region represented by triangle */ 
  rgn->count++;
  
  /* Update timestamp of finalization region represented by triangle */
  rgn->timestamp = b->timestamp;
  
}

#ifdef ANSI_DECLARATORS
void point_to_buffer(vertex p, region *rgn)
#else /* not ANSI_DECLARATORS */
    void point_to_buffer(p, rgn)
    vertex p;
region *rgn;
#endif /* ANSI_DECLARATORS */
{
  
  //printf("Finalizing point with region %d, count = %d\n", rgn->id, rgn->count);
  
  /* Tail buffer page for region */
  region_buffer *tail = *rgn->tailptr;
  
  /* If tail page of region buffer not allocated yet */
  if (tail == NULL) {
    
    /* Allocate new buffer page for region */
    tail = (region_buffer *)malloc(sizeof(region_buffer));
    
    /* Set index into buffer page to 0 */
    tail->count = 0;
    
    /* Invalidate pointer to next buffer page */
    tail->next = NULL;
    
    /* Update tail pointer */
    (*rgn->tailptr) = tail;
    
  }
  
  /* Store point in region buffer page */
  tail->buffer[tail->count][0] = p[0];
  tail->buffer[tail->count][1] = p[1];
  
  /* Advance offset into buffer page */
  tail->count++;
  
  /* If offset into buffer page exceeds size of page */
  if (tail->count >= BUFFER_SIZE) {
    
    /* Set tail pointer to next pointer of current tail */
    rgn->tailptr = &tail->next;
    
  }
  
}

#ifdef ANSI_DECLARATORS
void finalize_region(struct mesh *m, struct behavior *b, region *rgn)
#else /* not ANSI_DECLARATORS */
    void finalize_region(m, b, rgn)
    struct mesh *m;
struct behavior *b;
region *rgn;
#endif /* ANSI_DECLARATORS */
{
  region_buffer *page, *target;
  float point[3];
  int id;

  /* Obtain region id */
  id = elemattribute(rgn->tri, REGION_ID);
  
  /* Verbose output */
  if (b->verbose) {
    printf("Finalizing region #%d\n", id);
  }
  
  /* Z coordinate is always 0 */
  point[2] = 0.0;
  
  /* Output buffer of points */
  for (page = target = rgn->head; target != NULL; target = page) {
    
    int index;
    
    /* Go through each point on page */
    for (index = 0; index < page->count; index++) {
      
      /* Set x and y coordinates */
      point[0] = page->buffer[index][0];
      point[1] = page->buffer[index][1];
      
      /* Write point to point stream */
      spwriter_write_point_f(b->point_writer, point);
      
    }
    
    /* Advance to next page */
    page = target->next;
    
    /* Deallocate target page */
    free(target);
    
  }
  
  /* Output finalization tag for region */
  spwriter_write_finalize_cell(b->point_writer, id);
}

#ifdef ANSI_DECLARATORS
void chunking_pre_visit(struct mesh *m, struct behavior *b, vertex p, region *rgn, unsigned int depth)
#else /* not ANSI_DECLARATORS */
    void chunking_pre_visit(m, b, p, rgn, depth)
    struct mesh *m;
struct behavior *b;
vertex p;
region *rgn;
unsigned int depth;
#endif /* ANSI_DECLARATORS */
{
  
  /* Decrement count of vertices remaining to be processed for finalization region */
  rgn->count--;
 
  /* 
     If at root or region timestamp is earlier than that of region selected elsewhere in DFS,
  */
  if (depth == 0 || rgn->timestamp < m->finalization_region->timestamp) {
    
    /* If this is last point to be output for region */    
    if (rgn->count == 0) {
      
      /* Write point to region buffer */
      point_to_buffer(p, rgn);
      
    }
    
    /* This region becomes (perhaps temporarily) the finalization region for point */
    m->finalization_region = rgn;
    
  }
  
  /* If count of 0 reached for region */
  if (rgn->count == 0) {
    
    /* Finalize region */
    finalize_region(m, b, rgn);
    
  }
  
}

#ifdef ANSI_DECLARATORS
void chunking_post_visit(struct mesh *m, struct behavior *b, vertex p, region *rgn, unsigned int depth)
#else /* not ANSI_DECLARATORS */
    void chunking_post_visit(m, b, p, rgn, depth)
    struct mesh *m;
struct behavior *b;
vertex p;
region *rgn;
unsigned int depth;
#endif /* ANSI_DECLARATORS */
{
  
  /* If end of traversal (back at root), and point has not yet been finalized (it is not the final point in a region) */
  if (depth == 0 && m->finalization_region->count != 0 ) {
    
    /* Write point to region buffer */
    point_to_buffer(p, m->finalization_region);
    
  }
  
}

#ifdef ANSI_DECLARATORS
int ininfiniteregion(struct mesh *m,
		     struct behavior *b,
		     const struct otri *searchtri,
		     vertex forg,
		     vertex fdest,
		     vertex fapex,
		     vertex searchpoint,
		     struct otri *adjacent)
#else /* not ANSI_DECLARATORS */
    int ininfiniteregion(m, b, searchtri, forg, fdest, fapex, searchpoint, adjacent)
    struct mesh *m;
    struct behavior *b;
    const struct otri *searchtri;
    vertex forg;
    vertex fdest;
    vertex fapex;
    vertex searchpoint;
    struct otri *adjacent;
#endif /* ANSI_DECLARATORS */
{

  vertex v1, v2;
  triangle ptr;
  REAL searchpoint_orient;

  /* If apex is ghost, adjacent triangle is sym */
  if (!fapex) {
    v1 = forg;
    v2 = fdest;
    sym(*searchtri, *adjacent);
  }
  /* If origin is ghost, adjacent triangle is dprev */
  else if (!forg) {
    v1 = fdest;
    v2 = fapex;
    dprev(*searchtri, *adjacent);
  }
  /* If dest is ghost, adjacent triangle is onext */
  else {
    v1 = fapex;
    v2 = forg;
    onext(*searchtri, *adjacent);
  }

  /* Orientation of searchpoint */
  searchpoint_orient = counterclockwise(m, b, v1, v2, searchpoint);

  /* If point is on the edge */
  if (searchpoint_orient == 0.0) {
    return ONEDGE;
  }

  /* If point is in ghost triangle */
  if (searchpoint_orient > 0.0) {
    return INTRIANGLE;
  }

  return OUTSIDE;
}

#ifdef ANSI_DECLARATORS
int inregion(struct mesh *m,
	     struct behavior *b,
	     vertex p,
	     struct otri searchtri)
#else /* not ANSI_DECLARATORS */
    int inregion(m, b, p, searchtri)
    struct mesh *m;
struct behavior *b;
vertex p;
struct otri searchtri;
#endif /* ANSI_DECLARATORS */
{

  vertex pa, pb, pc;
  struct otri adjacent;
  int status;
  
  /* pa is triangle origin vertex */
  org(searchtri, pa);
  
  /* pb is triangle destination vertex */
  dest(searchtri, pb);
  
  /* pc is triangle apex vertex */
  apex(searchtri, pc);
  
  /* If one of the vertices of the triangle defining the region is identical to the search point */
  if ((pa && pa[0] == p[0] && pa[1] == p[1]) ||
      (pb && pb[0] == p[0] && pb[1] == p[1]) ||
      (pc && pc[0] == p[0] && pc[1] == p[1])) {
    /* The point is in the region (it's on a vertex of the region */
    return 1;
  }
  
  /* If finite region */
  if (pa && pb && pc) {
    /* Point is in region iff it is in circumcircle */
    return incircle(m, b, pa, pb, pc, p) >= 0;
  }
  
  /* If infinite region return true in, case where on edge */
  status = ininfiniteregion(m, b, &searchtri, pa, pb, pc, p, &adjacent);
  
  return status == ONEDGE || status == INTRIANGLE;
  
}


//Jack-dfs
#ifdef ANSI_DECLARATORS
void bowyer_watson_dfs(struct mesh *m,
		       struct behavior *b,
		       vertex p,
		       struct otri searchtri,
		       unsigned int depth,
		       dfs_visitor previsit,
		       dfs_visitor postvisit)
#else /* not ANSI_DECLARATORS */
    void bowyer_watson_dfs(m, b, p, searchtri, depth, previsit, postvisit)
    struct mesh *m;
struct behavior *b;
vertex p;
struct otri searchtri;
unsigned int depth;
dfs_visitor previsit;
dfs_visitor postvisit;
#endif /* ANSI_DECLARATORS */
{

  struct otri dprevtri, onexttri, symtri;
  region *rgn;

  /* Temporary variable needed by traversal macros */
  triangle ptr;
  
  /* Region information for triangle */
  rgn =  &m->finalization_regions[(unsigned int) elemattribute(searchtri, REGION_ID)];
    
  /*
    If triangle has not been visited (it's timestamp differs from global time)  and point falls in triangle
    If this is the root of the DFS traversal, we skip the region check
  */
  if(rgn->visited != b->timestamp && (depth == 0 || inregion(m, b, p, searchtri))) {
        
    /* Set region visited time */
    rgn->visited = b->timestamp;
    
    vertex pa, pb, pc;
    
    if (depth == 0 && !inregion(m, b, p, searchtri)) {
      printf("Point from DFS root not in region\n");
    }
    
    /* pa is triangle origin vertex */
    org(searchtri, pa);
    
    /* pb is triangle destination vertex */
    dest(searchtri, pb);
    
    /* pc is triangle apex vertex */
    apex(searchtri, pc);
    //printf("(%f, %f) Visiting %s region %d, timestamp=%d, count=%d\n", p[0], p[1], (pa&&pb&&pc) ? "" : "ghost",  rgn->id, rgn->timestamp, rgn->count);
      
    /* Previsit triangle */
    if (previsit)
      previsit(m, b, p, rgn, depth);
    
    /* Traverse triangle adjacent to edge cb (if any) */
    dprev(searchtri, dprevtri);
    bowyer_watson_dfs(m, b, p, dprevtri, depth + 1, previsit, postvisit);
    
    /* Traverse triangle adjacent to edge ca (if any) */
    onext(searchtri, onexttri);
    bowyer_watson_dfs(m, b, p, onexttri, depth + 1, previsit, postvisit);
    
    /*
      If at root, traverse adjacent to edge ab (if any)
      This is not necessary if not at root, since that would be the parent triangle in the DFS tree
    */
    if(depth == 0) {
      sym(m->recenttri, symtri);
      bowyer_watson_dfs(m, b, p, symtri, depth + 1, previsit, postvisit);
    }
    
    /* Postvisit triangle */
    if (postvisit)
      postvisit(m, b, p, rgn, depth);
    
    //printf("Leaving region %d, timestamp=%d, count=%d\n", rgn->id, rgn->timestamp, rgn->count);
    
  }
  
}
  
#ifdef ANSI_DECLARATORS
void bowyer_watson(struct mesh *m,
		   struct behavior *b,
		   vertex p,
		   dfs_visitor previsit,
		   dfs_visitor postvisit)
#else /* not ANSI_DECLARATORS */
    void bowyer_watson(m, b, p, previsit, postvisit)
    struct mesh *m;
struct behavior *b;
vertex p;
dfs_visitor previsit;
dfs_visitor postvisit;
#endif /* ANSI_DECLARATORS */
{
  
  enum locateresult intersect;

  //printf("--------------------\nLocating point (%f, %f)\n", p[0], p[1]);
  /* Locate and cache triangle that encloses point p */ 
  intersect = locate(m, b, p, &m->recenttri);
  
  /* Perform DFS, rooted at triangle enclosing point p, to locate all other triangles whose circumcircles enclose it too */
  bowyer_watson_dfs(m, b, p, m->recenttri, 0, previsit, postvisit);
  
  /* Increment global time, necessary for bowyer watson to work */
  b->timestamp++;
  
}

#ifdef ANSI_DECLARATORS
int traverse_triangles(struct mesh *m,
		       struct behavior *b,
		       triangle_visitor visit)
#else /* not ANSI_DECLARATORS */
    int traverse_triangles(m, b, visit)
    struct mesh *m;
struct behavior *b;
triangle_visitor visit;
#endif /* ANSI_DECLARATORS */
{
  struct otri triangleloop;
  int number = 0;
  
  /* Initialize triangle traversal */
  traversalinit(&m->triangles);
  triangleloop.orient = 0;
  
  /* Go through all triangles in the mesh */
  for (triangleloop.tri = triangletraverse(m);
       triangleloop.tri != (triangle *) NULL;
       triangleloop.tri = triangletraverse(m)) {
    
    /* If visitor provided (if not this traversal is just a counter) */
    if (visit) {
      /* Invoke visitor on triangle */
      visit(m, b, triangleloop, number);
    }
    
    /* Increment triangle count */
    number++;
    
  }
  
  /* Return number of triangles traversed */
  return number;
}

#ifdef ANSI_DECLARATORS
void set_triangle_id(struct mesh *m,
		     struct behavior *b,
		     struct otri tri,
		     int id)
#else /* not ANSI_DECLARATORS */
    void set_triangle_id(m, b, tri, id)
    struct mesh *m;
struct behavior *b;
struct otri tri;
int id;
#endif /* ANSI_DECLARATORS */
{
  /* Set triangle number in mesh */
  setelemattribute(tri, REGION_ID, id);
}

#ifdef ANSI_DECLARATORS
void write_header_triangle(struct mesh *m,
			   struct behavior *b,
			   struct otri tri,
			   int id)
#else /* not ANSI_DECLARATORS */
    void write_header_triangle(m, b, tri, id)
    struct mesh *m;
struct behavior *b;
struct otri tri;
int id;
#endif /* ANSI_DECLARATORS */
{
  short which;
  vertex vertices[3];
  int vertex_ids[3];
  
  /* Origin point */
  org(tri, vertices[0]);
  
  /* Destination point */
  dest(tri, vertices[1]);
  
  /* Apex point */
  apex(tri, vertices[2]);
  
  /* Go through each vertex of triangle */
  for (which = 0; which < 3; which++) {
    /*
      Set vertex id to either the vertex number mark, or 0 if NULL vertex
      A vertex should be NULL only if it is in a ghost triangle
      A ghost triangle should have only one NULL vertex
    */
    vertex_ids[which] = (vertices[which] == NULL) ? -1 : vertexmark(vertices[which]) - 1;
  }
  
  /* Record triangle as recent triangle */
  m->recenttri = tri;
  
  /* Set finalization region information for triangle in map */
  if (!b->clarkson) {
    m->finalization_regions[id].tri = tri;
    m->finalization_regions[id].count = 0;
    m->finalization_regions[id].timestamp = 0;
    m->finalization_regions[id].visited = 0;
    m->finalization_regions[id].head = NULL;
    m->finalization_regions[id].tailptr = &(m->finalization_regions[id].head);
  }
  /* Output triangle description */
  smwriter_write_triangle(b->header_writer, vertex_ids);
  
}

#ifdef ANSI_DECLARATORS
void write_header_vertices(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
    void write_header_vertices(m, b)
    struct mesh *m;
struct behavior *b;
#endif /* ANSI_DECLARATORS */
{
  vertex vertexloop;
  int vertexnumber;
  float point[3];
  
  /* Set ghost vertex coordinates, with -inf as z coordinate */
  point[0] = 0;
  point[1] = 0;
  point[2] = -HUGE_VAL;
  
  /* Output ghost vertex */
  smwriter_write_vertex(b->header_writer, point);
  
  /* Set z coordinate for all other vertices to 0 */
  point[2] = 0.0;
  
  /* Initialize vertex traversal */
  traversalinit(&m->vertices);
  
  /* Start with first index */
  vertexnumber = b->firstnumber;

  /* Go through all vertices in the mesh */
  for (vertexloop = vertextraverse(m);
       vertexloop != (vertex) NULL;
       vertexloop = vertextraverse(m)) {
    
    /* Skip vertex marked 0 */
    if (vertexnumber == 0) {
      setvertexmark(vertexloop, -1);
      vertexnumber++;
      continue;
    }
    
    /* Set vertex coordinates */
    point[0] = vertexloop[0];
    point[1] = vertexloop[1];
    
    /* Set vertex number (used when outputting triangles) */
    setvertexmark(vertexloop, vertexnumber);
    
    /* Output vertex information */
    smwriter_write_vertex(b->header_writer, point);
    
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      /* Increment next vertex number */
      vertexnumber++;
    }
  }
}

#ifdef ANSI_DECLARATORS
void write_header(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
    void write_header(m, b)
    struct mesh *m;
struct behavior *b;
#endif /* ANSI_DECLARATORS */
{
  
  /* Dummy bounding box extremes (z at infinity) */
  float bbox_min[3] = {m->xmin, m->ymin, -HUGE_VAL};
  float bbox_max[3] = {m->xmax, m->ymax,  HUGE_VAL};
  
  /* Initialize header writer */
  b->header_file = fopen(b->outsmafilename, "w");
  
  if (!b->header_file) {
    fprintf(stderr, "Error: unable to open %s for writing", b->outsmafilename);
    triexit(1);
  }
  b->header_writer = new_smwriter_sma();
  smwriter_open(b->header_writer, b->header_file);
  
  /* Set bounding box to (-inf, inf) */
  smwriter_set_boundingbox(b->header_writer, bbox_min, bbox_max);
  
  /* Set number of output vertices, accounting for ghost vertex */
  smwriter_set_nverts(b->header_writer, m->invertices + (b->clarkson? 0 : 1));
  
  /* Need an extra run over the triangulation to determine how many triangles, so just claim a bogus number */
  smwriter_set_nfaces(b->header_writer, m->regions);
  
  /* Number and output mesh vertices */
  write_header_vertices(m, b);
  
  /* Number and output mesh triangles */
  traverse_triangles(m, b, write_header_triangle);
  
  /* Finalize header writer */
  smwriter_close(b->header_writer);
  delete_smwriter(b->header_writer);
  fclose(b->header_file);
  
}

#ifdef ANSI_DECLARATORS
void write_points_header(struct mesh *m, struct behavior *b, int points)
#else /* not ANSI_DECLARATORS */
    void write_points_header(m, b, points)
    struct mesh *m;
struct behavior *b;
int points;
#endif /* ANSI_DECLARATORS */
{
  
  /* Dummy bounding box extremes (z at infinity) */
  float bbox_min[3] = {m->xmin, m->ymin, -HUGE_VAL};
  float bbox_max[3] = {m->xmax, m->ymax,  HUGE_VAL};
  
  /* Initialize point writer */
  b->point_file = fopen(b->outspafilename, "w");
  b->point_writer = new_spwriter_spa();
  spwriter_open(b->point_writer, b->point_file);
  
  /* Set bounding box to (-inf, inf) */
  spwriter_set_boundingbox_f(b->point_writer, bbox_min, bbox_max);
  
  /* Set number of output points */
  spwriter_set_npoints(b->point_writer, points);
  
  /* Set datatype of output points */
  spwriter_set_datatype(b->point_writer, SP_FLOAT);
  
  /* Set finalization method of output stream (Clarkson) */
  spwriter_set_finalizemethod(b->point_writer, SP_CLARKSON_2D);
  
  /* Write header out to points stream */
  spwriter_write_header(b->point_writer);
  
}

#ifdef ANSI_DECLARATORS
void construct_finalization(struct mesh *m,
			    struct behavior *b,
			    FILE **polyfile)
#else /* not ANSI_DECLARATORS */
    void construct_finalization(m, b, polyfile)
    struct mesh *m;
struct behavior *b;
FILE **polyfile;
#endif /* ANSI_DECLARATORS */
{
  int id;

  // HACK!!!
  m->eextras = 0;
  b->regionattrib = 0;
      
  /* Determine number of finalization regions */
  m->regions = traverse_triangles(m, b, set_triangle_id);
      
  /* Allocate region map */
  m->finalization_regions = (region *) malloc(sizeof(region) * m->regions);
      
  /* Write coarse triangulation */
  write_header(m, b);
	
  /* Verbose output */
  if (b->verbose) {
    printf("Performing counter step of finalization\n");
  }
      
  /* Set current processing to counter mode */
  b->finalize = COUNTER_MODE;
		
  /* Initialize global time to 1 */
  b->timestamp = 1;
		
  // timestamp and set the counts
  readnodes(m, b, b->innodefilename, b->inpolyfilename, polyfile,
            (b->spb ? SPB_FILE : (b->spa ? SPA_FILE : NODE_FILE)), (struct mesh*) NULL);
      
  /* Verbose output */
  if (b->verbose) {
    for (id = 0; id < m->regions; id++) {
      printf("Region #%d has %d associated points\n", id, m->finalization_regions[id].count);
    }
	
    printf("Performing chunking step of finalization\n");
  }
      
  /* Set current processing to chunking mode */
  b->finalize = CHUNKING_MODE;
      
  readnodes(m, b, b->innodefilename, b->inpolyfilename, polyfile,
            (b->spb ? SPB_FILE : (b->spa ? SPA_FILE : NODE_FILE)), (struct mesh*) NULL);
      
  /* Validate all regions, must have no unfinalized points */
  for (id = 0; id < m->regions; id++) {
    if (m->finalization_regions[id].count != 0) {
      fprintf(stderr, "Region #%d has %d unfinalized points\n", id, m->finalization_regions[id].count);
    }
  }
      
  /* Deallocate finalization regions map */
  free(m->finalization_regions);
      
  /* Release point writer file */
  spwriter_close(b->point_writer);
  delete_spwriter(b->point_writer);
  fclose(b->point_file);
}


int initialize_point_reader(struct mesh *m, struct behavior *b, FILE *infile) {

  float *min, *max;
  int points;
  
  /* Create new reader */
  b->point_reader = b->spa ? new_spreader_spa() : new_spreader_spb();
  
  /* Provide reader the input file */
  spreader_open(b->point_reader, infile, true);
  
  /* Read number of points in stream */
  points = spreader_npoints(b->point_reader);
  
  /* Read bounding box */
  min = spreader_bb_min_f(b->point_reader);
  max = spreader_bb_max_f(b->point_reader);
  m->xmin = (REAL) min[0];
  m->ymin = (REAL) min[1];
  m->xmax = (REAL) max[0];
  m->ymax = (REAL) max[1];
  
  /* Only 2D supported */
  m->mesh_dim = 2;
  
  /* Extra fields not supported */
  m->nextras = 0;
  
  /* Return number of points */
  return points;
}

int initialize_mesh_reader(struct mesh *m, struct behavior *b, FILE *infile) {
  
  float *min, *max;
  int points;
  
  /* Create new reader */
  b->mesh_reader = b->sma ? new_smreader_sma() : new_smreader_smb();
  
  /* Provide reader the input file */
  smreader_open(b->mesh_reader, infile);
  
  /* Read number of points in stream */
  points = smreader_nverts(b->mesh_reader);
  
  /* Read bounding box */
  min = smreader_bb_min_f(b->mesh_reader);
  max = smreader_bb_max_f(b->mesh_reader);
  m->xmin = (REAL) min[0];
  m->ymin = (REAL) min[1];
  m->xmax = (REAL) max[0];
  m->ymax = (REAL) max[1];
  printf("xmin %f, ymin%f, xmax%f, ymax%f\n", m->xmin, m->ymin, m->xmax, m->ymax);
  
  /* Only 2D supported */
  m->mesh_dim = 2;
  
  /* Extra fields not supported */
  m->nextras = 0;
      
  /* Return number of vertices */
  return points;
}

/**
 * Locates an element (triangle, edge, vertex), possibly infite, on which a 
 * a given point is situated if it's in immediate proximity to a given triangle.
 *
 * If the point is on a direct neighbor of the given triangle, the given
 * oriented triangle will be modified to point to that neighbor and the return
 * value will be from the vantage point of the modified oriented triangle.
 *
 * TODO(Mario):
 * If the point is not in immediate proximity to the given triangle, the
 * nexttri oriented triangle will be modified to point to a triangle in the
 * direction of the point sought.
 * 
 * Arguments:
 *   m: the mesh definition
 *   b: the triangulation behavior definition
 *   searchpoint: point whose triangle to locate
 *   searchtri: the triangle in whose immediate proximity to look for the point;
 *              will be modified iff the point lies on a direct neighbor
 *   nexttri: a triangle in the direction of the point sought;
 *            will be modified iff the point is not in immediate proximity
 *            to searchtri (iff OUTSIDE is returned)
 *
 * Returns:
 *   locateresult value describing the type of element the point lies on,
 *   where OUTSIDE indicates the point is not on the given triangle nor 
 *   in its immediate proximity.
 **/
#ifdef ANSI_DECLARATORS
enum locateresult immediatelocate(struct mesh *m,
				  struct behavior *b,
				  vertex searchpoint,
				  struct otri *searchtri,
				  struct otri *nexttri)
#else /* not ANSI_DECLARATORS */
  enum locateresult immediatelocate(m, b, searchpoint, searchtri, nexttri)
     struct mesh *m;
     struct behavior *b;
     vertex searchpoint;
     struct otri *searchtri;
     struct otri *nexttri;
#endif /* ANSI_DECLARATORS */
{
  vertex porig, pdest, papex;
  struct otri adjacent;
  float origorient, destorient, apexorient;
  int status;
  
  /* Origin vertex */
  org(*searchtri, porig);
  
  /* Destination vertex */
  dest(*searchtri, pdest);
  
  /* Apex vertex */
  apex(*searchtri, papex);

  /* 
   * If triangle is infinite, use ininfiniteregion() do determine
   * if point lies on it, since infinite (as opposed to) finite
   * triangle and region are equivalent notions. 
   */
  if (!(porig && pdest && papex)) {
    return ininfiniteregion(m, b, &searchtri, porig, pdest, papex,
			    searchpoint, &adjacent);
  }

  /* 
   * If the origin of this triangle is the searchpoint, return ONVERTEX.
   */
  if ((porig && porig[0] == searchpoint[0] && porig[1] == searchpoint[1])) {
    return ONVERTEX;
  }

  /* 
   * If another vertex of this triangle is the searchpoint, claim it's
   * outside the triangle to be consistent with preciselocate().
   */
  if (pdest[0] == searchpoint[0] && pdest[1] == searchpoint[1]) {
    /* TODO(Mario): modify searchtri and return ONVERTEX */
    return ONEDGE;
  }

  if (papex[0] == searchpoint[0] && papex[1] == searchpoint[1]) {
    /* TODO(Mario): modify searchtri and return ONVERTEX */
    return OUTSIDE;
  }

  /* 
   * Orientations of the searchpoint relative to the edges opposing the 
   * origin, destination and apex respectively
   */
  origorient = counterclockwise(m, b, papex, pdest, searchpoint);
  destorient = counterclockwise(m, b, porig, papex, searchpoint);
  apexorient = counterclockwise(m, b, pdest, porig, searchpoint);
  
  /* If the point lies on the outer side of either edge, return OUTSIDE */
  if (origorient > 0.0 || destorient > 0.0 || apexorient > 0.0) {
    /* TODO(Mario): set nexttri */
    return OUTSIDE;
  }

  /*
   * If point is on the edge opposing the apex, i.e. the primary edge,
   * claim it's ONEDGE
   */
  if (apexorient == 0.0) {
    return ONEDGE;
  }

  /* 
   * If point is on the edge opposing either the origin or destination, 
   * modify searchtri to point to the adjacent triangle and return ONEDGE
   */
  if (origorient == 0.00) {
    lnextself(*searchtri);
    return ONEDGE;
  }
  if (destorient == 0.00) {
    lprevself(*searchtri);
    return ONEDGE;
  }

  /* 
   * If the point is on the interior of this (now finite) triangle,
   * return INTRIANGLE.
   */
  return INTRIANGLE;
}

/**
 * Locates an element (triangle, edge, vertex), possibly infite, on which a 
 * a given point is situated by exhaustively searching through every triangle
 * in the triangulation.
 * This function is very expensive and should only be invoked if a hole is
 * encountered in the triangulation.
 *
 * Arguments:
 *   m: the mesh definition
 *   b: the triangulation behavior definition
 *   searchpoint: point whose triangle to locate
 *   searchtri: triangle which will contain searchpoint upon return
 *
 * Returns:
 *   locateresult value describing the type of element the point lies on
 **/
#ifdef ANSI_DECLARATORS
enum locateresult exhaustivelocate(struct mesh *m, struct behavior *b,
				   vertex searchpoint, struct otri *searchtri)
#else /* not ANSI_DECLARATORS */
  enum locateresult exhaustivelocate(m, b, searchpoint, searchtri)
     struct mesh *m;
     struct behavior *b;
     struct otri *searchtri;
#endif /* not ANSI_DECLARATORS */
{
  /* Initialize triangle traversal */
  traversalinit(&m->triangles);
  searchtri->orient = 0;

  enum locateresult location = OUTSIDE;
  /* Go through all live triangles in the mesh */
  for (searchtri->tri = triangletraverse(m);
       (location == OUTSIDE) && (searchtri->tri != (triangle *) NULL);
       searchtri->tri = triangletraverse(m)) {
    /* Locate point inside it and possibly its immediate neighbors */
    location = immediatelocate(m, b, searchpoint, searchtri, NULL);
  }
  /* 
   * Assuming infinite triangles, a point should never be on the outside
   * of a triangulation.
   */
  assert(location != OUTSIDE);
  return location;
}
