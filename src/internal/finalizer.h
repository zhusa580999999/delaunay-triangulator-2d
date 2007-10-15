#ifndef INTERNAL_FINALIZER_H
#define INTERNAL_FINALIZER_H

struct mesh;
struct behavior;
struct vertex;

#define BUFFER_SIZE 4096/sizeof(vertex)


/* region buffer to store points in the chunking step */
typedef struct region_buffer {
  unsigned int count;
  float buffer[BUFFER_SIZE][2];
  struct region_buffer *next;
} region_buffer;


/* a region */
typedef struct region {
  
  struct otri tri;
  
  unsigned long long timestamp;
  unsigned long long visited;
  
  unsigned int count;
  
  region_buffer *head;
  region_buffer **tailptr;

} region;

/* finalzier */
typedef struct finalizer {

  region *finalization_region;  /* Finalization region for point (temporary) */
  region *finalization_regions;                   /* Finalization region map */
  unsigned long long timestamp;                            /* Global timestamp */
	struct mesh coarse_mesh;
} finalizer;

// argc and argv are for writing the command line to the output files

void finalizer_finalize(struct mesh *m, struct behavior *b, int argc, char **argv);


#endif // INTERNAL_FINALIZER_H
