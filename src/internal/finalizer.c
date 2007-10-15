/* Internal data structures and function declarations */
#include "internal/triangle.h"
#include "internal/finalizer.h"
#include <assert.h>

typedef void (*dfs_visitor)(finalizer *f, struct behavior *b, vertex p, region *rgn, unsigned int depth);
typedef void (*triangle_visitor)(finalizer *f, struct behavior *b, struct otri tri, int number);

/* finalizer algorithms */
void finalizer_init(finalizer *f, struct behavior *b);
void finalizer_do_coarse(finalizer *f, struct behavior *b);
void finalizer_do_counter(finalizer *f, struct behavior *b);
void finalizer_do_chunking(finalizer *f, struct behavior *b);
void finalizer_deinit(finalizer *f, struct behavior *b);

/* bowyer_watson algorithm */
void finalizer_bowyer_watson(finalizer *f, struct behavior *b, vertex p, dfs_visitor previsit, dfs_visitor postvisit);
void finalizer_bowyer_watson_dfs(finalizer *f, struct behavior *b, vertex p, struct otri searchtri, unsigned int depth, dfs_visitor previsit, dfs_visitor postvisit);

/* bowyer_watson DFS visitors */
void finalizer_counter_visit(finalizer *f, struct behavior *b, vertex p, region *rgn, unsigned int depth);
void finalizer_chunking_pre_visit(finalizer *f, struct behavior *b, vertex p, region *rgn, unsigned int depth);
void finalizer_chunking_post_visit(finalizer *f, struct behavior *b, vertex p, region *rgn, unsigned int depth);
void finalizer_set_triangle_id(finalizer *f, struct behavior *b, struct otri tri, int id);

/* writes coarse mesh */
void finalizer_write_coarsemesh(finalizer *f, struct behavior *b); /* write_header() */
void finalizer_write_coarsemesh_vertices(finalizer *f, struct behavior *b); /* write_header_vertices() */
void finalizer_write_coarsemesh_triangles(finalizer *f, struct behavior *b); /* write_header_triangles() */

/* writes point stream */
void finalizer_write_pointstream_header(finalizer *f, struct behavior *b, int points); /* write_points_header() */
void finalizer_write_pointstream_finalized_point(finalizer *f, struct behavior *b, region *rgn); /* spwriter_write_point_f() */
void finalizer_write_pointstream_finalized_tag(finalizer *f, struct behavior *b, region *rgn); /* spwriter_write_finalize_cell() */


/****************************************************************
    Finalizer's main functions. Do coarse, counter, and chunking steps to
    generate coarse mesh (.sma file) and fine point streams with finalization tags
    (.spa file)
*/
void finalizer_finalize(struct mesh *m, struct behavior *b, int argc, char **argv)
{


}
