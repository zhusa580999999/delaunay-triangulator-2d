#ifndef INTERNAL_IO_CIO_H
#define INTERNAL_IO_CIO_H

typedef enum {false, true} bool;

typedef enum {
  SP_VOID = -1,
  SP_FLOAT = 0,
  SP_INT = 1,
  SP_DOUBLE = 2
} SPdatatype;

typedef enum {
  SP_NONE = 0,
  SP_QUAD_TREE = 1,
  SP_OCT_TREE = 2,
  SP_CLARKSON_2D = 3,
  SP_CLARKSON_3D = 4
} SPfinalizemethod;


typedef enum {
  SM_ERROR = -1,
  SM_EOF = 0,
  SM_VERTEX = 1,
  SM_TRIANGLE = 2,
  SM_FINALIZED = 3
} SMevent;

typedef enum {
  SP_ERROR = -1,
  SP_EOF = 0,
  SP_POINT = 1,
  SP_FINALIZED_CELL = 2
} SPevent;

typedef struct SMwriter SMwriter;
typedef struct SPwriter SPwriter;

typedef struct SMreader SMreader;
typedef struct SPreader SPreader;

/*
  Vector copy
*/
void veccopy3iv(int v[3], const int a[3]);
void veccopy3fv(float v[3], const float a[3]);

/*
  SMwriter
*/
SMwriter *new_smwriter_sma();
SMwriter *new_smwriter_smb();
void delete_smwriter(SMwriter *writer);

void smwriter_set_nverts(SMwriter *writer, int nverts);
void smwriter_set_nfaces(SMwriter *writer, int nfaces);
void smwriter_set_boundingbox(SMwriter *writer, const float *bb_min_f, const float *bb_max_f);

void smwriter_write_vertex(SMwriter * writer, const float *v_pos_f);
void smwriter_write_triangle(SMwriter *writer, const int *t_idx);

bool smwriter_open(SMwriter *writer, FILE* file);
void smwriter_close(SMwriter *writer);


/*
  SPwriter
*/
SPwriter *new_spwriter_spa();
SPwriter *new_spwriter_spb();
void delete_spwriter(SPwriter *writer);

void spwriter_set_npoints(SPwriter *writer, int npoints);
void spwriter_set_boundingbox_f(SPwriter *writer, const float *bb_min_f, const float *bb_max_f);
void spwriter_set_datatype(SPwriter *writer, SPdatatype datatype);
void spwriter_set_finalizemethod(SPwriter *writer, SPfinalizemethod finalizemethod);

void spwriter_write_header(SPwriter *writer);
void spwriter_write_point_f(SPwriter *writer, const float *p_pos_f);
void spwriter_write_finalize_cell(SPwriter *writer, int idx);

bool spwriter_open(SPwriter *writer, FILE* file);
void spwriter_close(SPwriter *writer);

/*
  SMreader
*/
SMreader *new_smreader_sma();
SMreader *new_smreader_smb();
void delete_smreader(SMreader *reader);

int smreader_nverts(SMreader *reader);
int smreader_nfaces(SMreader *reader);
float *smreader_bb_min_f(SMreader *reader);
float *smreader_bb_max_f(SMreader *reader);

SMevent smreader_read_element(SMreader *reader);  
SMevent smreader_read_event(SMreader *reader);

float *smreader_v_pos_f(SMreader *reader);
int *smreader_t_idx(SMreader *reader);

bool smreader_open(SMreader *reader, FILE* file);
void smreader_close(SMreader *reader);

/*
  SPreader
*/
SPreader *new_spreader_spa();
SPreader *new_spreader_spb();
void delete_spreader(SPreader *reader);

int spreader_npoints(SPreader *reader);
float *spreader_bb_min_f(SPreader *reader);
float *spreader_bb_max_f(SPreader *reader);
SPfinalizemethod spreader_finalizemethod(SPreader *reader);

SPevent spreader_read_event(SPreader *reader);

float *spreader_p_pos_f(SPreader *reader);
int spreader_final_idx(SPreader *reader);
  
bool spreader_open(SPreader *reader, FILE* file, bool skip_finalize_header);
void spreader_close(SPreader *reader);


#endif  /* INTERNAL_IO_CIO_H */
