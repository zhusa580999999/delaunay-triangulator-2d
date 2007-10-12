#include "smwriter_sma.h"
#include "smwriter_smb.h"

#include "spwriter_spa.h"
#include "spwriter_spb.h"

#include "smreader_sma.h"
#include "smreader_smb.h"

#include "spreader_spa.h"
#include "spreader_spb.h"

#include "vec3iv.h"
#include "vec3fv.h"

extern "C" {
  
  /*
    Vector copy
  */
  
  void veccopy3iv(int v[3], const int a[3]) {
    VecCopy3iv(v, a);
  }
    
  void veccopy3fv(float v[3], const float a[3]) {
    VecCopy3fv(v, a);
  }
  
  /*
    SMwriter
  */
  
  SMwriter *new_smwriter_sma() {
    return new SMwriter_sma();
  }
  
  SMwriter *new_smwriter_smb() {
    return new SMwriter_smb();
  }
  
  void delete_smwriter(SMwriter *writer) {
    delete writer;
  }
  
  void smwriter_set_nverts(SMwriter *writer, int nverts) {
    writer->set_nverts(nverts);
  }
    
  void smwriter_set_nfaces(SMwriter *writer, int nfaces) {
    writer->set_nfaces(nfaces);
  }
  
  void smwriter_set_boundingbox(SMwriter *writer, const float *bb_min_f, const float *bb_max_f) {
    writer->set_boundingbox(bb_min_f, bb_max_f);
  }
  
  void smwriter_write_vertex(SMwriter *writer, const float *v_pos_f) {
    writer->write_vertex(v_pos_f);
  }
  
  void smwriter_write_triangle(SMwriter *writer, const int *t_idx) {
    writer->write_triangle(t_idx);
  }
  
  bool smwriter_open(SMwriter *writer, FILE* file) {
    return writer->open(file);
  }
  
  void smwriter_close(SMwriter *writer) {
    writer->close();
  }
  
  /*
    SPwriter
  */
  
  SPwriter *new_spwriter_spa() {
    return new SPwriter_spa();
  }
  
  SPwriter *new_spwriter_spb() {
    return new SPwriter_spb();
  }
  
  void delete_spwriter(SPwriter *writer) {
    delete writer;
  }
  
  void spwriter_set_npoints(SPwriter *writer, int npoints) {
    writer->set_npoints(npoints);
  }
  
  void spwriter_set_boundingbox_f(SPwriter *writer, const float *bb_min_f, const float *bb_max_f) {
    writer->set_boundingbox(bb_min_f, bb_max_f);
  }
  
  void spwriter_set_datatype(SPwriter *writer, SPdatatype datatype) {
    writer->set_datatype(datatype);
  }
  
  void spwriter_set_finalizemethod(SPwriter *writer, SPfinalizemethod finalizemethod) {
    writer->set_finalizemethod(finalizemethod);
  }
  
  void spwriter_write_header(SPwriter *writer) {
    writer->write_header();
  }
  
  void spwriter_write_point_f(SPwriter *writer, const float *p_pos_f) {
    writer->write_point(p_pos_f);
  }
  
  void spwriter_write_finalize_cell(SPwriter *writer, int idx) {
    writer->write_finalize_cell(idx);
  }
  
  bool spwriter_open(SPwriter *writer, FILE* file) {
    return writer->open(file);
  }
  
  void spwriter_close(SPwriter *writer) {
    writer->close();
  }
  
  /*
    SMreader
  */
  
  SMreader *new_smreader_sma() {
    return new SMreader_sma();
  }
  
  SMreader *new_smreader_smb() {
    return new SMreader_smb();
  }
  
  void delete_smreader(SMreader *reader) {
    delete reader;
  }
  
  int smreader_nverts(SMreader *reader) {
    return reader->nverts;
  }
  
  int smreader_nfaces(SMreader *reader) {
    return reader->nfaces;
  }

  float *smreader_bb_min_f(SMreader *reader) {
    return reader->bb_min_f;
  }
  
  float *smreader_bb_max_f(SMreader *reader) {
    return reader->bb_max_f;
  }
  
  float *smreader_v_pos_f(SMreader *reader) {
    return reader->v_pos_f;
  }

  int *smreader_t_idx(SMreader *reader) {
    return reader->t_idx;
  }
  
  SMevent smreader_read_element(SMreader *reader) {
    return reader->read_element();
  }
  
  SMevent smreader_read_event(SMreader *reader) {
    return reader->read_event();
  }

  
  bool smreader_open(SMreader *reader, FILE* file) {
    return reader->open(file);
  }
  
  void smreader_close(SMreader *reader) {
    reader->close();
  }
  
  /*
    SPreader
  */


  SPreader *new_spreader_spa() {
    return new SPreader_spa();
  }

  SPreader *new_spreader_spb() {
    return new SPreader_spb();
  }

  void delete_spreader(SPreader *reader) {
    delete reader;
  }
  
  int spreader_npoints(SPreader *reader) {
    return reader->npoints;
  }

  float *spreader_bb_min_f(SPreader *reader) {
    return reader->bb_min_f;
  }

  float *spreader_bb_max_f(SPreader *reader) {
    return reader->bb_max_f;
  }
  
  SPfinalizemethod spreader_finalizemethod(SPreader *reader) {
    return reader->finalizemethod;
  }
   
  float *spreader_p_pos_f(SPreader *reader) {
    return reader->p_pos_f;
  }
  
  SPevent spreader_read_event(SPreader *reader) {
    return reader->read_event();
  }

  int spreader_final_idx(SPreader *reader) {
    return reader->final_idx;
  }
  
  bool spreader_open(SPreader *reader, FILE* file, bool skip_finalize_header) {
    return reader->open(file);
  }

  void spreader_close(SPreader *reader) {
    reader->close();
  }

}
