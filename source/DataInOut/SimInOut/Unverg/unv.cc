/******************************************************************************
 * Implementation of the Library unv.a
 *****************************************************************************/

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>

#include "simInputUnv.hh"
#include "unv.hh"

unv::~unv() {

  if (num_connections>0)
  {
    if (in_stream) {
      in_stream->close();
      delete in_stream;
    }
    else if (out_stream) {
      out_stream->close();
      delete out_stream;
    }
  }
}


int unv::goto_pos(const std::streampos pos) const {

  if (in_stream)
    return in_stream->seekg(pos).good();
  else if (out_stream)
    return out_stream->seekp(pos).good();
  else {
    Warning("In Method unv::goto_pos: No open stream!");
    return 0;
  }

}


int unv::find_delimiter(int* status) {

#ifdef DEBUG
  LOG_TRACE(simInputUNV) << "In find_delimiter!";
#endif

  if (!in_stream) {
    CoupledField::Warning("In Method unv::find_delimiter: No std::ifstream is open!");
    return 0;
  }

  char buffer[81];
  char *p_c;

  while (in_stream->getline(buffer, 81).good()) {
    if ((p_c=strtok(buffer, " "))==buffer+4 && strncmp(p_c, "-1", 2)==0
        && strtok(NULL, " ")==NULL) {
      odd_delimiter=!odd_delimiter;
   #ifdef DEBUG
      LOG_TRACE(simInputUNV) << "Delimiter found";
   #endif
      *status=UNVTrue;
      return 1;
    }
  }

  if (in_stream->eof()) {
    *status=UNVEOF;
    return 1;
  }
  else if (in_stream->fail())
    CoupledField::Warning("In Method unv::find_delimiter: Input line is too long!");

  return 0;

}


int unv::find_next_set(std::streampos* pos, int* status) {

#ifdef DEBUG
  LOG_TRACE(simInputUNV) << "In find_next_set!";
#endif

  if (!in_stream) {
    CoupledField::Warning("In Method unv::find_next_set: No std::ifstream is open!");
    return 0;
  }

  int stat;

  if (odd_delimiter)
    if (!find_delimiter(&stat)) return 0;

  if (!find_delimiter(&stat)) return 0;

  if (stat==UNVEOF) {
    *status=UNVEOF;
    return 1;
  }

  *in_stream >> *status;

#ifdef DEBUG
  LOG_TRACE(simInputUNV) << "Set Num: " << *status;
#endif

  if (!in_stream->good()) {
    CoupledField::Warning("In Method unv::find_next_set: Could not read dataset id!");
    return 0;
  }

  *pos=in_stream->tellg();
  return 1;

}


int unv::is_set_end(int* is_end) {

  if (!in_stream) {
    CoupledField::Warning("In Method unv::is_set_end: No std::ifstream is open!");
    return 0;
  }

  char  buffer[81];
  char *p_c;
  *is_end=UNVFalse;
  std::streampos pos=in_stream->tellg();

  do
    in_stream->getline(buffer, 81);
  while (strtok(buffer, " \n")==NULL && in_stream->good());

  if (in_stream->eof()) {
    CoupledField::Warning("In Method unv::is_set_end: Reached end of file!");
    return 0;
  }

  if (in_stream->fail()) {
    CoupledField::Warning("In Method unv::is_set_end: Input line is too long!");
    return 0;
  }

  if ((p_c=strtok(buffer, " "))==buffer+4 && strcmp(p_c, "-1")==0
      && strtok(NULL, " ")==NULL)
    *is_end=UNVTrue;

  if (!in_stream->seekg(pos).good()) {
    CoupledField::Warning("In Method unv::is_set_end: Could not reset old std::streampos!");
    return 0;
  }

  return 1;

}


int unv::open(const char* dir, const char* name, open_mode mode) {

  if (in_stream || out_stream) {
    CoupledField::Warning("In Method unv::open: A stream is already open!");
    return 0;
  }

  char* file=new char[strlen(dir)+strlen(name)+1];
  strcpy(file, dir);
  strcat(file, name);

  int status;

  if (mode==UNVRead) {
    status=(in_stream=new std::ifstream(file,std::ios::binary))->good();
    start=in_stream->tellg();
  }
  else if (mode==UNVWrite) {
    status=(out_stream=new std::ofstream(file))->good();
    start=out_stream->tellp();
  }
  else {
    status=(in_stream=new std::ifstream(file, std::ios::app))->good();
    start=in_stream->tellg();
  }

  delete[] file;
  num_connections=0;
  odd_delimiter=0;

  return status;

}


int unv::close(void) {

  if (num_connections>0) {
    std::ostringstream stream;
    stream << "In Method unv::close: There are still " << num_connections \
           << " connections open!";
    CoupledField::Warning(stream.str().c_str());
    return 0;
  }

  if (in_stream) {
    in_stream->close();
    delete in_stream;
    in_stream=NULL;
  }
  else if (out_stream) {
    out_stream->close();
    delete out_stream;
    out_stream=NULL;
  }
  else {
    CoupledField::Warning("In Method unv::close: No stream is open!");
    return 0;
  }

  return 1;

}


int unv::rewind(void) {

  odd_delimiter=0;

  if (in_stream) {
    in_stream->clear();
    return in_stream->seekg(start).good();
  }
  else if (out_stream) {
    out_stream->clear();
    return out_stream->seekp(start).good();
  }
  else {
    CoupledField::Warning("In Method unv::goto_pos: No stream is open!");
    return 0;
  }

}



int dataset::get_in_stream(std::ifstream** in) const {

  if (!u) {
    CoupledField::Warning("In Method dataset::get_in_stream: No connection to UNV-File!");
    return 0;
  }

  if (u->in_stream) {
    *in=u->in_stream;
    return 1;
  }

  CoupledField::Warning("In Method dataset::get_in_stream: No std::ifstream is open!");

  return 0;

}


int dataset::get_out_stream(std::ofstream** out) const {

  if (!u) {
    CoupledField::Warning("In Method dataset::get_out_stream: No connection to UNV-File!");
    return 0;
  }

  if (u->out_stream) {
    *out=u->out_stream;
    return 1;
  }

  CoupledField::Warning("In Method dataset::get_out_stream: No std::ofstream is open!");

  return 0;

}


int dataset::connect(unv* p_unv) {

  if (u!=NULL)
    --u->num_connections;

  u=p_unv;
  ++u->num_connections;

  if (u->in_stream) {
    mode=SETRead;
    return 1;
  }
  else if (u->out_stream) {
    mode=SETWrite;
    return 1;
  }

  mode=SETNot;
  CoupledField::Warning("In Method dataset::connect: No stream is open in UNV!");

  return 0;

}


int dataset::open_set(int *status) {

  std::ofstream* out;
  std::streampos pos;
  int stat;

  if (mode==SETRead) {
    if (set_num==0) {
      if (!u->find_next_set(&pos, &stat)) return 0;
    }
    else {
      do
        if (!u->find_next_set(&pos, &stat)) return 0;
      while (stat!=unv::UNVEOF && stat!=set_num);
    }
    if (stat==unv::UNVEOF)
      *status=SETEOF;
    else {
      int empty;

      if (!set_end(&empty)) return 0;
      set_start=pos;
      if (empty==unv::UNVTrue)
        *status=SETEmpty;
      else
        *status=stat;
    }
    return 1;
  }
  else if (mode==SETWrite && get_out_stream(&out)) {
    if (u->odd_delimiter) {
      CoupledField::Warning("In Method dataset::open_set: A set is already open!");
      return 0;
    };
    *out << "    -1\n";
    u->odd_delimiter=1;
    *out << std::setw(6) << set_num << std::endl;
    if (!out->good()) {
      CoupledField::Warning("In Method dataset::open_set: Could not write delimiter!");
      return 0;
    }
    set_start=out->tellp();
    return 1;
  }

  CoupledField::Warning("In Method dataset::open_set: No stream is open!");

  return 0;

}


int dataset::rewind_set(void) const {

  if (!u->odd_delimiter) {
    CoupledField::Warning("In Method dataset::rewind_set: There is no open set!");
    return 0;
  }

  return u->goto_pos(set_start);

}


int dataset::close_set() {

  int stat;

  if (mode==SETRead) {
    std::ifstream *in;

    if (!get_in_stream(&in)) {
      CoupledField::Warning("In Method dataset::close_set: Could not get std::ifstream!");
      return 0;
    }
    if (!u->odd_delimiter) {
      CoupledField::Warning("In Method dataset::close_set: There is no open set!");
      return 0;
    }
    --u->num_connections;
    return u->find_delimiter(&stat);
  }

  std::ofstream* out;

  if (!get_out_stream(&out)) {
    CoupledField::Warning("In Method dataset::close_set: Could not get std::ofstream!");
    return 0;
  }

  *out << "    -1\n";
  u->odd_delimiter=0;

  if (!out->good()) {
    CoupledField::Warning("In Method dataset::close_set: Could not write delimiter!");
    return 0;
  }

  return 1;

}


int dataset::read_header(void) {
  std::ostringstream stream;
  stream << "In Method dataset::read_header: For dataset " << set_num  \
         << " is no header-record defined!";
  CoupledField::Warning(stream.str().c_str());
  return 0;
}


int dataset::read_data(int* n) {
  std::ostringstream stream;
  stream << "In Method dataset::read_data: For dataset " << set_num  \
         << " is no data-record defined!";
  CoupledField::Warning(stream.str().c_str());
  return 0;
}


int dataset::write_header(void) {
  std::ostringstream stream;
  stream << "In Method dataset::write_header: For dataset " << set_num  \
         << " is no header-record defined!";
  CoupledField::Warning(stream.str().c_str());
  return 0;
}


int dataset::write_data(void) {
  std::ostringstream stream;
  stream << "In Method dataset::write_data: For dataset " << set_num  \
         << " is no data-record defined!";
  CoupledField::Warning(stream.str().c_str());
  return 0;
}



int dataset_15::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_15::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_15::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.node_label >> data.def_coord_sys_n >> data.dis_coord_sys_n \
      >> data.color >> data.x1 >> data.x2 >> data.x3;

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_15::read_data: Error by reading data!");
    return 0;
  }

  int is_end;

  if (!set_end(&is_end)) return 0;

  *stat=SETTrue;
  if (is_end==SETTrue)
    *stat=SETEmpty;

  return 1;

}


/*
  int dataset_15::write_data(void) {

  if (mode==SETRead || mode==SETNot) {
  CoupledField::Warning("In Method dataset_15::write_data: Open mode is not WRITE!");
  return 0;
  }

  std::ofstream* out;

  if (!get_out_stream(&out)) {
  CoupledField::Warning("In Method dataset_15::write_data: Could not get std::ofstream!");
  return 0;
  }

  out->form("%10i%10i%10i%10i\n", data.node_label, data.def_coord_sys_n, \
  data.dis_coord_sys_n, data.color);
  out->form("%13.5e%13.5e%13.5e\n", data.x, data.y, data.z);

  if (!out->good()) {
  CoupledField::Warning("In Method dataset_15::write_data: Error by writing data!");
  return 0;
  }

  return 1;

  }
*/



int dataset_55::read_header(void) {

  int i;
  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_55::read_header: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_55::read_header: Could not get std::ifstream!");
    return 0;
  }

  do {
    in->getline(&header.id[1][0], 81);
    strcpy(&header.id[0][0], &header.id[1][0]);
  } while (strtok(&header.id[1][0], " \n")==NULL);

  for (i=1; i<5; ++i)
    in->getline(&header.id[i][0], 81);

#ifdef DEBUG
  if (in->good())
    LOG_TRACE(simInputUNV) << "Header of dataset 55 (first ID-line):" \
                           << &header.id[0][0];
#endif

  *in >> header.model_type >> header.analysis_type >> header.data_charact \
      >> header.spec_data_type >> header.data_type \
      >> header.n_data_val_per_node >> header.n_int_data_val \
      >> header.n_real_data_val;

  for (i=0; i<header.n_int_data_val; ++i)
    *in >> header.type_spec_int_par[i];

  for (i=0; i<header.n_real_data_val; ++i)
    *in >> header.type_spec_double_par[i];

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_55::read_header: Error by reading header!");
    return 0;
  }

  return 1;

}


int dataset_55::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_55::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_55::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.node_num;

  int type_size=2;

  if (header.data_type==2) type_size=1;

  int n=header.n_data_val_per_node*type_size;

  if (data.node_data==NULL) {
    data.node_data=new double[n];
    curr_node_data_size=n;
  }
  else if (n>curr_node_data_size) {
    delete[] data.node_data;
    data.node_data=new double[n];
    curr_node_data_size=n;
  }

  for (int i=0; i<n; ++i)
    *in >> data.node_data[i];

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_55::read_data: Error by reading data!");
    return 0;
  }

  int is_end;

  if (!set_end(&is_end)) return 0;

  *stat=SETTrue;
  if (is_end==SETTrue)
    *stat=SETEmpty;

  return 1;

}



int dataset_56::read_header(void) {
  int i;
  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_56::read_header: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_56::read_header: Could not get std::ifstream!");
    return 0;
  }

  do {
    in->getline(&header.id[1][0], 81);
    strcpy(&header.id[0][0], &header.id[1][0]);
  } while (strtok(&header.id[1][0], " \n")==NULL);

  for (i=1; i<5; ++i)
    in->getline(&header.id[i][0], 81);

#ifdef DEBUG
  if (in->good())
    LOG_TRACE(simInputUNV) << "Header of dataset 56 (first ID-line): " \
                           << &header.id[0][0];
#endif

  *in >> header.model_type >> header.analysis_type >> header.data_charact \
      >> header.spec_data_type >> header.data_type \
      >> header.n_data_val_per_element_pos >> header.n_int_data_val \
      >> header.n_real_data_val;

  for (i=0; i<header.n_int_data_val; ++i)
    *in >> header.type_spec_int_par[i];

  for (i=0; i<header.n_real_data_val; ++i)
    *in >> header.type_spec_double_par[i];

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_56::read_header: Error by reading header!");
    return 0;
  }

  return 1;

}


int dataset_56::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_56::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_56::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.element_num >> data.n_val_on_element;

  int is_end;

  if (data.element_num>0) {

    int type_size=2;

    if (header.data_type==2) type_size=1;

    int n=data.n_val_on_element*type_size;

    if (data.element_data==NULL) {
      data.element_data=new double[n];
      curr_element_data_size=n;
    }
    else if (n>curr_element_data_size) {
      delete[] data.element_data;
      data.element_data=new double[n];
      curr_element_data_size=n;
    }

    for (int i=0; i<n; ++i)
      *in >> data.element_data[i];

    if (!in->good()) {
      CoupledField::Warning("In Method dataset_56::read_data: Error by reading data!");
      return 0;
    }


    if (!set_end(&is_end)) return 0;

    *stat=SETTrue;
    if (is_end==SETTrue)
      *stat=SETEmpty;

    return 1;
  }
  else {
    if (!set_end(&is_end)) return 0;
    char buffer[256];
    in->getline(buffer, 81);
    in->getline(buffer, 81);
    *stat=SETEmptySet;
    return 0;
  }

}

int dataset_752::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_752::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_752::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.group_n >> data.act_const_set_n >> data.act_rest_set_n \
      >> data.act_load_set_n >> data.act_dof_set_n >> data.n_entities;

  char buffer[41];

  do {
    in->getline(buffer, 41);
    strcpy(data.group_name, buffer);
  } while (strtok(buffer, " \n")==NULL);

  int n=data.n_entities;

  if (data.entities==NULL) {
    data.entities=new entity[n];
    curr_entity_size=n;
  }
  else if (n>curr_entity_size) {
    delete[] data.entities;
    data.entities=new entity[n];
    curr_entity_size=n;
  }

  for (int i=0; i<data.n_entities; ++i)
    *in >> data.entities[i].entity_type \
        >> data.entities[i].entity_tag;

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_752::read_data: Error by reading data!");
    return 0;
  }

  int is_end;

  if (!set_end(&is_end)) return 0;

  *stat=SETTrue;
  if (is_end==SETTrue)
    *stat=SETEmpty;

  return 1;

}



int dataset_780::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_780::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_780::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.element_label;

  int is_end;

  if (data.element_label>0) {
    *in >> data.fe_desc_id      >> data.phys_prob_tab_bin_n  \
        >> data.phys_prob_tab_n >> data.mat_prob_tab_bin_n \
        >> data.mat_prob_tab_n  >> data.color >> data.n_nodes;

    int n=data.n_nodes;

    if (data.node_labels==NULL) {
      data.node_labels=new int[n];
      curr_node_labels_size=n;
    }
    else if (n>curr_node_labels_size) {
      delete[] data.node_labels;
      data.node_labels=new int[n];
      curr_node_labels_size=n;
    }

    for (int i=0; i<data.n_nodes; ++i)
      *in >> data.node_labels[i];

    if (!in->good()) {
      CoupledField::Warning("In Method dataset_780::read_data: Error by reading data!");
      return 0;
    }


    if (!set_end(&is_end)) return 0;

    *stat=SETTrue;
    if (is_end==SETTrue)
      *stat=SETEmpty;

    return 1;
  }
  else {
    if (!set_end(&is_end)) return 0;
    char buffer[256];
    in->getline(buffer, 81);
    in->getline(buffer, 81);
    *stat=SETEmptySet;
    return 0;
  }

}


int dataset_781::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_781::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_781::read_data: Could not get std::ifstream!");
    return 0;
  }

  *in >> data.node_label >> data.def_coord_sys_n >> data.dis_coord_sys_n \
      >> data.color >> data.x1 >> data.x2 >> data.x3;

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_781::read_data: Error by reading data!");
    return 0;
  }

  int is_end;

  if (!set_end(&is_end)) return 0;

  *stat=SETTrue;
  if (is_end==SETTrue)
    *stat=SETEmpty;

  return 1;

}


int dataset_666::read_data(int* stat) {

  if (mode==SETWrite || mode==SETNot) {
    CoupledField::Warning("In Method dataset_666::read_data: Open mode is not READ!");
    return 0;
  }

  std::ifstream* in;

  if (!get_in_stream(&in)) {
    CoupledField::Warning("In Method dataset_666::read_data: Could not get std::ifstream!");
    return 0;
  }
  /*
   *in >> data.node_label >> data.def_coord_sys_n >> data.dis_coord_sys_n \
   >> data.color >> data.x1 >> data.x2 >> data.x3;
  */
  /*
    std::stringstream line;
    char linebuf[1024];
    in->getline(linebuf, sizeof(linebuf));
    line << linebuf;
  */

  int dummy;
  *in >> dummy >> dummy >> data.dim >> data.num_nodes >> data.num_elems;

  if (!in->good()) {
    CoupledField::Warning("In Method dataset_666::read_data: Error by reading data!");
    return 0;
  }

  int is_end;

  if (!set_end(&is_end)) return 0;

  *stat=SETTrue;
  if (is_end==SETTrue)
    *stat=SETEmpty;

  return 1;

}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
