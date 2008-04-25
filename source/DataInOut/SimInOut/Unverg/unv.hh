/*****************************************************************************
 * Header-File for the library libunv.a.                                     *
 *                                                                           *
 * Description: The library unv.a is used for manipulating UNIVERSAL-Files.  *
 *****************************************************************************/

#ifndef UNV_LIB_H
#define UNV_LIB_H

#include <fstream>

#include <def_use_unv.hh>

class unv {
  friend class dataset;
  enum UnvFlagType {
    UNVEOF=-1,
    UNVFalse=0,
    UNVTrue=1
  } ;
  std::ifstream* in_stream;
  std::ofstream* out_stream;
  std::streampos start;
  int       num_connections;
  int       odd_delimiter;
    
  int goto_pos(const std::streampos position) const;
  int find_delimiter(int* status);
  int find_next_set(std::streampos* position, int* status);
  int is_set_end(int* is_end);
public:
  unv() { in_stream=NULL, out_stream=NULL; }
  ~unv();
  enum open_mode {
    UNVRead,
    UNVWrite,
    UNVAppend
  };    
  int open(const char* dir, const char* name, open_mode =UNVRead);
  int close(void);
  int rewind(void);
};


class dataset {
  unv* u;
  std::streampos set_start;
protected:
  int set_num;
  enum open_mode { 
    SETRead,
    SETWrite,
    SETNot 
  };
  open_mode mode;

  int get_in_stream(std::ifstream** in_stream) const;
  int get_out_stream(std::ofstream** out_stream) const;
  int set_end(int* is_end) const { return u->is_set_end(is_end); }
public:
  dataset() { u=NULL, set_num=0; }
  virtual ~dataset() { }
  enum SetFlagType {
    SETEOF=-1,
    SETEmpty=0,
    SETTrue=1,
    SETEmptySet=2
  };
  int connect(unv* unv_file);
  void disconnect(void) { --u->num_connections, u=NULL, mode=SETNot; }
  int open_set(int* set_number);
  int rewind_set(void) const;
  int close_set(void);
  virtual int read_header(void);
  virtual int read_data(int* status);
  virtual int write_header(void);
  virtual int write_data(void);
};


struct data_15 {
  int   node_label;
  int   def_coord_sys_n;
  int   dis_coord_sys_n;
  int   color;
  double x1 ,x2, x3;
};


class dataset_15 : public dataset {
public:
  data_15 data;

  dataset_15() { set_num=15; }
  int read_data(int* status);
  //  int write_data(void);

};


struct header_55 {
  char   id[5][81];
  int    model_type;
  int    analysis_type;
  int    data_charact;
  int    spec_data_type;
  int    data_type;
  int    n_data_val_per_node;
  int    n_int_data_val;
  int    n_real_data_val;
  int    type_spec_int_par[10];
  double  type_spec_double_par[12];
};

struct data_55 {
  int    node_num;
  double* node_data;
};


class dataset_55 : public dataset {
  int curr_node_data_size;
public:
  header_55 header;
  data_55   data;
    
  dataset_55() { set_num=55, data.node_data=NULL; }
  virtual ~dataset_55() { delete[] data.node_data; } 
  int read_header(void);
  int read_data(int* status);

};


struct header_56 {
  char   id[5][81];
  int    model_type;
  int    analysis_type;
  int    data_charact;
  int    spec_data_type;
  int    data_type;
  int    n_data_val_per_element_pos;
  int    n_int_data_val;
  int    n_real_data_val;
  int    type_spec_int_par[10];
  double  type_spec_double_par[12];
};

struct data_56 {
  int    element_num;
  int    n_val_on_element;
  double* element_data;
};


class dataset_56 : public dataset {
  int curr_element_data_size;
public:
  header_56 header;
  data_56   data;
    
  dataset_56() { set_num=56, data.element_data=NULL; }
  virtual ~dataset_56() { delete[] data.element_data; }
  int read_header(void);
  int read_data(int* status);
  //  int write_header(void);
  //  int write_data(void);

};



struct entity {
  int entity_type;
  int entity_tag;
};


struct data_752 {
  int     group_n;
  int     act_const_set_n;
  int     act_rest_set_n;
  int     act_load_set_n;
  int     act_dof_set_n;
  int     n_entities;
  char    group_name[41];
  entity* entities;
};


class dataset_752 : public dataset {
  int curr_entity_size;
public:
  data_752 data;
    
  dataset_752() { set_num=752, data.entities=NULL; } 
  virtual ~dataset_752() { delete[] data.entities; } 
  int read_data(int* status);
  //  int write_data(void);

};


struct data_780 {
  int  element_label;
  int  fe_desc_id;
  int  phys_prob_tab_bin_n;
  int  phys_prob_tab_n;
  int  mat_prob_tab_bin_n;
  int  mat_prob_tab_n;
  int  color;
  int  n_nodes;
  int* node_labels;
};


class dataset_780 : public dataset {
  int curr_node_labels_size;
public:
  data_780 data;
    
  dataset_780() { set_num=780, data.node_labels=NULL; } 
  virtual ~dataset_780() { delete[] data.node_labels; } 
  int read_data(int* status);
  //  int write_data(void);

};


struct data_781 {
  int    node_label;
  int    def_coord_sys_n;
  int    dis_coord_sys_n;
  int    color;
  double x1, x2, x3;
};


class dataset_781 : public dataset {
public:
  data_781 data;
    
  dataset_781() { set_num=781; } 
  int read_data(int* status);
  //  int write_data(void);
};

// added by Simon Triebenbacher to get the dimension of the mesh

struct data_666 {
  int    num_nodes;
  int    num_elems;
  int    dim;
};

class dataset_666 : public dataset {
public:
  data_666 data;
    
  dataset_666() { set_num=666; } 
  int read_data(int* status);
  int get_dim() 
  {
    return data.dim;
  }
  //  int write_data(void);
};

#endif /* !UNV_LIB_H */

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
