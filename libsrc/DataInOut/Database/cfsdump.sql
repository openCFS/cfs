-- MySQL dump 9.10
--
-- Host: localhost    Database: cfs
-- ------------------------------------------------------
-- Server version	4.0.18

--
-- Table structure for table `Dataset55`
--

CREATE TABLE Dataset55 (
  idx int(10) unsigned NOT NULL auto_increment,
  analysis_type int(11) NOT NULL default '0',
  dataset55_trans_prob_idx int(11) NOT NULL default '0',
  dataset55_norm_prob_idx int(11) NOT NULL default '0',
  dataset55_freq_prob_idx int(11) NOT NULL default '0',
  data_characteristic int(11) NOT NULL default '0',
  data_type int(11) NOT NULL default '0',
  values_per_node int(11) NOT NULL default '0',
  dataset55_norm_param_idx int(11) NOT NULL default '0',
  dataset55_trans_param_idx int(11) NOT NULL default '0',
  dataset55_freq_param_idx int(11) NOT NULL default '0',
  result_idx int(11) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_freq_param`
--

CREATE TABLE Dataset55_freq_param (
  idx int(10) unsigned NOT NULL auto_increment,
  freq_result_type int(11) NOT NULL default '0',
  freq_step_no int(11) NOT NULL default '0',
  frequency double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_freq_prob`
--

CREATE TABLE Dataset55_freq_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  freq_step int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_norm_param`
--

CREATE TABLE Dataset55_norm_param (
  idx int(10) unsigned NOT NULL auto_increment,
  mode_no int(11) NOT NULL default '0',
  freq double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_norm_prob`
--

CREATE TABLE Dataset55_norm_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type char(100) NOT NULL default '',
  mode_no int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_result`
--

CREATE TABLE Dataset55_result (
  Dataset55_idx int(11) unsigned NOT NULL default '0',
  result float NOT NULL default '0',
  node_no mediumint(11) unsigned NOT NULL default '0',
  node_idx smallint(11) unsigned NOT NULL default '0',
  PRIMARY KEY  (Dataset55_idx,node_no,node_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_trans_param`
--

CREATE TABLE Dataset55_trans_param (
  idx int(10) unsigned NOT NULL auto_increment,
  time_step int(11) NOT NULL default '0',
  time double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset55_trans_prob`
--

CREATE TABLE Dataset55_trans_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  time_step int(11) NOT NULL default '0',
  time_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56`
--

CREATE TABLE Dataset56 (
  idx int(10) unsigned NOT NULL auto_increment,
  analysis_type int(11) NOT NULL default '0',
  dataset56_trans_prob_idx int(11) NOT NULL default '0',
  dataset56_norm_prob_idx int(11) NOT NULL default '0',
  dataset56_freq_prob_idx int(11) NOT NULL default '0',
  data_characteristic int(11) NOT NULL default '0',
  data_type int(11) NOT NULL default '0',
  values_per_node int(11) NOT NULL default '0',
  dataset56_norm_param_idx int(11) NOT NULL default '0',
  dataset56_trans_param_idx int(11) NOT NULL default '0',
  dataset56_freq_param_idx int(11) NOT NULL default '0',
  result_idx int(11) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_freq_param`
--

CREATE TABLE Dataset56_freq_param (
  idx int(10) unsigned NOT NULL auto_increment,
  freq_result_type int(11) NOT NULL default '0',
  freq_step_no int(11) NOT NULL default '0',
  frequency double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_freq_prob`
--

CREATE TABLE Dataset56_freq_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  freq_step int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_norm_param`
--

CREATE TABLE Dataset56_norm_param (
  idx int(10) unsigned NOT NULL auto_increment,
  mode_no int(11) NOT NULL default '0',
  freq double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_norm_prob`
--

CREATE TABLE Dataset56_norm_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type char(100) NOT NULL default '',
  mode_no int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_result`
--

CREATE TABLE Dataset56_result (
  Dataset56_idx int(11) NOT NULL default '0',
  result double NOT NULL default '0',
  node_no int(11) NOT NULL default '0',
  node_idx int(11) NOT NULL default '0',
  PRIMARY KEY  (Dataset56_idx,node_no,node_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_trans_param`
--

CREATE TABLE Dataset56_trans_param (
  idx int(10) unsigned NOT NULL auto_increment,
  time_step int(11) NOT NULL default '0',
  time double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset56_trans_prob`
--

CREATE TABLE Dataset56_trans_prob (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  time_step int(11) NOT NULL default '0',
  time_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset666`
--

CREATE TABLE Dataset666 (
  idx int(10) unsigned NOT NULL auto_increment,
  solution_type int(10) unsigned NOT NULL default '0',
  analysis_type int(10) unsigned NOT NULL default '0',
  model_dimension int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset780`
--

CREATE TABLE Dataset780 (
  idx int(10) unsigned NOT NULL auto_increment,
  elem_label int(11) unsigned NOT NULL default '0',
  elem_type_geo tinyint(11) unsigned NOT NULL default '0',
  subtype tinyint(11) unsigned NOT NULL default '0',
  elem_grp_no mediumint(11) unsigned NOT NULL default '0',
  result_idx bigint(20) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset780_node`
--

CREATE TABLE Dataset780_node (
  dataset780_idx int(11) unsigned NOT NULL default '0',
  node int(11) unsigned NOT NULL default '0',
  PRIMARY KEY  (dataset780_idx,node)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset781`
--

CREATE TABLE Dataset781 (
  idx bigint(20) unsigned NOT NULL auto_increment,
  result_idx bigint(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Dataset781_node`
--

CREATE TABLE Dataset781_node (
  idx int(11) unsigned NOT NULL auto_increment,
  dataset781_idx int(10) unsigned NOT NULL default '0',
  node_label mediumint(11) unsigned NOT NULL default '0',
  x_coord float NOT NULL default '0',
  y_coord float NOT NULL default '0',
  z_coord float NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `InputParam`
--

CREATE TABLE InputParam (
  idx bigint(20) unsigned NOT NULL auto_increment,
  filename text NOT NULL,
  date_modified datetime NOT NULL default '0000-00-00 00:00:00',
  file longtext NOT NULL,
  no_references int(10) unsigned NOT NULL default '1',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_conduc_lin`
--

CREATE TABLE Material_conduc_lin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v11 double NOT NULL default '0',
  v12 double NOT NULL default '0',
  v13 double NOT NULL default '0',
  v21 double NOT NULL default '0',
  v22 double NOT NULL default '0',
  v23 double NOT NULL default '0',
  v31 double NOT NULL default '0',
  v32 double NOT NULL default '0',
  v33 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_coupling_lin`
--

CREATE TABLE Material_coupling_lin (
  v34 double NOT NULL default '0',
  v35 double NOT NULL default '0',
  v36 double NOT NULL default '0',
  material_idx bigint(20) unsigned NOT NULL default '0',
  v11 double NOT NULL default '0',
  v12 double NOT NULL default '0',
  v13 double NOT NULL default '0',
  v14 double NOT NULL default '0',
  v15 double NOT NULL default '0',
  v16 double NOT NULL default '0',
  v21 double NOT NULL default '0',
  v22 double NOT NULL default '0',
  v23 double NOT NULL default '0',
  v24 double NOT NULL default '0',
  v25 double NOT NULL default '0',
  v26 double NOT NULL default '0',
  v31 double NOT NULL default '0',
  v32 double NOT NULL default '0',
  v33 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_damping`
--

CREATE TABLE Material_damping (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v1 double NOT NULL default '0',
  v2 double NOT NULL default '0',
  v3 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_dielec_lin`
--

CREATE TABLE Material_dielec_lin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v11 double NOT NULL default '0',
  v12 double NOT NULL default '0',
  v13 double NOT NULL default '0',
  v21 double NOT NULL default '0',
  v22 double NOT NULL default '0',
  v23 double NOT NULL default '0',
  v31 double NOT NULL default '0',
  v32 double NOT NULL default '0',
  v33 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_hysteresis`
--

CREATE TABLE Material_hysteresis (
  material_idx varchar(100) NOT NULL default '',
  parameter enum('coupling_E','coupling_S','dielec_E') NOT NULL default 'coupling_E',
  saturation double NOT NULL default '0',
  remanence double NOT NULL default '0',
  coercitivity double NOT NULL default '0'
) TYPE=MyISAM;

--
-- Table structure for table `Material_magnet_lin`
--

CREATE TABLE Material_magnet_lin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v1 double NOT NULL default '0',
  v2 double NOT NULL default '0',
  v3 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_mech_lin`
--

CREATE TABLE Material_mech_lin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v11 double NOT NULL default '0',
  v12 double NOT NULL default '0',
  v13 double NOT NULL default '0',
  v14 double NOT NULL default '0',
  v15 double NOT NULL default '0',
  v16 double NOT NULL default '0',
  v21 double NOT NULL default '0',
  v22 double NOT NULL default '0',
  v23 double NOT NULL default '0',
  v24 double NOT NULL default '0',
  v25 double NOT NULL default '0',
  v26 double NOT NULL default '0',
  v31 double NOT NULL default '0',
  v32 double NOT NULL default '0',
  v33 double NOT NULL default '0',
  v34 double NOT NULL default '0',
  v35 double NOT NULL default '0',
  v36 double NOT NULL default '0',
  v41 double NOT NULL default '0',
  v42 double NOT NULL default '0',
  v43 double NOT NULL default '0',
  v44 double NOT NULL default '0',
  v45 double NOT NULL default '0',
  v46 double NOT NULL default '0',
  v51 double NOT NULL default '0',
  v52 double NOT NULL default '0',
  v53 double NOT NULL default '0',
  v54 double NOT NULL default '0',
  v55 double NOT NULL default '0',
  v56 double NOT NULL default '0',
  v61 double NOT NULL default '0',
  v62 double NOT NULL default '0',
  v63 double NOT NULL default '0',
  v64 double NOT NULL default '0',
  v65 double NOT NULL default '0',
  v66 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_nonlin`
--

CREATE TABLE Material_nonlin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  param enum('coupling_E','coupling_S','dielec_E','conduc_T','conduc_B','perm_T','perm_B','mag_B','mech_S') NOT NULL default 'coupling_E',
  matrix_index enum('11','12','13','14','15','16','21','22','23','24','25','26','31','32','33','34','35','36','41','42','43','44','45','46','51','52','53','54','55','56','61','62','63','64','65','66') NOT NULL default '11'
) TYPE=MyISAM;

--
-- Table structure for table `Material_parameter`
--

CREATE TABLE Material_parameter (
  idx bigint(20) unsigned NOT NULL auto_increment,
  mat_no int(11) NOT NULL default '0',
  mat_name varchar(100) NOT NULL default '',
  mechanic enum('uninitialized','linear','nonlinear_S','hysteresis_S') NOT NULL default 'uninitialized',
  coupling enum('uninitialized','linear','nonlinear_E','nonlinear_S','hysteresis_E','hysteresis_S') NOT NULL default 'uninitialized',
  dielectric enum('uninitialized','linear','nonlinear_E','hysteresis_E') NOT NULL default 'uninitialized',
  density double NOT NULL default '0',
  conductivity enum('uninitialized','linear','nonlinear_T','nonlinear_B','hysteresis_B','hysteresis_T') NOT NULL default 'uninitialized',
  permeability enum('uninitialized','linear','nonlinear_T','nonlinear_B','hysteresis_T','hysteresis_B') NOT NULL default 'uninitialized',
  magnetisation enum('uninitialized','linear','nonlinear_B','hysteresis_B') NOT NULL default 'uninitialized',
  compression double NOT NULL default '0',
  alpha double NOT NULL default '0',
  beta double NOT NULL default '0',
  diffusity double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Material_perm_lin`
--

CREATE TABLE Material_perm_lin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  v11 double NOT NULL default '0',
  v12 double NOT NULL default '0',
  v13 double NOT NULL default '0',
  v21 double NOT NULL default '0',
  v22 double NOT NULL default '0',
  v23 double NOT NULL default '0',
  v31 double NOT NULL default '0',
  v32 double NOT NULL default '0',
  v33 double NOT NULL default '0',
  PRIMARY KEY  (material_idx)
) TYPE=MyISAM;

--
-- Table structure for table `NodeHistory`
--

CREATE TABLE NodeHistory (
  idx bigint(20) unsigned NOT NULL auto_increment,
  node_no int(11) NOT NULL default '0',
  result_idx bigint(20) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `NodeHistoryValue`
--

CREATE TABLE NodeHistoryValue (
  nodehistory_idx bigint(20) unsigned NOT NULL default '0',
  dof mediumint(9) NOT NULL default '0',
  time double NOT NULL default '0',
  value double NOT NULL default '0'
) TYPE=MyISAM;

--
-- Table structure for table `Result`
--

CREATE TABLE Result (
  idx int(10) unsigned NOT NULL auto_increment,
  dataset666_idx int(10) unsigned NOT NULL default '0',
  stamp timestamp(14) NOT NULL,
  input_param_idx bigint(20) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

