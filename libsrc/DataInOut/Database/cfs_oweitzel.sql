-- MySQL dump 9.10
--
-- Host: localhost    Database: cfs_oweitzel
-- ------------------------------------------------------
-- Server version	4.0.18

--
-- Table structure for table `Calculation`
--

CREATE TABLE Calculation (
  idx bigint(20) unsigned NOT NULL auto_increment,
  stamp timestamp(14) NOT NULL,
  inputparam_idx int(20) unsigned NOT NULL default '0',
  solution_type int(11) NOT NULL default '0',
  analysis_type int(11) NOT NULL default '0',
  model_dimension int(11) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element`
--

CREATE TABLE Element (
  idx bigint(10) unsigned NOT NULL auto_increment,
  elem_label int(11) unsigned NOT NULL default '0',
  elem_type_geo tinyint(11) unsigned NOT NULL default '0',
  subtype tinyint(11) unsigned NOT NULL default '0',
  elem_grp_no mediumint(11) unsigned NOT NULL default '0',
  result_idx int(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_nodes`
--

CREATE TABLE Element_nodes (
  localidx mediumint(9) unsigned NOT NULL default '0',
  element_idx bigint(11) unsigned NOT NULL default '0',
  node int(11) unsigned NOT NULL default '0',
  PRIMARY KEY  (element_idx,localidx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_result`
--

CREATE TABLE Element_result (
  idx int(10) unsigned NOT NULL auto_increment,
  result_idx int(11) unsigned NOT NULL default '0',
  analysis_type int(11) NOT NULL default '0',
  data_characteristic int(11) NOT NULL default '0',
  data_type int(11) NOT NULL default '0',
  values_per_node int(11) unsigned NOT NULL default '0',
  element_result_norm_idx int(11) unsigned NOT NULL default '0',
  element_result_trans_idx int(11) unsigned NOT NULL default '0',
  element_result_freq_idx int(11) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_result_freq`
--

CREATE TABLE Element_result_freq (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  freq_step int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_result_norm`
--

CREATE TABLE Element_result_norm (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type char(100) NOT NULL default '',
  mode_no int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_result_trans`
--

CREATE TABLE Element_result_trans (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  time_step int(11) unsigned NOT NULL default '0',
  time_value double unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Element_result_value`
--

CREATE TABLE Element_result_value (
  element_result_idx int(11) unsigned NOT NULL default '0',
  elem_no int(11) unsigned NOT NULL default '0',
  dof tinyint(11) unsigned NOT NULL default '0',
  result double NOT NULL default '0',
  PRIMARY KEY  (element_result_idx,elem_no,dof)
) TYPE=MyISAM;

--
-- Table structure for table `InputParam`
--

CREATE TABLE InputParam (
  idx int(20) unsigned NOT NULL auto_increment,
  filename text NOT NULL,
  date_modified datetime NOT NULL default '0000-00-00 00:00:00',
  file longtext NOT NULL,
  no_references int(10) unsigned NOT NULL default '1',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Nodal_result`
--

CREATE TABLE Nodal_result (
  idx bigint(10) unsigned NOT NULL auto_increment,
  result_idx int(11) NOT NULL default '0',
  analysis_type int(11) NOT NULL default '0',
  data_characteristic int(11) NOT NULL default '0',
  data_type int(11) NOT NULL default '0',
  values_per_node int(11) NOT NULL default '0',
  nodal_result_norm_idx int(11) NOT NULL default '0',
  nodal_result_trans_idx int(11) NOT NULL default '0',
  nodal_result_freq_idx int(11) NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Nodal_result_freq`
--

CREATE TABLE Nodal_result_freq (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  freq_step int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Nodal_result_norm`
--

CREATE TABLE Nodal_result_norm (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type char(100) NOT NULL default '',
  mode_no int(11) NOT NULL default '0',
  freq_value double NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Nodal_result_trans`
--

CREATE TABLE Nodal_result_trans (
  idx int(10) unsigned NOT NULL auto_increment,
  result_type varchar(100) NOT NULL default '',
  time_step int(11) unsigned NOT NULL default '0',
  time_value double unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `Nodal_result_value`
--

CREATE TABLE Nodal_result_value (
  nodal_result_idx int(11) unsigned NOT NULL default '0',
  node_no int(11) unsigned NOT NULL default '0',
  dof tinyint(11) unsigned NOT NULL default '0',
  result double NOT NULL default '0',
  PRIMARY KEY  (nodal_result_idx,node_no,dof)
) TYPE=MyISAM;

--
-- Table structure for table `Node`
--

CREATE TABLE Node (
  idx int(20) unsigned NOT NULL auto_increment,
  result_idx int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `NodeHistory`
--

CREATE TABLE NodeHistory (
  idx bigint(20) unsigned NOT NULL auto_increment,
  node_no int(11) unsigned NOT NULL default '0',
  result_idx int(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

--
-- Table structure for table `NodeHistoryValue`
--

CREATE TABLE NodeHistoryValue (
  nodehistory_idx bigint(20) unsigned NOT NULL default '0',
  dof mediumint(9) unsigned NOT NULL default '0',
  time double unsigned NOT NULL default '0',
  value double NOT NULL default '0',
  PRIMARY KEY  (nodehistory_idx,dof,time)
) TYPE=MyISAM;

--
-- Table structure for table `Node_coordinates`
--

CREATE TABLE Node_coordinates (
  node_idx int(10) unsigned NOT NULL default '0',
  node_label int(11) unsigned NOT NULL default '0',
  x_coord float NOT NULL default '0',
  y_coord float NOT NULL default '0',
  z_coord float NOT NULL default '0',
  PRIMARY KEY  (node_idx,node_label)
) TYPE=MyISAM;

--
-- Table structure for table `Result`
--

CREATE TABLE Result (
  idx int(10) unsigned NOT NULL auto_increment,
  calculation_idx int(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (idx)
) TYPE=MyISAM;

