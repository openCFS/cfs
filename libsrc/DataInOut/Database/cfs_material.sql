-- MySQL dump 9.10
--
-- Host: localhost    Database: cfs_material
-- ------------------------------------------------------
-- Server version	4.0.18

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
  x double NOT NULL default '0',
  y double NOT NULL default '0',
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

