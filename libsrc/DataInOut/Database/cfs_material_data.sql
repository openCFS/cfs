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
-- Dumping data for table `Material_conduc_lin`
--

INSERT INTO Material_conduc_lin VALUES (3,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_conduc_lin VALUES (8,1000000,0,0,0,1000000,0,0,0,1000000);
INSERT INTO Material_conduc_lin VALUES (10,10000000,0,0,0,10000000,0,0,0,10000000);

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
-- Dumping data for table `Material_coupling_lin`
--

INSERT INTO Material_coupling_lin VALUES (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (5,0,0,0,0,12.3,0,0,0,0,12.3,0,0,-5.4,-5.4,15.2,0,0,0);
INSERT INTO Material_coupling_lin VALUES (6,0,0,0,0,12.3,0,0,0,0,12.3,0,0,-5.4,-5.4,15.2,0,0,0);
INSERT INTO Material_coupling_lin VALUES (7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (9,0,0,0,0,10.5,0,0,0,0,10.5,0,0,-4.1,-4.1,-4.1,0,0,0);
INSERT INTO Material_coupling_lin VALUES (10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
INSERT INTO Material_coupling_lin VALUES (12,0,0,0,0,10.5,0,0,0,0,10.5,0,0,-4.1,-4.1,-4.1,0,0,0);

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
-- Dumping data for table `Material_damping`
--


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
-- Dumping data for table `Material_dielec_lin`
--

INSERT INTO Material_dielec_lin VALUES (2,8.85419e-12,0,0,0,8.85419e-12,0,0,0,8.85419e-12);
INSERT INTO Material_dielec_lin VALUES (1,0.0885,0,0,0,0.0885,0,0,0,0.0885);
INSERT INTO Material_dielec_lin VALUES (5,7.836e-09,0,0,0,7.836e-09,0,0,0,6.188e-09);
INSERT INTO Material_dielec_lin VALUES (6,7.836e-09,0,0,0,7.836e-09,0,0,0,6.188e-09);
INSERT INTO Material_dielec_lin VALUES (7,6.19793e-11,0,0,0,6.19793e-11,0,0,0,6.19793e-11);
INSERT INTO Material_dielec_lin VALUES (3,8.85419e-12,0,0,0,8.85419e-12,0,0,0,8.85419e-12);
INSERT INTO Material_dielec_lin VALUES (9,7.12e-09,0,0,0,7.12e-09,0,0,0,5.84e-09);
INSERT INTO Material_dielec_lin VALUES (10,8.85419e-12,0,0,0,8.85419e-12,0,0,0,8.85419e-12);
INSERT INTO Material_dielec_lin VALUES (11,1,0,0,0,1,0,0,0,1);
INSERT INTO Material_dielec_lin VALUES (12,8.85419e-12,0,0,0,8.85419e-12,0,0,0,8.85419e-12);

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
-- Dumping data for table `Material_hysteresis`
--


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
-- Dumping data for table `Material_magnet_lin`
--

INSERT INTO Material_magnet_lin VALUES (3,0,0,0);
INSERT INTO Material_magnet_lin VALUES (8,0,0,0);
INSERT INTO Material_magnet_lin VALUES (10,0,0,0);

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
-- Dumping data for table `Material_mech_lin`
--

INSERT INTO Material_mech_lin VALUES (1,274400000000,108400000000,108400000000,0,0,0,108400000000,274400000000,108400000000,0,0,0,108400000000,108400000000,274400000000,0,0,0,0,0,0,83000000000,0,0,0,0,0,0,83000000000,0,0,0,0,0,0,83000000000);
INSERT INTO Material_mech_lin VALUES (2,165682000000,18409100000,18409100000,0,0,0,18409100000,165682000000,18409100000,0,0,0,18409100000,18409100000,165682000000,0,0,0,0,0,0,73636400000,0,0,0,0,0,0,73636400000,0,0,0,0,0,0,73636400000);
INSERT INTO Material_mech_lin VALUES (5,134600000000,57700000000,57700000000,0,0,0,57700000000,134600000000,57700000000,0,0,0,57700000000,57700000000,134600000000,0,0,0,0,0,0,38460000000,0,0,0,0,0,0,38460000000,0,0,0,0,0,0,38460000000);
INSERT INTO Material_mech_lin VALUES (6,134600000000,57700000000,57700000000,0,0,0,57700000000,134600000000,57700000000,0,0,0,0,0,134600000000,0,0,0,0,0,0,38460000000,0,0,0,0,0,0,38460000000,0,0,0,0,0,0,38460000000);
INSERT INTO Material_mech_lin VALUES (7,256944000000,90277800000,90277800000,0,0,0,90277800000,256944000000,90277800000,0,0,0,90277800000,90277800000,256944000000,0,0,0,0,0,0,83333300000,0,0,0,0,0,0,83333300000,0,0,0,0,0,0,83333300000);
INSERT INTO Material_mech_lin VALUES (3,256944000000,90277800000,90277800000,0,0,0,90277800000,256944000000,90277800000,0,0,0,90277800000,90277800000,256944000000,0,0,0,0,0,0,83333300000,0,0,0,0,0,0,83333300000,0,0,0,0,0,0,83333300000);
INSERT INTO Material_mech_lin VALUES (9,132000000000,71000000000,71000000000,0,0,0,71000000000,132000000000,71000000000,0,0,0,71000000000,71000000000,115000000000,0,0,0,0,0,0,26000000000,0,0,0,0,0,0,26000000000,0,0,0,0,0,0,26000000000);
INSERT INTO Material_mech_lin VALUES (10,165682000000,18409100000,18409100000,0,0,0,18409100000,165682000000,18409100000,0,0,0,18409100000,18409100000,165682000000,0,0,0,0,0,0,73636400000,0,0,0,0,0,0,73636400000,0,0,0,0,0,0,73636400000);
INSERT INTO Material_mech_lin VALUES (11,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2);
INSERT INTO Material_mech_lin VALUES (12,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1);

--
-- Table structure for table `Material_nonlin`
--

CREATE TABLE Material_nonlin (
  material_idx bigint(20) unsigned NOT NULL default '0',
  param enum('coupling_E','coupling_S','dielec_E','conduc_T','conduc_B','perm_T','perm_B','mag_B','mech_S') NOT NULL default 'coupling_E',
  matrix_index enum('11','12','13','14','15','16','21','22','23','24','25','26','31','32','33','34','35','36','41','42','43','44','45','46','51','52','53','54','55','56','61','62','63','64','65','66') NOT NULL default '11'
) TYPE=MyISAM;

--
-- Dumping data for table `Material_nonlin`
--


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
-- Dumping data for table `Material_parameter`
--

INSERT INTO Material_parameter VALUES (1,1,'steel','linear','linear','linear',7847.09,'uninitialized','uninitialized','uninitialized',0,10000,5e-10,0);
INSERT INTO Material_parameter VALUES (2,20,'silizium','linear','linear','linear',2300,'uninitialized','uninitialized','uninitialized',0,50,1e-07,0);
INSERT INTO Material_parameter VALUES (3,21,'air','linear','linear','linear',1.205,'linear','linear','linear',141767,0,0,0);
INSERT INTO Material_parameter VALUES (4,22,'water','uninitialized','uninitialized','uninitialized',1000,'uninitialized','uninitialized','uninitialized',2250000000,0,0,0);
INSERT INTO Material_parameter VALUES (5,23,'isoMat','linear','linear','linear',0,'uninitialized','uninitialized','uninitialized',0,0,0,0);
INSERT INTO Material_parameter VALUES (6,24,'P4','linear','linear','linear',0,'uninitialized','uninitialized','uninitialized',0,0,0,0);
INSERT INTO Material_parameter VALUES (7,25,'SiN','linear','linear','linear',2700,'uninitialized','uninitialized','uninitialized',0,0,0,0);
INSERT INTO Material_parameter VALUES (8,26,'iron','uninitialized','uninitialized','uninitialized',0,'linear','linear','linear',0,0,0,0);
INSERT INTO Material_parameter VALUES (9,27,'pzt-4','linear','linear','linear',7502.55,'uninitialized','uninitialized','uninitialized',0,0,0,0);
INSERT INTO Material_parameter VALUES (10,28,'alu','linear','linear','linear',2300,'linear','linear','linear',0,50,1e-07,0);
INSERT INTO Material_parameter VALUES (11,29,'springmat2','linear','linear','linear',7502.55,'uninitialized','uninitialized','uninitialized',0,80000,6e-10,0);
INSERT INTO Material_parameter VALUES (12,30,'springmat','linear','linear','linear',7502.55,'uninitialized','uninitialized','uninitialized',0,80000,6e-10,0);

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
-- Dumping data for table `Material_perm_lin`
--

INSERT INTO Material_perm_lin VALUES (3,1.25664e-06,0,0,0,1.25664e-06,0,0,0,1.25664e-06);
INSERT INTO Material_perm_lin VALUES (8,0.000125664,0,0,0,0.000125664,0,0,0,0.000125664);
INSERT INTO Material_perm_lin VALUES (10,1.25664e-06,0,0,0,1.25664e-06,0,0,0,1.25664e-06);

