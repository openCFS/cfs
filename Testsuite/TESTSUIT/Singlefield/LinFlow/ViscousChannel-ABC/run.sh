#!/bin/sh


# mesh
trelis -batch -nographics -nojournal Thv_channel_PML.jou

#Thv-acou-PML
cfs_master Thermoviscous_PML

#Thv-acou-ABC 
cfs_master  Thermoviscous_acous_ABC

#Thv-ABC
cfs_master  ViscousChannel-ABC


