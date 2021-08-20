# -*- coding: utf-8 -*-
"""
Created on Wed May 26 08:17:44 2021

@author: Harald_Schuerz \n


The purpose of this script is to convert exported meshes from COMSOL (.mphtxt) to hdf5 for openCFS usage.
A short tutorial on how to use it can be found in the openCFS userdocu.
Tested with COMSOL 5.6, might also work with previous versions.


To convert the mesh open Anaconda prompt (or a Terminal), go to the directory where both the script (comsol_meshconvert.py)
and the meshfile (mesh_comsol.mphtxt) are and type "python comsol_meshconvert.py mesh_comsol.mphtxt".
After processing, the converted meshfile should be in the folder aswell, named like the .mphtxt-file,
but with a .cfs extension.


"""

import h5py
import numpy as np
import argparse

def read_node_coord(file):
    
    """
    function to read the node coordinates
    
    parameters:
        
        file: meshfile that was converted to a list
        
    returns:
        
        node_coord: array with the coordinates for all nodes
        
        coord_index_start: entry in file, at which the coordinate-section starts
        
        coord_index_stop: entry in file, at which the coordinate-section starts
    
    """
    #list for saving the coordinates
    node_coord=[]
    #read_nodes=0 means that the current line is not within the coordinate section  
    read_nodes=0                                
    
    
    
    for line in range(len(file)):
        
        if file[line]=="# Mesh vertex coordinates":
        #line is now in the coordinate section, read_nodes is set to 1                        
            coord_index_start=line+1
            read_nodes=1
        #number of vertices is read
            num_vert=int(file[line-3].split('#')[0].rstrip())
            coord_index_stop=coord_index_start+num_vert
        if read_nodes==1:
            for n in range(num_vert):
                #all lines of coordinates are read
                temp=[]
                for i in (file[n+coord_index_start].strip().split(" ")):
                    temp.append(float(i))
                #in 2D-cases we also need an additional z-coodinate (z=0)
                if len(temp)==2:
                    temp.append(0)
                node_coord.append(temp)
            read_labels=0
            vtx_found=0
            node_label_index_end=line-1 
            break
    #change the type to an array to directly write it to the hdf5 later
    node_coord=np.array(node_coord)             
    return node_coord,coord_index_start,coord_index_stop


def read_connectivity(file,start,regiondict):
    
    """
    
    function to read the connectivity for all supported element types
    
    parameters:
        
        file: meshfile that was converted to a list
        
        start: line, at which the reading starts
        
        regiondict: dict with a dict for each region that contains:
            name of the region
            index assigned by COMSOL
            dimension
            list of elements the region contains(empty)
            list of all nodes the region contains(empty)
    
    
    returns:
    
        connectivity: array of all connectivities (all elements)
        
        eltype_dict: dictionary with the openCFS-specific element-typenumber
        
        regiondict: dict with a dict for each region
        (in this function the list of elements each region contains and 
        list of all nodes each region contains is added)
        
        linepointer: the line that was last read 
        (so we do not have to start reading from the top again)
        
        maxdim: maximum dimension of the used elements
    
    
    """
    

    #dict of all CFS-elementtype-indexes for each element
    eltype_dict={}                                                          
    eltype_dict["cfs_typenum"]=[]                                           
    # the index that COMSOL gives the enteties of each selections (domains)
    # is not unique, there can be the same number
    # for entities in different selections, but the indexes are unique
    # for all selections with the same dimension.
    # so we need to take the dimension into account to get global uniqueness 
    eltype_dict["dimension"]=[]                                  
    #connectivity for each element
    connectivity=[]                
                                                                            
    maxnodes=0
    
    #line is currently not in the connectivity section
    read_conn=0
    
    #linepointer is updated
    linepointer=start
    #absolute element counter
    abs_elcount=0                    
    #highest order of dimension will be saved here
    maxdim=0
    
    
    for line in range(len(file)):
        #only start after the coordinate section
        if line >= linepointer:  
            #line is now at the type name                                          
            if (file[line].rstrip()[-9:]=="type name"): 
                #read the type (hex,quad)                   
                el_name=file[line].split('#')[0].rstrip().split(" ")[1].rstrip()
                #read number of nodes per element
                num_nodes=file[line+3].split('#')[0].rstrip()
                #complete element type (quad4,..)                                     
                el_type=el_name+num_nodes 
                # naming in COMSOL is different than in openCFS
                # -> find the corresponding element type in openCFS
                
                #4-node quad
                if el_type=="quad4":                                         #TODO: implement all types
                    dim=2
                    cfs_num=6
                    cfs_name="QUAD4"
                
                #3-node tria    
                elif el_type=="tri3":
                    dim=2
                    cfs_num=4
                    cfs_name="TRIA3"   
                
                #4-node tetraeder
                elif el_type=="tet4":
                    dim=3
                    cfs_num=9
                    cfs_name="TET4"
                
                #8-node hexaeder
                elif el_type=="hex8":
                    dim=3
                    cfs_num=11
                    cfs_name="HEXA8"
                #2 node line    
                elif el_type=="edg2":
                    dim=1
                    cfs_num=2
                    cfs_name="LINE2"
                #3 node line
                elif el_type=="edg3":
                    dim=1
                    cfs_num=3
                    cfs_name="LINE3"
                
                #6 node wedge (CFS), prism (COMSOL)
                elif el_type=="prism6":
                    dim=3
                    cfs_num=16
                    cfs_name="WEDGE6"
                    
                #5 node pyramid
                elif el_type=="pyr5":
                    dim=3
                    cfs_num=14
                    cfs_name="PYRA5"
                
                #elements that are not implemented (yet) are undefined    
                else:
                    dim=0
                    cfs_num=0
                    cfs_name="UNDEF"
                #we need to fill the shorter connectivities with zeros to have
                #equal length so we need to find the longest connectivity
                if int(num_nodes)>maxnodes:                                
                    maxnodes=int(num_nodes)                                 
                if dim>maxdim:
                    maxdim=dim
                #read number of elements
                num_els=file[line+4].split('#')[0].rstrip()
                #create dics entry for the type
                eltype_dict[el_type]=[]  
                #start of connectivity is reached, reading can start                                         
                read_conn=1   
                #6 lines after "type name" the values start                                              
                conn_index_start=line+6   
                #tells us at which line we soonest should look again
                #for "type name"                                  
                linepointer=conn_index_start+int(num_els)                   
            
            #only elements that have a CFS-equivalent are added
            if read_conn==1 and cfs_num!=0:                                 
                for i in range(int(num_els)):
                    #write the cfs_number for the current element down
                    eltype_dict["cfs_typenum"].append(cfs_num)
                    #write the dimension of the current element
                    eltype_dict["dimension"].append(dim)
                    #add the element number to the dict #TODO: is this needed? 
                    eltype_dict[el_type].append(abs_elcount)                             
                    #temporary list for the connectivity
                    temp=[]
                    #add the connectivity to the list
                    for v in file[i+conn_index_start].split(" "):           
                        #int(v)+1, because COMSOL numbers start at 0 which
                        #would be considered empty so the index is raised 
                        #by 1 for the node labels
                        temp.append(int(v)+1)
                        
                    #numbering convention of COMSOL is unusual, some entries 
                    #need to be switched to fit the openCFS convention
                    if cfs_num==6:
                        temp[2],temp[3]=temp[3],temp[2]
                    if cfs_num==11:
                       temp[2],temp[3],temp[6],temp[7]=temp[3],temp[2],temp[7],temp[6]
                    if cfs_num==14:
                        temp[2],temp[3]=temp[3],temp[2]
                    connectivity.append(temp)
                    #since the lists inside the list are not equally long, we need to convert each list                               
                    for k in regiondict.keys():
                        #if the element belongs to the current selection (k) and the dimension is equal
                        
                        for j in range(len(regiondict[k]['index'])):
                            if int(file[i+conn_index_start+int(num_els)+3].rstrip())==regiondict[k]['index'][j] and regiondict[k]['dimension']==dim:                #corresponding region index
                                #abs_elcount+1, because COMSOL numbers start at 0 which would be considered empty
                                #so the index is raised by 1 for the element labels
                                regiondict[k]['elementlist'].append(abs_elcount+1)
                                regiondict[k]['nodelist'].extend(temp)
                            
                    abs_elcount += 1
                read_conn=0
            elif read_conn==1:
                read_conn=0
    #Fill entries up with zeros
    for i in range(len(connectivity)):
        connectivity[i]=connectivity[i]+[0]*(maxnodes-len(connectivity[i]))     
                                                                                                                            
    #sorting the nodelists of the regions,
    #deleting double entries and convert to array
    for k in regiondict.keys():
        regiondict[k]['nodelist']=list(dict.fromkeys(regiondict[k]['nodelist']))
        #if the entry is the same as the one before it will be deleted
        
        regiondict[k]['nodelist']=np.array(regiondict[k]['nodelist'])    
    connectivity=np.array(connectivity) 
    eltype_dict["cfs_typenum"]=np.array(eltype_dict["cfs_typenum"])    
    return connectivity,eltype_dict,regiondict,linepointer,maxdim



def read_regions(file,start):
    
    """
    reads the regions that were defined in COMSOL
    
    parameters:
        
         file: meshfile that was converted to a list
        
         start: line, at which the reading starts
    
    returns:
        
        regions: dict with a dict for each region that contains:
            name of the region
            index assigned by COMSOL
            dimension
            list of elements the region contains (empty, filled in read_connectivity)
            list of all nodes the region contains (empty, filled in read_connectivity)
     
    """
    
    linepointer=start
    regions={}
    regnum=0
    for line in range(len(file)):
        if line >= linepointer:
            indexline=file[line].split('#')[0].rstrip().split(' ')
            #search for the assignments of the sections
            if indexline[-1]=="Selection":                              
                num_indices=int(file[line+5].rstrip().split('#')[0].rstrip())
                index=[]
                for i in range(num_indices):
                    index.append(int(file[line+7+i].rstrip()))
                dimension=int(file[line+4].rstrip().split('#')[0].rstrip())
                region_name_list=file[line+2].split('#')[-2].rstrip().split(' ')
                region_name=""
                for reg in region_name_list[1:]:
                    region_name+=reg.rstrip()
                #dict in dict
                regions[regnum]={}                                      
                regions[regnum]["region_name"]=region_name
                regions[regnum]["index"]=index
                regions[regnum]["dimension"]=dimension
                regions[regnum]["elementlist"]=[]
                regions[regnum]["nodelist"]=[]
                linepointer=line
                regnum+=1
                          
    return regions                                              
    
                       

def write_hdf5(filename,node_coord,connectivity,eltype_dict,regiondict,maxdim):
    
    """
    
    creates a hdf5-file and write all data to the hdf5-file
    
    parameters:
        
        filename: desired name of the hdf5-file (similar to the original name)
        
        node_coord: x,y,z - coordinates for each node
        
        connectivity: connectivity for all elements
        
        eltype_dict: dictionary with the openCFS-specific element-typenumber
        
        regiondict: dict with a dict for each region that contains:
            name of the region
            index assigned by COMSOL
            dimension
            list of elements the region contains
            list of all nodes the region contains
            
        maxdim: maximum dimension of the used elements
        
        returns:
            
            myfile: hdf5 - meshfile
    
    """
    
    #create hdf5-file
    new_name=filename.split('.')[0]+".cfs"
    myfile=h5py.File(new_name,"a")
    
    mesh=myfile.create_group("Mesh")
  
    mesh.attrs.create('Dimension', maxdim , dtype='uint32')       

    results=myfile.create_group("Results")                            
    userdata=myfile.create_group("UserData")
    fileinfo=myfile.create_group("FileInfo")
    
    fileinfo.create_dataset('Content', shape=(1,), dtype=np.uint32,data=np.array([1]), chunks=True, maxshape=(5,))    
    fileinfo.create_dataset('Creator',data='SH')
    fileinfo.create_dataset('Date',data='today')
    fileinfo.create_dataset('Version',data='0.9')
    
    #counters for elements
    num1D=0
    num2D=0
    num3D=0
    numelems=len(connectivity)
    numhex8=0
    numquad4=0
    numline2=0
    numline3=0
    numtet4=0
    numtria3=0
    numprism6=0
    numpyr5=0
    
    #counting elements
    for r in range(len(eltype_dict["cfs_typenum"])):
        
        if eltype_dict["cfs_typenum"][r]==2:
            numline2+=1
            num1D+=1
            
        elif eltype_dict["cfs_typenum"][r]==3:
            numline3+=1
            num1D+=1   
            
        elif eltype_dict["cfs_typenum"][r]==6:
            numquad4+=1
            num2D+=1  
            
        elif eltype_dict["cfs_typenum"][r]==11:
            numhex8+=1
            num3D+=1 
            
        elif eltype_dict["cfs_typenum"][r]==4:
            numtria3+=1
            num2D+=1 
            
        elif eltype_dict["cfs_typenum"][r]==9:
            numtet4+=1
            num3D+=1 
            
        elif eltype_dict["cfs_typenum"][r]==16:
            numprism6+=1
            num3D+=1 
            
        elif eltype_dict["cfs_typenum"][r]==9:
            numpyr5+=1
            num3D+=1 
        
            
    elements=mesh.create_group("Elements")
    elements.attrs.create('Num1DElems',num1D,dtype='uint32')                          #TO-DO: read numbers of each element fro COMSOl-output
    elements.attrs.create('Num2DElems',num2D,dtype='uint32')  
    elements.attrs.create('Num3DElems',num3D,dtype='uint32')  
    elements.attrs.create('NumElems',numelems,dtype='uint32')  
    elements.attrs.create('QuadraticElems',0)  
    elements.attrs.create('Num_HEXA8',numhex8,dtype='uint32')  
    elements.attrs.create('Num_QUAD4',numquad4,dtype='uint32')  
    elements.attrs.create('Num_HEXA20',0,dtype='uint32')  
    elements.attrs.create('Num_HEXA27',0,dtype='uint32')  
    elements.attrs.create('Num_LINE2',numline2,dtype='uint32')  
    elements.attrs.create('Num_LINE3',numline3,dtype='uint32') 
    elements.attrs.create('Num_POINT',0,dtype='uint32')  
    elements.attrs.create('Num_POLYGON',0,dtype='uint32')  
    elements.attrs.create('Num_POLYHEDRON',0,dtype='uint32')  
    elements.attrs.create('Num_PYRA13',0,dtype='uint32')  
    elements.attrs.create('Num_PYRA14',0,dtype='uint32')  
    elements.attrs.create('Num_PYRA5',numpyr5,dtype='uint32')  
    elements.attrs.create('Num_QUAD8',0,dtype='uint32')  
    elements.attrs.create('Num_QUAD9',0,dtype='uint32')  
    elements.attrs.create('Num_TET10',0,dtype='uint32')
    elements.attrs.create('Num_TET4',numtet4,dtype='uint32')
    elements.attrs.create('Num_TRIA3',numtria3,dtype='uint32')
    elements.attrs.create('Num_TRIA6',0,dtype='uint32')
    elements.attrs.create('Num_UNDEF',0,dtype='uint32')
    elements.attrs.create('Num_WEDGE15',0,dtype='uint32')
    elements.attrs.create('Num_WEDGE18',0,dtype='uint32')
    elements.attrs.create('Num_WEDGE6',numprism6,dtype='uint32')

    
    elements.create_dataset("Connectivity", data=connectivity,dtype='uint32')
    elements.create_dataset("Types", data=eltype_dict["cfs_typenum"],dtype='uint32')
    
    groups=mesh.create_group("Groups")
    
    nodes=mesh.create_group("Nodes")
    
    nodes.attrs.create('NumNodes',len(node_coord),dtype='uint32')                           
    coordinates=nodes.create_dataset("Coordinates",data=node_coord,dtype='float64')
    
    regions=mesh.create_group("Regions")
    for k in regiondict.keys():
        g=regions.create_group(regiondict[k]['region_name'])
        g.create_dataset("Elements",data=np.array(regiondict[k]["elementlist"]),dtype='uint32')
        g.create_dataset("Nodes",data=np.array(regiondict[k]["nodelist"]),dtype='uint32')
        g.attrs.create('Dimension',regiondict[k]["dimension"],dtype='uint32')
    myfile.close()
    return myfile

          
def read_file(filename):
    
    """
    
    read the meshfile
    
    parameters:
        
        filename: name of the exported mphxtx-file
        
    returns:
        
        templist: list where every entry is a line of the file
    
    """
    f=open(filename,'r')
    templist=[]
    for line in f:
        templist.append(line.rstrip())
    f.close()
    return templist
 


            


if __name__ == "__main__":
    
    
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filename")
    args = parser.parse_args()
    filename = args.filename
    print(filename)
   
    #read file
    file=read_file(filename)
    
    #read coordinates
    node_coordinates,coord_index_start,coord_index_stop=read_node_coord(file)
    
    #read regions and connectivity
    regiondict=read_regions(file,coord_index_stop)
    connectivity,eltype_dict,regiondict,conn_index_end,maxdim=read_connectivity(file,coord_index_start,regiondict)
 
    #create hdf5-file  
    write_hdf5(filename,node_coordinates,connectivity,eltype_dict,regiondict,maxdim)
    

   
    