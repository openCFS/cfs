#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:38:49 2021

@author: alex

This converts multiharmonic openCFS-results into the timedomain.

This is used in 2 different ways:
    a) You can create a new CFS-file with time results out a MH-Result file, if you want to visulize it e.g. in Paraview
    For that tpye in terminal "mh2t.py --help"
    b) You can read in the mh-resultfile and obtain a class to work with that, if you want to postprocess in python
    For some example look into TESTSUIT/Singlefield/Electrostatics/compareSequenceStep.ipynb
"""
import sys
import os
# This file should be in /share/python/ where also hdf5_tools is, which is needed here
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import hdf5_tools as h5
import h5py
#%% Functions and subclasses

def f2t(X,Nt,return_Nh=False):
    '''
    Parameters
    ----------
    X : 1d array
    the harmonics vector.

    Nt : int
    Number of points in the time series.

    return_Nh : Boolean
    returns the number of harmonics computed by (len(X)-1)/2 as second output

    Returns
    -------
    x : 1d array, len(Nt)
    Time signal.

    Nh : int
    Number of harmonics in X
    '''
    from numpy import linspace, exp, outer, pi, arange
    t = linspace(0,2*pi,Nt,endpoint=False)
    if (X.shape[0]-1) % 2 > 0 :
        # print("Optimized version assumed!")
        X_new = np.zeros((2*X.shape[0]-1,), dtype=complex)
        # print(X_new)
        i=0
        for entry in X:
            #print(X)
            # print(i)
            #print("hello")
            # print(entry)
            # print(*entry)
            # print(entry[0]+entry[1])
            X_new[i] = entry
            i = i + 2
            
        X = X_new
    #     print('error harmonics: X.shape[0]-1 must be even. X.shape=',X.shape)
    #     return
    # else :
    Nh = int( (X.shape[0]-1)/2 )
    Xt = (exp(1j*outer(t,arange(-Nh,Nh+1))) @ X).real
    if return_Nh:
        return Xt, Nh
    else :
        return Xt

class Container:
    """
    A container saves either the frequency result (doFFT = False)
    Or saves the time result (doFFT = True)
    
    Container for frequeny domain or Time domain
    Frequency domain: [Nh,Elem]
    Time domain: [Nt,Elem]
    """
    def __init__(self,Name,NumDOFs, file,Nt, doFFT = False, region = None, multistep = 1):
        """
        Parameters
        ----------
        Name : String
            Used to obtain the result with hdf5_tools.get_results().
        NumDOFs : Integer
            Should be 1 (NodeResult) or 3 (ElementResult).
        file : string
            Path to file (for hdf5_tools.get_results()).
        Nt : Integer
            Number of points in the time series. Used for the FFT (f2t()) .
        doFFT : Boolean, optional
            Triggers the fft for the obtained frequency results . The default is False.
        region : string, optional
            For what region the result should be obtained. The default is None.
        multistep : string, optional
            For what multistep the result should be obtained. The default is 1.

        """

        self.__NumDOFs = NumDOFs
        self.__file = file
        self.__Nt=Nt
        self.__doFFT = doFFT
        
        if not self.__doFFT:
            
            if self.__NumDOFs ==1:
                self.val= h5.get_result(self.__file, Name, step = "all", region = region, multistep = multistep)
                self.__Nh= self.val.shape[0]
                
            elif self.__NumDOFs ==3:
                if len(h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep).shape) == 2:
                    # print("Exists only of 1 element")
                    self.x=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,0]
                    self.y=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,1]
                    self.z=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,2]
                    self.__Nh= self.x.shape[0]
                else:
                    # print("More elements in this region")
                    self.x=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,0]
                    self.y=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,1]
                    self.z=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,2]
                    self.__Nh= self.x.shape[0]

        else:
            if self.__NumDOFs ==1:
                val= h5.get_result(self.__file, Name, step = "all", region = region, multistep = multistep)
                try:
                    self.__nodeNum=val.shape[1]
                except:
                    self.__nodeNum=val.shape[0]
                    
                self.val=np.zeros(shape=(self.__Nt,self.__nodeNum))
                for node in range(self.__nodeNum):
                    try:
                        self.val[:,node] = f2t(val[:,node],self.__Nt)
                    except:
                        self.val[:,node] = f2t(val,self.__Nt)

                    
            elif self.__NumDOFs ==3:
                
                if len(h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep).shape) == 2:
                    x=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,0]
                    y=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,1]
                    z=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,2]
                    self.__elemNum=1
                    
                    self.x = np.zeros(shape=(self.__Nt, self.__elemNum))
                    self.y = np.zeros(shape=(self.__Nt, self.__elemNum))
                    self.z = np.zeros(shape=(self.__Nt, self.__elemNum))
                    # print(self.x.shape)
                    # print(x.shape)
                    for elem in range(self.__elemNum):
                        self.x[:,elem] = f2t(x,self.__Nt)
                        self.y[:,elem] = f2t(y,self.__Nt)
                        self.z[:,elem] = f2t(z,self.__Nt)

                else:
                    x=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,0]
                    y=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,1]
                    z=h5.get_result(self.__file,Name,step="all", region = region, multistep = multistep)[:,:,2]
                    self.__elemNum=x.shape[1]
                    
                    self.x = np.zeros(shape=(self.__Nt, self.__elemNum))
                    self.y = np.zeros(shape=(self.__Nt, self.__elemNum))
                    self.z = np.zeros(shape=(self.__Nt, self.__elemNum))                    
                
                    for elem in range(self.__elemNum):
                        self.x[:,elem] = f2t(x[:,elem],self.__Nt)
                        self.y[:,elem] = f2t(y[:,elem],self.__Nt)
                        self.z[:,elem] = f2t(z[:,elem],self.__Nt)
                        
    def Info(self):
        """
        Prints out how the values are stored in the container.
        """
        if self.__doFFT:
            if self.__NumDOFs==3:
                print(f"[Timesteps, Elem] = [{self.__Nt}, {self.__elemNum}]")
            else:
                print(f"[Timesteps, Nodes] = [{self.__Nt}, {self.__nodeNum}]")
        else:
            if self.__NumDOFs==3:
                print(f"[Nh, Elem] = [{self.__Nh}, {self.__elemNum}]")
            else:
                print(f"[Nh, Nodes] = [{self.__Nh}, {self.__nodeNum}]")

class Result:
    """
    Each Result gets a time-container and a frequency-container.
    The time-container is the fft from the frequency conatiner.
    """
    def __init__(self,Name,NumDOFs, DefinedOn, Nt, file, region = None, multistep = 1):
        self.__Name = Name
        self.__NumDOFs = NumDOFs
        self.__DefinedOn = DefinedOn
        self.__Nt = Nt
        self.__file=file
        self.freq = Container(self.__Name, self.__NumDOFs, self.__file , self.__Nt, region = region, multistep = multistep)
        self.time = Container(self.__Name, self.__NumDOFs, self.__file , self.__Nt, doFFT= True, region = region, multistep = multistep)

    def getInfo(self, printInfo=True):
        if printInfo:
            print(f"NumDOFS = {self.__NumDOFs}")
            print(f"Defined on {self.__DefinedOn}")
        return self.__NumDOFs, self.__DefinedOn

#%% Main class
class MH2T:
    def __init__(self, filePath , Nt = 128, multistep = 1, region = None):
        self.__multistep = multistep
        self.__filePath=filePath
        self.__f=h5py.File(self.__filePath,"r")
        self.__getNumberOfHarmonics()
        self.Nt = Nt   
        self.__region = region
        self.__getAllRegions()
        self.__readAllResults()
        
        
    def __getAllRegions(self):
        """
        Checks if only 1 region is in result file, otherwise error occurs
        """
        regionList3d =[]
        if self.__region == None:
            regionList = self.__f["Mesh/Regions"].keys()
        else:
            regionList = [self.__region]
            
        for region in regionList:
            if self.__f[f"Mesh/Regions/{region}"].attrs["Dimension"] == 3:
                regionList3d.append(region)
                
        if len(regionList3d) < 1:
            raise Exception("Please specify at least 1 Volumeregion")
        elif len(regionList3d) > 1:
            raise Exception("No region specified but more than one region present. Please choose 1 region.")
        else:
            self.__region3d = regionList3d[0]

    def __getNumberOfHarmonics(self):
        """
        gets the Number of harmonics. Is not needed for further calculations.
        """
        Nh = len([*self.__f[f"Results/Mesh/MultiStep_{self.__multistep}/"].keys()])-1
        self.Nh=Nh
        
    def __readResultInfo(self, result):
        """
        Searches for the Result and returns DOFs and where its definded on
        """
        
        try:
            # This function is need for the visit() function
            # Its kinda confusing but it works.
            def find_Key(Key):
                 """ Find first object with 'foo' anywhere in the name """
                 if result in Key:
                     return Key
                 
            path="Results/"+self.__f["Results"].visit(find_Key)
            NumDOFs = self.__f[path]["NumDOFs"][0]
        
            DefinedOn = self.__f[path]["DefinedOn"][0]
            if DefinedOn == 1:
                DefinedOn="node"
            elif DefinedOn == 4:
                DefinedOn="element"
            
            return NumDOFs, DefinedOn
        except:
            print(f"Result {result} not found. Be careful, its casesensitive!")
            
    def __readAllResults(self):
        """
        This function reads all results and creates for each result a class named like the Result).
        """
        # This function is need for the visit() function
        def find_ResultDescription(Key):
            """ Find first object with 'foo' anywhere in the name  for specified multistep"""
            # print(f"MultiStep_{self.__multistep}/ResultDescription")
            if f"MultiStep_{self.__multistep}/ResultDescription" in Key:
                return Key
        
        # used the h5.visit() function to find all results and store it in a list
        self.__AllResults = [*self.__f[self.__f.visit(find_ResultDescription)]]
        for result in self.__AllResults:
            NumDOFs, DefinedOn = self.__readResultInfo(result)
            tmpClass = Result(result, NumDOFs, DefinedOn, self.Nt, self.__f, region= self.__region3d, multistep = self.__multistep)
            setattr(self, result, tmpClass)
            
            
    def printAllResults(self):
        """
        Lists all results available
        """
        print(self.__AllResults)

    def createCFSTimeOutput(self):
        """
        This creates a new Result-File where the results are stored in Timedomain.
        """
        newFile = self.__filePath.split(".")[0]+"_2Time.cfs"
        try:
            ft=h5py.File(newFile,"a")
            for groups in ["Results", "Mesh"]:
                print(f"Start cloning {groups}.")
                self.__f.copy(groups, ft)
                
            print("Deleting exsisting harmonic steps")
            for harmonics in range(self.Nh):
                del ft[f"Results/Mesh/MultiStep_1/Step_{harmonics}"]
    
            print("Adjusting Result Description")
            
            string = "transient"
            byteString = bytes(string, 'ascii')
            ft["Results/Mesh/MultiStep_1"].attrs["AnalysisType"] = byteString
            ft["Results/Mesh/MultiStep_1"].attrs["LastStepNum"] = self.Nt
            ft["Results/Mesh/MultiStep_1"].attrs["LastStepValue"] = 1.0
            
            for result in self.__AllResults:
                del ft[f"Results/Mesh/MultiStep_1/ResultDescription/{result}/StepNumbers"]
                del ft[f"Results/Mesh/MultiStep_1/ResultDescription/{result}/StepValues"]
                
                ft.create_dataset(f"Results/Mesh/MultiStep_1/ResultDescription/{result}/StepNumbers", data=np.arange(1,self.Nt+1,1))
                ft.create_dataset(f"Results/Mesh/MultiStep_1/ResultDescription/{result}/StepValues", data=np.linspace(1/(self.Nt+1), 1,self.Nt))
    
            print("Creating timesteps")
            for step in range(1,self.Nt+1):
                ft.create_group(f"Results/Mesh/MultiStep_1/Step_{step}")
                ft[f"Results/Mesh/MultiStep_1/Step_{step}"].attrs["StepValue"]=1/self.Nt * step
                
            for result in self.__AllResults:
                print(f"Start copying {result} for {self.Nt} timesteps.")
                for step in range(0,self.Nt):
                    Obj = getattr(self, result)
                    numDOFs, DefinedOn = Obj.getInfo(printInfo=False)
                    if numDOFs == 1:
                        if DefinedOn.lower() == "element":
                            group = f"Results/Mesh/MultiStep_1/Step_{step+1}/{result}/Vol/Elements"
                            ft.create_group(group)
                            ft.create_dataset(group + "/Real", data = Obj.time.val[step,:].transpose())
                        else:
                            group = f"Results/Mesh/MultiStep_1/Step_{step+1}/{result}/Vol/Nodes"
                            ft.create_group(group)
                            ft.create_dataset(group + "/Real", data = Obj.time.val[step,:].transpose())
                    else:
                        if DefinedOn.lower() == "element":
                            group = f"Results/Mesh/MultiStep_1/Step_{step+1}/{result}/Vol/Elements"
                            ft.create_group(group)
                            ft.create_dataset(group + "/Real", data = np.array((Obj.time.x[step,:], Obj.time.y[step,:], Obj.time.z[step,:])).transpose())
                        else:
                            group = f"Results/Mesh/MultiStep_1/Step_{step+1}/{result}/Vol/Nodes"
                            ft.create_group(group)
                            ft.create_dataset(group + "/Real", data = np.array((Obj.time.x[step,:], Obj.time.y[step,:], Obj.time.z[step,:])).transpose())
                            
            print(f"Conversion finished. Saved file in same location as source under: \n {newFile}")
        except:
            print("#"*20)
            print("-"*20)
            print(f"Canceled the Process! Converted file already exists! Please delete {newFile} and try again.")
#%%

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser=argparse.ArgumentParser(description="Converts a multiharmonic CFS-result into a transient CFS-result. \r\nCurrently supports only one sequenceStep!", formatter_class=RawTextHelpFormatter)
    parser.add_argument("filename", help="Should be something like that: \npath/to/FileName.cfs or \n./FileName.cfs")
    parser.add_argument("-Nt", nargs="?", type=int ,default=128 , help="Numbers of timesteps for the FFT")
    args = parser.parse_args()
    filename = args.filename
    Nt = args.Nt
    print(f"Reading this file: {filename}")
    if filename == None:
        print("No filename specified")
    Obj = MH2T(filename,Nt=Nt)
    Obj.createCFSTimeOutput()