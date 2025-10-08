import numpy as np
import h5py
import shutil
import os

def PostProcessHarmonic25d(file_name, copy_file_name, pos_list = np.array([0.0]), create_copy = True, rho=1.205, K=1.41767e5):
    """
        Function to evaluate the inverse fourier transform at given position

        Input:
            file_name ... CFS result file that contains the wavenumber spectrum (from 2.5d harmonic simulation)
            copy_file_name ... Name of copied result file which is going to store the postprocessed results
            pos_list ... A list of positions for inverse fourier transform evaluation
            create_copy ... Store evaluated results in a new copied result file
            rho ... Density of air, default: 1.205 kg/m^3
            K ... Compression Modulus, default: 141767 Pa
    """
    ######### Calculate Speed of Sound based on Material Properties ########
    c0 = np.sqrt(K / rho)

    ######### Create a copy of original result file if True ########
    if create_copy:
        try:
            os.remove(copy_file_name)
        except OSError:
            pass
        shutil.copyfile(file_name, copy_file_name)
    ######### Create a copy of original result file if True########

    ######## open the file that should store postprocessed data ########
    h5_file = copy_file_name if create_copy else file_name
    h5_pos = h5py.File(h5_file, mode = 'r+') # open file with read write
    ######## open the file that should store postprocessed data ########

    ######## we have to rename the AnalysisType attribute to trick the paraview cfsreader ########
    myatr = h5_pos['Results/Mesh/MultiStep_1'].attrs
    if myatr.__getitem__('AnalysisType') == 'harmonic25d':
        myatr.__setitem__('AnalysisType',value='harmonic') # set attribute name to 'harmonic'
    ######## we have to rename the AnalysisType attribute to trick the paraview cfsreader ########

    ######## We create a new analysis step called 'MultiStep_2' to store the physical sound pressure results ########
     # first create a new group
    MultiStep2 = h5_pos['Results/Mesh'].create_group("MultiStep_2")
    # after creating the group, we also have to set the right attributes which will be read by the cfsreader in paraview
    MultiStep2.attrs.__setitem__('AnalysisType', value='harmonic')
    MultiStep2.attrs.__setitem__('LastStepNum', value=len(pos_list))
    MultiStep2.attrs.__setitem__('LastStepValue', value=pos_list[-1])
    # then we want to copy the ResultDescription part to the MultiStep_2
    # first get all the dataset that exists in MultiStep_1
    # we only want to copy the acouPressure part
    MultiStep2.create_group("ResultDescription/acouPressure")
    Datasetname = [name for name in h5_pos['Results/Mesh/MultiStep_1/ResultDescription/acouPressure']]
    for name in Datasetname:
        src = h5_pos['Results/Mesh/MultiStep_1/ResultDescription/acouPressure/%s' % name]
        dest = h5_pos['Results/Mesh/MultiStep_2/ResultDescription/acouPressure']
        MultiStep2.copy(source=src, dest = dest, name = name)
        #we also have to modify the StepNumbers and StepValues dataset
        if name == "StepNumbers":
            stepnum = h5_pos['Results/Mesh/MultiStep_2/ResultDescription/acouPressure/StepNumbers']
            stepnum.resize((len(pos_list),))
            stepnum[:] = np.arange(1,len(pos_list)+1)
        elif name == "StepValues":
            stepvalue = h5_pos['Results/Mesh/MultiStep_2/ResultDescription/acouPressure/StepValues']
            stepvalue.resize((len(pos_list),))
            stepvalue[:] = pos_list
    ######## We create a new analysis step called 'MultiStep_2' to store the physical sound pressure results ########

    ######## Perform Inverse Fourier Transform for each position z in the pos_list and store ########
    NumSteps = myatr.__getitem__("LastStepNum") # get the total number of wavenumber steps
    ### Preparation work for calculating inverse fourier transform ###
    # array storing the stepvalue_idx for reading pressure information
    stepvalue_idx = np.arange(1,NumSteps+1)
    # Extract wavenumber step values
    wavenumbers = np.empty(shape=(NumSteps,1))
    for i in range(1,NumSteps+1):
        step = h5_pos['Results/Mesh/MultiStep_1/Step_%d' % i]
        freq = step.attrs.__getitem__("StepValue") # we stored the frequency value in hdf5
        wavenumbers[i-1] = (2 * np.pi * freq) / c0 # we have to calculate the actual wavenumbers k = 2*pi*f / c0 c0...speed of sound
    # wavenumbers = np.concatenate((-1*np.flip(wavenumbers),wavenumbers))
    WaveNumberStepSize = wavenumbers[1] - wavenumbers[0]
    wavenumbers = wavenumbers.ravel()

    ### Loop through Regions
    for region in list(h5_pos['Results/Mesh/MultiStep_1/Step_1/acouPressure'].keys()):
        
        #### Integration using Trapezoidal Rule ####
        
        coeff = np.cos(pos_list[:, np.newaxis] * wavenumbers)

        FirstStep = h5_pos['Results/Mesh/MultiStep_1/Step_%d' % stepvalue_idx[0]]
        LastStep = h5_pos['Results/Mesh/MultiStep_1/Step_%d' % stepvalue_idx[-1]]
        PressureRealFirstStep = FirstStep['acouPressure/%s/Nodes/Real' % region][:].ravel()
        PressureRealLastStep = LastStep['acouPressure/%s/Nodes/Real' % region][:].ravel()
        PressureImagFirstStep = FirstStep['acouPressure/%s/Nodes/Imag' % region][:].ravel()
        PressureImagLastStep = LastStep['acouPressure/%s/Nodes/Imag' % region][:].ravel()

        SolReal = PressureRealFirstStep*coeff[:,0][:,np.newaxis] + PressureRealLastStep*coeff[:,-1][:,np.newaxis]
        SolImag = PressureImagFirstStep*coeff[:,0][:,np.newaxis] + PressureImagLastStep*coeff[:,-1][:,np.newaxis]
        
        for i in np.arange(1,len(stepvalue_idx)-1):
            step = h5_pos['Results/Mesh/MultiStep_1/Step_%d' % stepvalue_idx[i]]
            PressureReal = step['acouPressure/%s/Nodes/Real' % region][:].ravel()
            PressureImag = step['acouPressure/%s/Nodes/Imag' % region][:].ravel()
            SolReal += 2*PressureReal*coeff[:,i][:,np.newaxis]
            SolImag += 2*PressureImag*coeff[:,i][:,np.newaxis]

        SolReal *= WaveNumberStepSize / 2 / np.pi
        SolImag *= WaveNumberStepSize / 2 / np.pi

        #### Write Results To File ####
        for step_num, z in enumerate(pos_list,start=1):
            actStep = h5_pos['Results/Mesh/MultiStep_2'].create_group('Step_%d' % step_num) # create a new group
            actStep.attrs.__setitem__(name='StepValue', value = z) # store the z value to group attribute
            
            result_real = actStep.create_dataset(name = 'acouPressure/%s/Nodes/Real' % region, shape=SolReal[0].shape)
            result_imag = actStep.create_dataset(name = 'acouPressure/%s/Nodes/Imag' % region, shape=SolImag[0].shape)
            result_real[:] = SolReal[step_num-1]
            result_imag[:] = SolImag[step_num-1]
    ######## Perform Inverse Fourier Transform for each position z in the pos_list and store ########

    ######## don't forget to close the file afterwards ########
    h5_pos.close()
    ######## don't forget to close the file afterwards ########