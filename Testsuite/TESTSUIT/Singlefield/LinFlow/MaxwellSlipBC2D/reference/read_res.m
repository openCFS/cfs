function [t,val_cell,coords_cell] = read_res(fname,coord,coord_val,res_type,seqStep)
% This function reads hdf5 files and returns the results for a certain line
% specified with the coordinate as well as the value of this coordinate

    fig_width = 16;
    fig_height = 9;
    

    % general stuff
    seqStepStr = append("MultiStep_",string(seqStep));
    
    steps = h5readatt(fname,['/Results/Mesh/' char(seqStepStr)],'LastStepNum');
    StepNumbers = h5read(fname,['/Results/Mesh/' char(seqStepStr) '/ResultDescription/' res_type '/StepNumbers']);
    StepValues = h5read(fname,['/Results/Mesh/' char(seqStepStr) '/ResultDescription/' res_type '/StepValues']);

    regions = h5info(fname, ['/Mesh/Regions/']);

    regionGroups = regions.Groups;
    
    
    % search for the corresponding nodes
    
    if coord=="x"
        eval_row = 1;    
    elseif coord=="y"
        eval_row = 2;
    elseif coord=="z"
        eval_row = 3;
    else
        error("Unknown coordinate type!")
    end
    
    coords = transpose(h5read(fname,['/Mesh/Nodes/Coordinates']));
    
    nodes_to_read = find(abs(coords(:,eval_row)-coord_val)<1e-8);
    
    val_cell = cell(length(nodes_to_read),1);
    coords_cell = cell(length(nodes_to_read),1);
    
    % read in the valid regions for evaluation
    
    valid_regions = h5info(fname, ['/Results/Mesh/MultiStep_1/Step_1/' res_type]);

    valid_regions_groups = valid_regions.Groups;

    % read the results for each node match
    
    for pp=1:length(nodes_to_read)
        for ii=1:length(regionGroups)
            regionName = regionGroups(ii).Name;
            regionNodes = h5read(fname,[regionName '/Nodes']);

            nodeInd = find(regionNodes==nodes_to_read(pp));
            node_valid = 0;
            if ~isempty(nodeInd)
                disp(['Found node ' num2str(nodes_to_read(pp)) ' in region ' regionName])
                regionSplit = strsplit(regionName,'/');
                regionName_ = regionSplit{end};
                readInd = nodes_to_read(pp);
                
                for uu=1:length(valid_regions_groups)
                    validRegionName = valid_regions_groups(uu).Name;
                    validRegionSplit = strsplit(validRegionName,'/');
                    validRegionName_= validRegionSplit{end};
                    if strcmp(validRegionName_,regionName_)
                        node_valid = 1;
                    end
                end
            end
            
            if node_valid==1
                break
            end
        end

        if ~isempty(nodeInd) && node_valid==1
            disp(['Coordinates of the node: x=' num2str(coords(readInd,1)) 'm, y=' num2str(coords(readInd,2)) 'm, z=' num2str(coords(readInd,3)) 'm'])
        else
            error('Node not found!')
        end
        
        coords_cell{pp} = coords(readInd,:);

        disp('Reading result...')
        t = StepValues;
        val = zeros(length(StepValues),1);
        for ii=1:length(StepValues)
            stepVal = h5read(fname,['/Results/Mesh/MultiStep_1/Step_' num2str(ii) '/' res_type '/' regionName_ '/Nodes/Real']);
            val(ii) = stepVal(eval_row,nodeInd);
        end

        val_cell{pp} = val;
        disp('Finished!')

    end

end
