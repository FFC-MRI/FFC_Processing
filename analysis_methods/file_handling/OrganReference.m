classdef OrganReference < handle
    %OrganReference is a class designed to store reference data from
    %specific organs in order to automatise the selection of evolution
    %fields and time before starting an acquisition. 
    %   Detailed explanation goes here
    
    properties
        tissueList
        dataList                    % list of data used to generate the model
        dispersionModel             % fit object, provides T1 as a function of Bevo
        modelName                   % name of the model, for reference 
        fileName                    % path and name of the file
        userData
    end
    
    methods
        % constructor method:
        % The list of tissues is a Nx1 cell array containing the name of
        % the tissues inside the organ.
        % Each entry to the list of tissues corresponds to an entry in the
        % data list, which contains a list of dispersion curves obtained at
        % 37C for that tissue. Each dispersion curve is described by a 3xN
        % array where the first line is the magnetic field strength in T,
        % the second is the corresponding T1 in sec and the third is the
        % error on T1 in sec (0 if not known).
        % example:
        % organ.tissueList = {'Muscle';'Fat';'Blood';,...};
        % organ.dataList = {[8.5e6 5e6 2e6; 0.2 0.1 0.06; 0.01 0.005
        % 0.002],...};
        function organ = OrganReference
            organ.tissueList = {};
            organ.dataList = {};
            organ.dispersionModel = @(x) zeros(size(x));
            organ.modelName = 'empty model';
            organ.fileName = '';
            organ.userData = {};
        end
        
        % add a tissue type to the object
        function AddTissueType(organ)
            
        end
        
        % remove a tissue type from the object
        function RemoveTissueType(organ)
            
        end
        
        % add a dispersion curve to a tissue dataset
        function AddDispersion(organ)
            
        end
        
        % remove a dispersion curve from a tissue dataset
        function RemoveDispersion(organ)
            
        end
        
        function T1est = GetT1Estimation(organ,Bevo)
            T1est = organ.dispersionModel(Bevo);
        end
        
        % interrogate the object for the optimium evolution parameters to
        % look at a tissue type
        function [evolutionFieldList, evolutionTimeList] = GetBestEvolutionParam(organ,tissue,fieldNumber,timeNumber)
            
        end
        
    end
    
end

