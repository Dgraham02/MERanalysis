% Purpose: allow for code testing/use on personal computers and workstation 
% interchangably without having to re-type all paths. 

function [file_path] = load_data_2018_03(USER)

    if contains(USER,'workstation')
        % Workstation Paths
        patient_folder = 'G:\Dakota\dbs-study\patient-2018-03';
        function_folder = 'G:\Dakota\dbs-study\global-functions';
        data_folder = 'G:\DBS Study Data';
        file_path = 'G:\DBS Study Data\2018_03_Ipsilateral.mat';
    elseif contains(USER,'dakota') 
        % Dakota's Paths
        patient_folder = 'C:\Users\Dakota\dbs-study\patient-2018-03';
        function_folder = 'C:\Users\Dakota\dbs-study\global-functions';
        data_folder = 'D:\DBS Study\DBS DATA';
        file_path = 'D:\DBS Study\DBS DATA\2018_03_Ipsilateral.mat';
    elseif contains(USER,'bradley')
        % Bradley's Paths
%         patient_folder =
%         function_folder =
%         data_folder =
%         file_path = 
    elseif contains(USER,'ashley')
        % Ashley's Paths
%         patient_folder =
%         function_folder =
%         data_folder =
%         file_path = 
    end

    % Generate paths to git repository & data 
    path_patient_folder = genpath(patient_folder);
    path_global_functions = genpath(function_folder);
    path_data_files = genpath(data_folder);

    % Add folders & subfolders to working path 
    addpath(path_patient_folder);
    addpath(path_global_functions);
    addpath(path_data_files);  
    
end 