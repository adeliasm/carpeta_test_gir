%% DIBUJAR MATRICES CONECTIVIDAD SUJETO INDIVIDUAL
% añado esta nueva línea para probar
% Deletes everything in the workspace.
clc
clear
close all
fclose all;

%% Configurable parameters

% Defines the paths.

config.path.plv     = '../data/MatricesConectividadDCLMEG/';
config.path.patt     = '*.mat'; 

config.subj = 18; %'NEMOS-001';

%% Script

% Load the connectivity matrices

% Gets the list of files.
files = dir(sprintf('%s%s',config.path.plv,config.path.patt));

if isempty(files)
    fprintf ( 1, 'No files found with pattern: "%s%s".\n',config.path.plv, config.path.patt );
end

% Load data for specific subject

findex = config.subj;  % findex = file index

    % Loads the data for this iteration's subject.
    plvdata = load(sprintf('%s%s',config.path.plv,files(findex).name));

    % Creates a large structure to load all subjects, generates the 'subject' field and fills it with the name for this iteration's subject.
    plvdatas(findex).subject = plvdata.subject;
    
    % In the first iteration of the loop, it creates a variable with the band names found in the original plv file. Also replaces '-' with '_' to avoid naming errors in Matlab in the struct fields.
     if ~ exist("band_list","var")
        band_list = strrep(cat(1,{plvdata.band.name}),'-','_');
                     % strrep -> find and replace substrings: strrep(str,old,new)
                     % cat -> concatenate arrays: cat(dim,A,B)
    end

    % Goes through each band inside this iteration's subject.
    for bindex = 1 : numel ( band_list )

        % Loads the plv matrix for this iteration's band.
            banddata  = plvdata.band( bindex ).plv_rms;
            heatmap(banddata, 'Colormap', parula) % draw graph connectivity matrix
            
          
    end
