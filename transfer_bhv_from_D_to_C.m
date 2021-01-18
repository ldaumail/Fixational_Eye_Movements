%this script is intended to transfer all bhv files corresponding to the
%single units data analyzed in the previous analysis

indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
metafilename = load(strcat(indexdir, 'single_units_ns6_metadata.mat'));


Unique_sessions =uniqueStrCell(metafilename.ns6FileName);
 for i = 1:length(Unique_sessions)-1 
     %Big Drobo 2 - Backup\Drobo2\data\rig021
     %r1a\maierlab\DATA\NEUROPHYS\rig021\
      underscore = strfind(Unique_sessions{i+1}, '_');
      session = Unique_sessions{i+1};
      session =  session(1:underscore(2)-1);
     myFolder=strcat('D:\LGN_data_TEBA\rig022_2\', session);
     if exist(myFolder, 'dir')
mkdir(strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\bhv_selected_units\', session))
sourcedir=strcat('D:\LGN_data_TEBA\rig022_2\', session, '\', Unique_sessions{i+1});
destdir=strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\bhv_selected_units\',session);
     
status = copyfile( strcat(sourcedir,'.bhv'),  destdir);

     end
 end

 
 
 %% Find missing files 160609_I_cinterocdrft013 and 180809_I_cinterocdrft002
 
 %%160609_I_cinterocdrft013
 %session = '160609_I_cinterocdrft012';
 BRdatafile = '160609_I_cinterocdrft012';
 directory = 'D:\LGN_data_TEBA\rig021_2\160609_I\';
 filename   = [directory BRdatafile]; 
 eye_info = concatBHV(strcat(filename,'.bhv'));

 %directory = strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\bhv_selected_units\',session,'\');
 ns6dir = strcat(filename, '.ns6');
 bhvdir = strcat(filename,'.bhv');
[eye_dva,bhv,eye_ana] = getAnalogEyeinDVA(ns6dir,bhvdir);

BRdatafile = '160611_I_cinterocdrft001';
 directory = 'D:\LGN_data_TEBA\rig021_2\160611_I\';
 filename   = [directory BRdatafile]; 
 eye_info = concatBHV(strcat(filename,'.bhv'));

%[eye_dva,bhv,eye_ana] = getAnalogEyeinDVA(ns6dir,bhvdir);


 
 
 
 
