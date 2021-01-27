%this script is intended to transfer all bhv files from one hard drive to another.
%.bhv files correspond to the
%single units data analyzed in the previous analysis
%the file names of thos bhv files was retrieved using the script
%get_concat_filenames.m
%Last updated: 1/22/2021 by Loic Daumail
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
metafilename = load(strcat(indexdir, 'concat_filenames_completenames.mat'));

clusters = fieldnames(metafilename);

 for i = 1:length(fieldnames(metafilename))
     %Big Drobo 2 - Backup\Drobo2\data\rig021
     %r1a\maierlab\DATA\NEUROPHYS\rig021\
      %
      xcluster = clusters{i};
      underscore = strfind(xcluster, '_');
      session =  xcluster(2:underscore(2)-1);
     myFolder=strcat('D:\LGN_data_TEBA\rig022_2\', session);
     if exist(myFolder, 'dir')
mkdir(strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\', xcluster(2:end)))
        for nf =1:length(metafilename.(xcluster))
            xbhvfilename = metafilename.(xcluster){nf};
            bhvfilename = xbhvfilename(2:end);
            sourcedir=strcat('D:\LGN_data_TEBA\rig022_2\', session, '\', bhvfilename);
            destdir=strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\', xcluster(2:end));

            status = copyfile( strcat(sourcedir,'.bhv'),  destdir);
            disp(bhvfilename)
        end
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


 
 
 
 
