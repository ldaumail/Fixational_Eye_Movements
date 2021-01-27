%this script was written to extract the names of the concatenated ns6 files in
%which single units where isolated.
%Indeed, those .ns6 files bear the same name as the bhv files that contain
%the eye movements data

ns6dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\';

%ssdir
ssdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\kilosorted_files\';
%ssdir content
selected_name_list = dir(ssdir);  

%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']);

SingUnitConcatFiles = struct();
for pnt = 1:length(selected_trials_idx.logicals)
        if ~isempty(selected_trials_idx.logicals(pnt).idx)
        penetration_name = erase(selected_trials_idx.logicals(pnt).penetration, 'matDE50_NDE0');
        underscore = strfind(penetration_name, '_');
            for ss =1:length(selected_name_list)
                if contains(erase(selected_name_list(ss).name,'_ss.mat'),penetration_name(1:underscore(3)-1))
                    ss_file = load(strcat(ssdir,selected_name_list(ss).name));
                    ss_file_fieldnames = fieldnames(ss_file.ss);
                    cinterocdrft_names = ss_file_fieldnames(contains(ss_file_fieldnames,'cinterocdrft'));
                    %SingUnitConcatFiles.(strcat('x',penetration_name(1:underscore(2)-1))) = cinterocdrft_names;
                    SingUnitConcatFiles.(strcat('x',penetration_name)) = cinterocdrft_names;
             
                end
            end
        end
end
length(fieldnames(SingUnitConcatFiles))
save(strcat(indexdir,'concat_filenames_completenames.mat'),'-struct', 'SingUnitConcatFiles')