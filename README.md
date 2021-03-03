# Fixational_Eye_Movements
## This repository is intended to analyze small fixational eye movements and the impact of these on rapid neural adaptation
### microsaccade data was first transferred on daily work harddrive from the overall LGN data on storage harddrive with the script "transfer_bhv_from_D_to_C.m"
### "Microsaccades_analysis.m" was the first data processing file created. It requires the clustering algorithm toolbox from Otero-Millan, 2014 (JOV), and is intended to isolate microsaccades (onsets/offsets/amplitude/velocity/acceleration...)
         -Data location:-
         #### indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
         #### selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices of the relevant files (with filenames)
         #### concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames (obtained with selected_trials_idx file)
         #### newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
         #### trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']); %neural data
         #### C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
 
        


### "trigger_neuraldat_to_msacc.m" was subsequently created to plot the neural data around microsaccade onset times 
