%This script attempts to align the neural data to microsaccade onset/offset times
%Written by Loic Daumail -02-04-2021

monodatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
neurfilename = [monodatadir 'su_peaks_03032020_corrected\all_units\clean_origin_sup_50'];
neuralDat = load(strcat(neurfilename));
indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices
% peak locations
locsfilename = [monodatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_locs'];
peakLocs = load(strcat(locsfilename));

condSelectedTrialsIdx = struct();
condNeuralDat =struct();
condPeakLocs = struct();
for n =1:length(neuralDat.clean_origin_data)
    if ~isempty(neuralDat.clean_origin_data(n).unit)
        condSelectedTrialsIdx.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = selected_trials_idx.logicals(n).idx;
        %use selected_trial_idx to get cluster filename for each single
        %unit so the neural data is identifiable
        condNeuralDat.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = neuralDat.clean_origin_data(n).unit;
        condPeakLocs.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = peakLocs.peaks_locs(n).locs;
    end
end
eyeMovDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\data\';
eyeMovDat = load( [eyeMovDir, 'all_eye_movement_data']); %eye movement data
 
%% neural activity on the peaks before and after microsaccade onset


%1) get eye movement locations of interest and store them relative to each
%peak
xfilenames = fieldnames(eyeMovDat);
allSelectSaccLocs = struct();
for i =1:length(fieldnames(eyeMovDat))
    xcluster = xfilenames{i};
    trialindex = condSelectedTrialsIdx.(xcluster);
    selectSaccLocs = nan(6,6,length(trialindex));
    negSaccLocs = nan(6,6,length(trialindex));
    posSaccLocs = nan(6,6,length(trialindex));
    %ControlselectSaccLocs = nan(6,6,length(trialindex));
    for tr = 1:length(trialindex)
        clear trSaccLocs trPeakLocs
        if isfield(eyeMovDat.(xcluster),sprintf('t%d',trialindex(tr)))
            %get eye movement location (specify eye movement below 2 degrees..?)
            saccades = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).saccades;
            enum = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).enum;
            trSaccLocs = saccades(:,enum.startIndex)+200;
            %select eye movement locations within the peak spans (lets start
            %with +-125)
            trPeakLocs = condPeakLocs.(xcluster)(:,tr);%-200; %need to subtract 200 to peak location (or add 200 to msacc location, e.g. trSaccLocs) as the peak locations were estimated from trials that start 200 ms before stim onset (in get_clean_peaks_and_data.m)
            for pl =1:length(trPeakLocs)
                cntneg =0;
                cntpos =0;
                cnt =0;
                amplitudes = saccades(:,enum.amplitude);
                for sl=1:length(trSaccLocs)
                    if amplitudes(sl) >= 0.5 && amplitudes(sl) <= 2
                        
                        if (trPeakLocs(pl)-75 <= trSaccLocs(sl)) && (trSaccLocs(sl)<= trPeakLocs(pl)+75)
                            cnt =cnt+1;
                            selectSaccLocs(cnt,pl,tr) = trSaccLocs(sl); %dim1 =msacc location, dim2 = peak number, dim3=trial
                            %control
                            % rand = randi([trPeakLocs(pl)-75, trPeakLocs(pl)+75],1,1);
                            %ControlselectSaccLocs(cnt,pl,tr) = rand;
                        end
                        
                        if (trPeakLocs(pl)-75 <= trSaccLocs(sl)) && (trSaccLocs(sl)<= trPeakLocs(pl))
                            cntneg =cntneg+1;
                            
                            negSaccLocs(cntneg,pl,tr) = trSaccLocs(sl); %dim1 =msacc location, dim2 = peak number, dim3=trial
                            
                        end
                        if (trSaccLocs(sl)> trPeakLocs(pl)) && (trSaccLocs(sl)<= trPeakLocs(pl)+75)
                            cntpos =cntpos+1;
                            
                            posSaccLocs(cntpos,pl,tr) = trSaccLocs(sl); %dim1 =msacc location, dim2 = peak number, dim3=trial
                            
                        end
                    end
                end
                
            end
        end
        
    end
    allSelectSaccLocs.(xcluster).all = selectSaccLocs;
    allSelectSaccLocs.(xcluster).neg = negSaccLocs;
    allSelectSaccLocs.(xcluster).pos = posSaccLocs;
    if isfield(eyeMovDat.(xcluster), 'cellclass')
        allSelectSaccLocs.(xcluster).cellclass =  eyeMovDat.(xcluster).cellclass;
    else
        allSelectSaccLocs.(xcluster).cellclass = [];
    end
end

%allSelectSaccLocs.(xcluster) = ControlselectSaccLocs;



%2) Align neural data on microsaccade onset from  peak1 to peak 4
%Here, neural data is aligned on every microsaccade onset. This means that
%the whole trial is entirely stored for each peak respectively. If multiple
%saccades occur during a single peak, the trial is also stored multiple
%times (once per microsaccade)
aligned_trials = struct();
for i =1:length(fieldnames(eyeMovDat))
    xcluster = xfilenames{i};
    if find(~all(isnan(allSelectSaccLocs.(xcluster).pos)))
        nDat = condNeuralDat.(xcluster)(401:1900,:);
        trialindex = condSelectedTrialsIdx.(xcluster);
        
        %MsaccLocked(:,ms,pn,n) =
        
        up_dist_trials = nan(6,4,length(trialindex));%dim1 =msacc location, dim2 = peak number, dim3=trial
        clear pn ms
        for pn = 1:4
            for ms =1:6
                locs_msacc = allSelectSaccLocs.(xcluster).pos(ms,pn, :);
                up_dist_trials(ms,pn,:)= length(nDat(:,1))- locs_msacc;
            end
        end
        %get the max distance between the saccalign and the stimulus onset
        max_low_dist_unit = max(allSelectSaccLocs.(xcluster).pos,[],'all');
        %create new matrix with the length(max(d)+max(xabs - d))
        new_dist_unit = max_low_dist_unit + max(up_dist_trials,[],'all');
        fp_locked_trials = nan(new_dist_unit,length(nDat(1,:)),6,4);
        %filtered_fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)),4);
        clear n ms pn
        for pn =1:4
            for ms =1:6
                for n = 1:length(nDat(1,:))
                    if ~isnan(allSelectSaccLocs.(xcluster).pos(ms,pn,n))
                        lower_unit_bound =max_low_dist_unit-allSelectSaccLocs.(xcluster).pos(ms,pn,n)+1;
                        upper_unit_bound =max_low_dist_unit-allSelectSaccLocs.(xcluster).pos(ms,pn,n)+length(nDat(:,1));
                        fp_locked_trials(lower_unit_bound:upper_unit_bound,n,ms,pn) = nDat(:,n);
                    end
                end
                
                %mean_origin_dSUA.(xcluster)(ms,pn)=  nanmean(fp_locked_trials(:,:,ms,pn),2);
            end
        end
        %get the aligned data if it exists for the unit
        aligned_trials.(xcluster)= fp_locked_trials;
        max_low_dist.(xcluster) = max_low_dist_unit;
        mat_mld(i) = max_low_dist_unit;
    end
end


%%store all trials of a given peak, in a matrix, across all units
%clear aligned_trials
trials_dat = nan(251,6,4,41, length(fieldnames(aligned_trials)));%Dim1 =length data that we wanna look at (+-125ms for each msacc) on 4 peaks,Dim2: 41 = max number of trials found in a unit in the variable aligned_trials
%control
control_trials_dat = nan(251,6,4,41, length(fieldnames(aligned_trials)));
min_mld = min(mat_mld);
max_mld = max(mat_mld);

clear i
for i = 1:length(fieldnames(aligned_trials))
    xcluster = xfilenames{i};
    if isfield(aligned_trials, xcluster)
        %if ~isnan(max_low_dist.(xcluster))
        for tr = 1:length(aligned_trials.(xcluster)(1,:,1,1))
            for pn = 1:4
                for ms =1:6
                    %if ~all(isnan(aligned_trials.(xcluster)(:,tr,ms,pn)))
                    if ~isempty(find(aligned_trials.(xcluster)(:,tr,ms,pn),1))
                        % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
                        trials_dat(1:251,ms,pn,tr,i)= aligned_trials.(xcluster)(max_low_dist.(xcluster)-1-124:max_low_dist.(xcluster)+125,tr,ms,pn);
                        %control
                        %rand = randi([min_mld, 1375],1,1);
                        %trials_dat(1:251,ms,pn,tr,i)= aligned_trials.(xcluster)(rand-1-124:rand+125,tr,ms,pn);
                    end
                end
            end
        end
        %normalizing with max and min of each unit
        %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
        
    end
    
end

%% Plots
%plot mean of single units across all peaks
reshTrialsDat = reshape(trials_dat,[251, 6*4*41,length(fieldnames(aligned_trials))]);
for i =1:42

xcluster = xfilenames{i};
%msMean = squeeze(nanmean(trials_dat,2));
%pnMean = squeeze(nanmean(msMean,2));
%trMean = squeeze(nanmean(pnMean,2));
overallMean = nanmean(reshTrialsDat(:,:,i),2);
hi_ci = overallMean + 1.96*std(reshTrialsDat(:,:,i),[],2,'omitnan')/sqrt(length(~all(isnan(reshTrialsDat(1,:,i)))));
low_ci = overallMean - 1.96*std(reshTrialsDat(:,:,i),[],2,'omitnan')/sqrt(length(~all(isnan(reshTrialsDat(1,:,i)))));
figure();
plot([-125:125],overallMean)
hold on
h1= ciplot( hi_ci, low_ci,[-125:125],[40/255 40/255 40/255],0.1);
set(h1, 'edgecolor','none') 
xlabel('time (ms)')
ylabel('Spike rate (spikes/sec)')
h2 = vline(0);
set(h2(1),'linewidth',1);
title(strcat(xfilenames{i}, sprintf(' Cell class: %s',allSelectSaccLocs.(xcluster).cellclass)),'Interpreter', 'none')

%saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\overall_mean_',sprintf('%s.png',xcluster(2:9))));
end


%just to compare the data processing that gathers only 251ms of the trials
%together
for i =1:42
figure();
xcluster = xfilenames{i};
msMean = squeeze(nanmean(aligned_trials.(xcluster),3));
pnMean = squeeze(nanmean(msMean,3));
trMean = nanmean(pnMean,2);
plot(trMean(max_low_dist.(xcluster)-1-124:max_low_dist.(xcluster)+125))
xlabel('time (ms)')
ylabel('Spike rate (spikes/sec)')
vline(125)
title(xfilenames{i},'Interpreter', 'none')
end
%% Plot overall mean
%reshTrialsDat = reshape(control_trials_dat,[251, 6*4*41,length(fieldnames(aligned_trials))]);
NeurDat2D = reshape(reshTrialsDat(:,:,24:end),[251,6*4*41*19]); %only 2018 and over units
%NeurDat2D = reshape(reshTrialsDat(:,:,24:end),[251,6*4*41*18]);%after
%peaks, we loose one unit
%NeurDat2D = reshape(reshTrialsDat,[251,6*4*41*42]); %all units
overallMean = nanmean(NeurDat2D,2);
hi_ci = overallMean + 1.96*std(NeurDat2D,[],2,'omitnan')/sqrt(length(~all(isnan(NeurDat2D(1,:)))));
low_ci = overallMean - 1.96*std(reshTrialsDat(:,:,i),[],2,'omitnan')/sqrt(length(~all(isnan(NeurDat2D(1,:)))));
figure();
plot([-125:125],overallMean)
hold on
h1= ciplot( hi_ci, low_ci,[-125:125],[40/255 40/255 40/255],0.1);
set(h1, 'edgecolor','none') 
xlabel('time (ms)')
ylabel('Spike rate (spikes/sec)')
h2 = vline(0);
set(h2(1),'linewidth',1);
title('Overall mean across all units, Msaccs After Peak','Interpreter', 'none')
saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\overall_mean_across_all_units_post2018_msaccs_afterpeak.png'));


%% Test code for just one unit
% neural activity on the peaks before and after microsaccade onset


%1) get eye movement locations of interest and store them relative to each
%peak
xfilenames = fieldnames(eyeMovDat);
%allSelectSaccLocs = struct();
%for i =1:length(fieldnames(eyeMovDat))
xcluster = xfilenames{31};
trialindex = condSelectedTrialsIdx.(xcluster);
selectSaccLocs = nan(6,6,length(trialindex));

for tr = 1:length(trialindex)
    clear trSaccLocs trPeakLocs
    if isfield(eyeMovDat.(xcluster),sprintf('t%d',trialindex(tr)))
        %get eye movement location (specify eye movement below 2 degrees..?)
        saccades = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).saccades;
        enum = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).enum;
        trSaccLocs = saccades(:,enum.startIndex)+200;
        %select eye movement locations within the peak spans (lets start
        %with +-125)
        trPeakLocs = condPeakLocs.(xcluster)(:,tr);%-200; %need to subtract 200 to peak location (or add 200 to msacc location, e.g. trSaccLocs) as the peak locations were estimated from trials that start 200 ms before stim onset (in get_clean_peaks_and_data.m)
        for pl =1:length(trPeakLocs)
            cnt =0;
            for sl=1:length(trSaccLocs)
                if (trPeakLocs(pl)-75 <= trSaccLocs(sl)) && (trSaccLocs(sl)<= trPeakLocs(pl)+75)
                 cnt =cnt+1;   
                    selectSaccLocs(cnt,pl,tr) = trSaccLocs(sl); %dim1 =msacc location, dim2 = peak number, dim3=trial
                end
            end
        end

    end

end
allSelectSaccLocs.(xcluster) = selectSaccLocs;
%end


%2) Align neural data on microsaccade onset from  peak1 to peak 4
aligned_trials = struct();
%for i =1:length(fieldnames(eyeMovDat))
   % xcluster = xfilenames{i};
    nDat = condNeuralDat.(xcluster)(601:1900,:);
    trialindex = condSelectedTrialsIdx.(xcluster);
    
    %MsaccLocked(:,ms,pn,n) = 
    
    up_dist_trials = nan(6,4,length(trialindex));%dim1 =msacc location, dim2 = peak number, dim3=trial
    clear pn ms
    for pn = 1:4
        for ms =1:6
            locs_msacc = allSelectSaccLocs.(xcluster)(ms,pn, :);
            up_dist_trials(ms,pn,:)= length(nDat(:,1))- locs_msacc;
        end
    end
    %get the max distance between the saccalign and the stimulus onset
    max_low_dist_unit = max(allSelectSaccLocs.(xcluster),[],'all');
    %create new matrix with the length(max(d)+max(xabs - d))
    new_dist_unit = max_low_dist_unit + max(up_dist_trials,[],'all'); 
    fp_locked_trials = nan(new_dist_unit,length(nDat(1,:)),6,4);
    %filtered_fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)),4);
     clear n ms pn
     for pn =1:4
         for ms =1:6
           for n = 1:length(nDat(1,:))
               if ~isnan(allSelectSaccLocs.(xcluster)(ms,pn,n))
                  lower_unit_bound =max_low_dist_unit-allSelectSaccLocs.(xcluster)(ms,pn,n)+1;
                  upper_unit_bound =max_low_dist_unit-allSelectSaccLocs.(xcluster)(ms,pn,n)+length(nDat(:,1));
                  fp_locked_trials(lower_unit_bound:upper_unit_bound,n,ms,pn) = nDat(:,n); 
               end
           end

     %mean_origin_dSUA.(xcluster)(ms,pn)=  nanmean(fp_locked_trials(:,:,ms,pn),2);
         end
     end
    %get the aligned data if it exists for the unit 
    aligned_trials.(xcluster)= fp_locked_trials;
    max_low_dist.(xcluster) = max_low_dist_unit;
    
    
%end


%%store all trials of a given peak, in a matrix, across all units
%clear aligned_trials
trials_dat = nan(251,6,4,41, length(fieldnames(aligned_trials)));%Dim1 =length data that we wanna look at (+-125ms for each msacc) on 4 peaks,Dim2: 41 = max number of trials found in a unit in the variable aligned_trials
clear i 
for i = 1:length(fieldnames(aligned_trials))
    %xcluster = xfilenames{31};
    %if ~isnan(max_low_dist.(xcluster))
        for tr = 1:length(aligned_trials.(xcluster)(1,:,1,1))
            for pn = 1:4
                for ms =1:6
                    if ~isempty(find(aligned_trials.(xcluster)(:,tr,ms,pn),1))
                        % aligned_trials(250*(pn-1)+1:250*pn+1,i)= mean(suas_trials(layer_idx(i)).aligned(max_low_dist(layer_idx(i))-125:max_low_dist(layer_idx(i))+125,:,pn),2);
                        trials_dat(1:251,ms,pn,tr,i)= aligned_trials.(xcluster)(max_low_dist.(xcluster)-1-124:max_low_dist.(xcluster)+125,tr,ms,pn);
                    end
                end
            end
        end
        %normalizing with max and min of each unit
        %norm_aligned_trials(:,i) = (aligned_trials(:,i) - min(aligned_trials(:,i)))/(max(aligned_trials(:,i))-min(aligned_trials(:,i)));
        
    %end
end

for i =1:42


msMean = squeeze(nanmean(trials_dat,2));
pnMean = squeeze(nanmean(msMean,2));
trMean = squeeze(nanmean(pnMean,2));

figure();
plot(trMean)
vline(125)
title(xfilenames{31},'Interpreter', 'none')
end

