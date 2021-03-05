%This script was developped to plot specific trials including microsaccades
%and see whether neural activity is modulated by microsaccades on specific
%trials


indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']); %neural data


%selected trials without spaces
condSelectedTrialsIdx = struct();

for n =1:length(selected_trials_idx.logicals)
    if ~isempty(selected_trials_idx.logicals(n).idx)
      condSelectedTrialsIdx.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = selected_trials_idx.logicals(n).idx;
    end 
end

xfilenames = fieldnames(condSelectedTrialsIdx);
cnt = 0;
eyeMovData = struct();

for i=27:length(fieldnames(condSelectedTrialsIdx))
    try
        xcluster = xfilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
        
        xBRdatafiles = concat_filenames.(xcluster);
        eye_info =[];
        all_codes = [];
        all_times = [];
        all_analogData =[];
        for fn =1:length(xBRdatafiles)
            xBRdatafile = xBRdatafiles{fn};
            filename   = [directory xBRdatafile(2:end)];
            if exist(strcat(filename, '.bhv'),'file')
                eye_info.(strcat(xBRdatafile,'_bhvfile')) = concatBHV(strcat(filename,'.bhv'));
                all_codes = [all_codes, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeNumbers];
                all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
            end
            
        end
        samplerate = 1000;
        
        %penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
        %STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
        trialindex = condSelectedTrialsIdx.(xcluster);
        
        ampl = [];
        veloc = [];
        ntr =0; %number of trials in which microsaccades were detected
        
        for tr = 1:length(trialindex)
            
            %if trialindex(tr) <= length(eye_info.CodeNumbers)
            codes                 = all_codes{trialindex(tr)};
            times                 = all_times{trialindex(tr)};
            
            
            %timestamps = zeros(length( eye_info.AnalogData{1,1}.EyeSignal(:,1)),1);
            %timestamps(trialindex) = trialindex;
            if nnz(find( codes == 23))
                samples = [];
                samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
                %(-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{1,1}.EyeSignal(:,1)) - times(codes == 23)); %timestamps of the recording in miliseconds
                %1:length(eye_info.AnalogData{1,1}.EyeSignal(:,1)); %timestamps of the recording in miliseconds
                if ~isempty(samples)
                    samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                    samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                    samples(:,4) = nan();
                    samples(:,5) = nan();
                    blinks = zeros(length(samples(:,1)),1);
                    recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                    
                    % Runs the saccade detection
                    [saccades, stats] = recording.FindSaccades();
                    
                    % Plots a main sequence
                    enum = ClusterDetection.SaccadeDetector.GetEnum;
                    ampl = [ampl; saccades(:,enum.amplitude)];
                    veloc = [veloc; saccades(:,enum.peakVelocity)];
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                    if ~all(isnan((saccades(:,enum.startIndex))))
                        ntr =ntr+1;
                    end
                    
                end
            end
        end
        eyeMovData.(xcluster).cellclass = trialsTraces.peak_aligned_trials.(xcluster).cellclass;
        %catch
        %     cnt = cnt+1;
        %     disp(xBRdatafile)
        % end
        
        
        figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        subplot(4,2,1)
        plot(ampl,veloc,'o')
        xlabel('Saccade amplitude (deg)');
        ylabel('Saccade peak velocity (deg/s)');
        set(gca,'box','off');
        set(gca, 'linewidth',1)
        xlim([0 2])
        title(strcat(sprintf('Microsaccades detected in %d trials',ntr)),'Interpreter', 'none')
        % Plots the traces with the labeled microsaccades
        
        if ~isempty(samples)
            subplot(4,2,[3:4])
            plot(samples(:,1), samples(:,2:3),'linewidth',1);
            hold
            yl = get(gca,'ylim');
            u1= zeros(size(samples(:,1)))+yl(1);
            u2= zeros(size(samples(:,1)))+yl(1);
            u1((saccades(:,enum.startIndex))) = yl(2);
            u2(saccades(:,enum.endIndex)) = yl(2);
            uOri = (cumsum(u1)-cumsum(u2));
            u =double(ischange(uOri));
            u(u ==0) = NaN;
            idc = find(~isnan(u));
            %u(idc-1) =0;
            u(idc+1) =0;
            plot(samples(:,1), u.*uOri./uOri,'k','linewidth',1)
            hold on
            uOri(uOri ==0) = NaN;
            plot(samples(:,1),uOri./uOri,'k','linewidth',1)
            xlim([-200 1500])
            ylim([-5 5])
            set(gca,'XTickLabel',[]);
            set(gca,'box','off');
            set(gca, 'linewidth',1)
            ylabel('Eye Position (deg)');
            legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
            
            
            subplot(4,2,[5:6])
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,2:3),'linewidth',1)
            hold
            yl = get(gca,'ylim');
            u3= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
            u4= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
            u3((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.startIndex))) = yl(2);
            u4(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.endIndex)) = yl(2);
            uOri = (cumsum(u3)-cumsum(u4));
            u =double(ischange(uOri));
            u(u ==0) = NaN;
            idc = find(~isnan(u));
            %u(idc-1) =0;
            u(idc+1) =0;
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
            hold on
            uOri(uOri ==0) = NaN;
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uOri./uOri,'k','linewidth',1)
            xlim([-200 1500])
            ylim([-5 5])
            set(gca,'box','off');
            set(gca,'XTickLabel',[]);
            set(gca, 'linewidth',1)
            ylabel('Eye Position (deg)');
            legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
            
            subplot(4,2,[7:8])
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,2:3),'linewidth',1)
            hold
            yl = get(gca,'ylim');
            u5= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
            u6= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
            u5((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.startIndex))) = yl(2);
            u6(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.endIndex)) = yl(2);
            uOri = (cumsum(u5)-cumsum(u6));
            uOri(uOri > 0)= 1;
            u =double(ischange(uOri*4));
            u(u ==0) = NaN;
            idc = find(~isnan(u));
            %u(idc-1) =0;
            u(idc+1) =0;
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
            hold on
            uOri(uOri ==0) = NaN;
            plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uOri./uOri,'k','linewidth',1)
            xlim([-200 1500])
            ylim([-5 5])
            set(gca,'box','off');
            set(gca, 'linewidth',1)
            xlabel('Time from stimulus onset (ms)');
            ylabel('Eye Position (deg)');
            
            legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
        end
        
        %saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.png'));
        %saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.svg'));
        
    catch
        cnt = cnt+1;
        disp(xBRdatafile)
    end
end


%% Simply plot a main sequence

figure('Renderer', 'painters', 'Position', [10 10 1200 1200]);
plot(ampl,veloc,'o','MarkerFaceColor', [146/255 197/255 222/255], 'MarkerEdgeColor', 'none')
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');
set(gca,'box','off');
set(gca, 'linewidth',1)
xlim([0 2])
title(strcat(sprintf('Main sequence in %d trials',ntr)),'Interpreter', 'none')
filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_',xcluster); 
saveas(gcf,strcat(filename,'.png'));
saveas(gcf,strcat(filename,'.svg'));

%% Plot a few eye movement trials
 figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
 % Plots the traces with the labeled microsaccades
 
 if ~isempty(samples)
     subplot(3,2,[1:2])
     plot(samples(:,1), samples(:,2:3),'linewidth',1);
     hold
     yl = get(gca,'ylim');
     u1= zeros(size(samples(:,1)))+yl(1);
     u2= zeros(size(samples(:,1)))+yl(1);
     u1((saccades(:,enum.startIndex))) = yl(2);
     u2(saccades(:,enum.endIndex)) = yl(2);
     uOri = (cumsum(u1)-cumsum(u2));
     u =double(ischange(uOri));
     u(u ==0) = NaN;
     idc = find(~isnan(u));
     %u(idc-1) =0;
     u(idc+1) =0;
     plot(samples(:,1), u.*uOri./uOri,'k','linewidth',1)
     hold on
     uOri(uOri ==0) = NaN;
     plot(samples(:,1),uOri./uOri,'k','linewidth',1)
     xlim([-200 1500])
     ylim([-5 5])
     set(gca,'XTickLabel',[]);
     set(gca,'box','off');
     set(gca, 'linewidth',1)
     ylabel('Eye Position (deg)');
     legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
     
     
     subplot(3,2,[3:4])
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,2:3),'linewidth',1)
     hold
     yl = get(gca,'ylim');
     u3= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
     u4= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
     u3((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.startIndex))) = yl(2);
     u4(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.endIndex)) = yl(2);
     uOri = (cumsum(u3)-cumsum(u4));
     u =double(ischange(uOri));
     u(u ==0) = NaN;
     idc = find(~isnan(u));
     %u(idc-1) =0;
     u(idc+1) =0;
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
     hold on
     uOri(uOri ==0) = NaN;
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uOri./uOri,'k','linewidth',1)
     xlim([-200 1500])
     ylim([-5 5])
     set(gca,'box','off');
     set(gca,'XTickLabel',[]);
     set(gca, 'linewidth',1)
     ylabel('Eye Position (deg)');
     legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
     
     subplot(3,2,[5:6])
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,2:3),'linewidth',1)
     hold
     yl = get(gca,'ylim');
     u5= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
     u6= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
     u5((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.startIndex))) = yl(2);
     u6(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.endIndex)) = yl(2);
     uOri = (cumsum(u5)-cumsum(u6));
     uOri(uOri > 0)= 1;
     u =double(ischange(uOri*4));
     u(u ==0) = NaN;
     idc = find(~isnan(u));
     %u(idc-1) =0;
     u(idc+1) =0;
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
     hold on
     uOri(uOri ==0) = NaN;
     plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uOri./uOri,'k','linewidth',1)
     xlim([-200 1500])
     ylim([-5 5])
     set(gca,'box','off');
     set(gca, 'linewidth',1)
     xlabel('Time from stimulus onset (ms)');
     ylabel('Eye Position (deg)');
     
     legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
 end
 filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eye_movement_examples_',xcluster); 
saveas(gcf,strcat(filename,'.png'));
saveas(gcf,strcat(filename,'.svg'));

%% Plot electro-oculograms together with neurophysiological data

%load neural data
monodatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
neurfilename = [monodatadir 'su_peaks_03032020_corrected\all_units\clean_origin_sup_50'];
neuralDat = load(strcat(neurfilename));

%load selected trial indices
indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices

%Store neural data and trial indices in new structure without empty space
condSelectedTrialsIdx = struct();
condNeuralDat =struct();
for n =1:length(neuralDat.clean_origin_data)
    if ~isempty(neuralDat.clean_origin_data(n).unit)
        condSelectedTrialsIdx.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = selected_trials_idx.logicals(n).idx;
        %use selected_trial_idx to get cluster filename for each single
        %unit so the neural data is identifiable
        condNeuralDat.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = neuralDat.clean_origin_data(n).unit;
    end
end

% eye movement data
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
trialsTraces =load([newdatadir 'all_orig_bs_zscore_trials']); %neural data

xfilenames = fieldnames(condSelectedTrialsIdx);
cnt = 0;
eyeMovData = struct();
for i=27:length(fieldnames(condSelectedTrialsIdx))
    
    try
        xcluster = xfilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
        
        xBRdatafiles = concat_filenames.(xcluster);
        eye_info =[];
        all_codes = [];
        all_times = [];
        all_analogData =[];
        for fn =1:length(xBRdatafiles)
            xBRdatafile = xBRdatafiles{fn};
            filename   = [directory xBRdatafile(2:end)];
            if exist(strcat(filename, '.bhv'),'file')
                eye_info.(strcat(xBRdatafile,'_bhvfile')) = concatBHV(strcat(filename,'.bhv'));
                all_codes = [all_codes, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeNumbers];
                all_times = [all_times, eye_info.(strcat(xBRdatafile,'_bhvfile')).CodeTimes];
                all_analogData = [all_analogData,eye_info.(strcat(xBRdatafile,'_bhvfile')).AnalogData];
                xBaseline = eye_info.(strcat(xBRdatafile,'_bhvfile')).ScreenXresolution/4/eye_info.(strcat(xBRdatafile,'_bhvfile')).PixelsPerDegree; %since the screen monitor is split in 2 parts with the stereoscope, the center for each eye becomes the center for each side of the stereoscope (half of the half, justifying dividing by 4
            end
            
        end
        samplerate = 1000;
        
        %penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
        %STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
        trialindex = condSelectedTrialsIdx.(xcluster);
        
        ampl = [];
        veloc = [];
        ntr =0; %number of trials in which microsaccades were detected
        
        for tr = 1:length(trialindex)
            
            %if trialindex(tr) <= length(eye_info.CodeNumbers)
            codes                 = all_codes{trialindex(tr)};
            times                 = all_times{trialindex(tr)};
            
            
            %timestamps = zeros(length( eye_info.AnalogData{1,1}.EyeSignal(:,1)),1);
            %timestamps(trialindex) = trialindex;
            if nnz(find( codes == 23))
                samples = [];
                samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
                %(-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{1,1}.EyeSignal(:,1)) - times(codes == 23)); %timestamps of the recording in miliseconds
                %1:length(eye_info.AnalogData{1,1}.EyeSignal(:,1)); %timestamps of the recording in miliseconds
                if ~isempty(samples)
                    samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1)+xBaseline; %horizontal position of the left eye in degrees baseline corrected
                    samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees
                    samples(:,4) = nan();
                    samples(:,5) = nan();
                    blinks = zeros(length(samples(:,1)),1);
                    recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);
                    
                    % Runs the saccade detection
                    [saccades, stats] = recording.FindSaccades();
                    
                    % Plots a main sequence
                    enum = ClusterDetection.SaccadeDetector.GetEnum;
                    ampl = [ampl; saccades(:,enum.amplitude)];
                    veloc = [veloc; saccades(:,enum.peakVelocity)];
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).saccades = saccades;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).enum = enum;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).stats = stats;
                    eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples = samples;
                    if ~all(isnan((saccades(:,enum.startIndex))))
                        ntr =ntr+1;
                    end
                    
                end
            end
        end
        eyeMovData.(xcluster).cellclass = trialsTraces.peak_aligned_trials.(xcluster).cellclass;
        
        
        %%Plot trials of both electro-oculograms and neural activity
        figure('Renderer', 'painters', 'Position', [10 10 1000 1200]);
        %%%Plots the traces with the labeled microsaccades %%%
        
        %if ~isempty(samples)
        subplot(6,2,[1:2])
        x = -600:1300;
        plot(x, condNeuralDat.(xcluster)(:,tr),'linewidth',1)
        xlim([-600 1300])
        ylabel('Spike rate (spikes/s)')
        set(gca, 'box','off');
        set(gca, 'linewidth',1)
        
        subplot(6,2,[3:4])
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples(520:end,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr))).samples(520:end,2:3),'linewidth',1);
        
        hold
        yl = get(gca,'ylim');
        u1= zeros(size(samples(:,1)))+yl(1);
        u2= zeros(size(samples(:,1)))+yl(1);
        u1((saccades(:,enum.startIndex))) = yl(2);
        u2(saccades(:,enum.endIndex)) = yl(2);
        uOri = (cumsum(u1)-cumsum(u2));
        u =double(ischange(uOri));
        u(u ==0) = NaN;
        idc = find(~isnan(u));
        %u(idc-1) =0;
        u(idc+1) =0;
        plot(samples(:,1), u.*uOri./uOri,'k','linewidth',1)
        hold on
        uOri(uOri ==0) = NaN;
        plot(samples(:,1),uOri./uOri,'k','linewidth',1)
        xlim([-600 1300])
        %ylim([-5 5])
        set(gca,'XTickLabel',[]);
        set(gca,'box','off');
        set(gca, 'linewidth',1)
        ylabel('Eye Position (deg)');
        %legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
        
        subplot(6,2,[5:6])
        x = -600:1300;
        plot(x, condNeuralDat.(xcluster)(:,tr-6),'linewidth',1)
        xlim([-600 1300])
        set(gca, 'box','off');
        set(gca, 'linewidth',1)
        
        subplot(6,2,[7:8])
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,2:3),'linewidth',1)
        hold
        yl = get(gca,'ylim');
        u3= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
        u4= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1)))+yl(1);
        u3((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.startIndex))) = yl(2);
        u4(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).enum.endIndex)) = yl(2);
        uOri = (cumsum(u3)-cumsum(u4));
        u =double(ischange(uOri));
        u(u ==0) = NaN;
        idc = find(~isnan(u));
        %u(idc-1) =0;
        u(idc+1) =0;
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
        hold on
        uOri(uOri ==0) = NaN;
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uOri./uOri,'k','linewidth',1)
        xlim([-600 1300])
        %ylim([-5 5])
        set(gca,'box','off');
        set(gca,'XTickLabel',[]);
        set(gca, 'linewidth',1)
        
        
        subplot(6,2,[9:10])
        x = -600:1300;
        plot(x, condNeuralDat.(xcluster)(:,tr-4),'linewidth',1)
        xlim([-600 1300])
        set(gca, 'box','off');
        set(gca, 'linewidth',1)
        
        subplot(6,2,[11:12])
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1),eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,2:3),'linewidth',1)
        hold
        yl = get(gca,'ylim');
        u5= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
        u6= zeros(size(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1)))+yl(1);
        u5((eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.startIndex))) = yl(2);
        u6(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).saccades(:,eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).enum.endIndex)) = yl(2);
        uOri = (cumsum(u5)-cumsum(u6));
        uOri(uOri > 0)= 1;
        u =double(ischange(uOri*4));
        u(u ==0) = NaN;
        idc = find(~isnan(u));
        %u(idc-1) =0;
        u(idc+1) =0;
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), u.*uOri./uOri,'k','linewidth',1)
        hold on
        uOri(uOri ==0) = NaN;
        plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uOri./uOri,'k','linewidth',1)
        xlim([-600 1300])
        %ylim([-5 5])
        set(gca,'box','off');
        set(gca, 'linewidth',1)
        xlabel('Time from stimulus onset (ms)');
        filename = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\electrooculogram_neurophys',xcluster);
        saveas(gcf,strcat(filename,'.png'));
        saveas(gcf,strcat(filename,'.svg'));
        %legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
    catch
        cnt = cnt+1;
        disp(xBRdatafile)
    end
end