%ns6_filename = '160602_I_cinterocdrft012';
ns6dir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\160602_I\160602_I_cinterocdrft012.ns6';
bhvdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\bhv_selected_files\160602_I_cinterocdrft012.bhv';

%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'selected_trials_idx']);

        penetration_name = erase(selected_trials_idx.logicals(1).penetration, 'matDE50_NDE0');
        STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
    
     %1)get the stim onset times for all trials of the given penetration
        STIM_onsets = STIM_file.STIM.photo_on;
     %only keep the selected trials onsets
        selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(1).idx));
        selected_STIM_onsets = selected_STIM_onsets(1:4:end);
     %contact index/letter
        contact = STIM_file.STIM.chan;
        
     %2) get appropriate ns6 file of penetration 
     %isolate session date/animal for ns6 folder name
        underscore = strfind(penetration_name, '_');
        session =  penetration_name(1:underscore(2)-1);

        %{
          ext          = 'ns6'; 

      electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext),electrode,'read','uV');
        DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!
%}
%convert data into DVA
[eye_dva,bhv,eye_ana] = getAnalogEyeinDVA(ns6dir,bhvdir);


%need to scale up the triggering onset and offset times
%for data sampled af FS = 15kHz, ==> 15 X more than FS
%at 1000 Hz
pre   = -500;
post  = 1500; 

trigX  = trigData(eye_ana(:,1),floor(selected_STIM_onsets./30),-pre,post); 
trigY  = trigData(eye_ana(:,2),floor(selected_STIM_onsets./30),-pre,post);

x = pre:post;

figure();
subplot(2,1,1)
plot(x,squeeze(trigX))
ylabel('X position')
xlabel('time from stimulus onset (ms)')

subplot(2,1,2)
plot(x,squeeze(trigY))
ylabel('Y position')
xlabel('time from stimulus onset (ms)')
sgtitle('Eye positions')

saveas(gcf,'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eye_positions_single_unit1.png')


%compute distance from center

d = sqrt(trigX.^2 +trigY.^2);

figure();
plot(x,squeeze(d))
ylabel('Position from center')
xlabel('time from stimulus onset (ms)')
sgtitle('Eye positions from center')
saveas(gcf,'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eye_positions_fromcenter_single_unit1.png')

%compute amplitude and velocity
%first triggered data
trigX_dva  = trigData(eye_dva(:,1),floor(selected_STIM_onsets./30),-pre,post); 
trigY_dva  = trigData(eye_dva(:,2),floor(selected_STIM_onsets./30),-pre,post);

amplitude = sqrt(trigX_dva(2:end).^2 +trigY_dva(2:end).^2)-sqrt(trigX_dva(1:end-1).^2 +trigY_dva(1:end-1).^2);

velocity = amplitude/(10^(-3));

figure();
plot(amplitude, velocity);
xlabel('Amplitude (ang deg)')
ylabel('Angular Velocity (ang deg/sec')
title('Velocity vs Amplitude')
saveas(gcf,'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eye_velocity_vs_amplitude_single_unit1.png')


%second: all the eye movements data


%% Using Otero-Millan 2014 toolbox to detect microsaccades



%
            %{
     %1)get the stim onset times for all trials of the given penetration
        STIM_onsets = STIM_file.STIM.photo_on;
     %only keep the selected trials onsets
        selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(1).idx));
        selected_STIM_onsets = selected_STIM_onsets(1:4:end);
        
     %contact index/letter
        contact = STIM_file.STIM.chan;
        
     %2) get appropriate ns6 file of penetration 
     %isolate session date/animal for ns6 folder name
        underscore = strfind(penetration_name, '_');
        session =  penetration_name(1:underscore(2)-1);
%}


%metadir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
%metafilename = load(strcat(metadir, 'single_units_ns6_metadata.mat'));

indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
newdatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\binocular_adaptation\all_units\';
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
for i = 26:length(fieldnames(condSelectedTrialsIdx))
    
        xcluster = xfilenames{i};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');

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
            end
            
        end
            samplerate = 1000;

             %penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
             %STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
             trialindex = condSelectedTrialsIdx.(xcluster);
             
             ampl = [];
             veloc = [];
             ntr =0; %number of trials in which microsaccades were detected
            try
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
                             samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1); %horizontal position of the left eye in degrees
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
            catch
                cnt = cnt+1;
                disp(xBRdatafile)
            end

       %{   
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
                    u = (cumsum(u1)-cumsum(u2))/2;
                    plot(samples(:,1), u,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
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
                    uu = (cumsum(u3)-cumsum(u4))/2;
                    plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uu,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
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
                    uuu = (cumsum(u5)-cumsum(u6))/2;
                    plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uuu,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
                    set(gca,'box','off');
                    set(gca, 'linewidth',1)
                    xlabel('Time from stimulus onset (ms)');
                    ylabel('Eye Position (deg)');

                    legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
                end
                saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.png'));
                saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.svg'));

            
%}
       
           
end
savedir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\data\';
%save( [savedir, 'all_eye_movement_data', '.mat'],'-struct', 'eyeMovData'); %save eye movement data 
    


%% Look up power spectrum density for 2016 files 

%metadir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
%metafilename = load(strcat(metadir, 'single_units_ns6_metadata.mat'));
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']);
concat_filenames = load( [indexdir, 'concat_filenames_completenames']); %cluster filenames
condSelectedTrialsIdx = struct();

for n =1:length(selected_trials_idx.logicals)
    if ~isempty(selected_trials_idx.logicals(n).idx)
      condSelectedTrialsIdx.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = selected_trials_idx.logicals(n).idx;
    end 
end

xfilenames = fieldnames(condSelectedTrialsIdx);

cnt = 0;
Ses = [];



for i = 1:length(fieldnames(condSelectedTrialsIdx))
    
        xcluster = xfilenames{i};
        
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        
        directory = strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');
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
            end
            
        end
        %BRdatafile = concat_filenames.(xcluster);
        %filename   = strcat(directory, BRdatafile(2:end)); 
        %if exist(strcat(filename, '.bhv'),'file') 
            %eye_info = concatBHV(strcat(filename,'.bhv'));
            
            samplerate = 1000;



             %penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
             %STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
             %trialindex = selected_trials_idx.logicals(i).idx;
             %cnt =0;
             ampl = [];
             veloc = [];
            trialindex = condSelectedTrialsIdx.(xcluster);
       
            %try
             for tr = 1:length(trialindex)
                 
                 if trialindex(tr) <= length(eye_info.CodeNumbers)
                    % codes                 = eye_info.(xcluster).CodeNumbers{trialindex(tr)};
                     %times                 = eye_info.(xcluster).CodeTimes{trialindex(tr)};
                     codes                 = all_codes{trialindex(tr)};
                     times                 = all_times{trialindex(tr)};
                     

                      if nnz(find( codes == 23))
                         samples = [];
                         %samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
                         samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(all_analogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
                        
                         if ~isempty(samples) 
                             %{
                             samples(:,2) = eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,1); %horizontal position of the left eye in degrees
                             samples(:,3) = eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees 
                            samples(:,4) = nan();
                             samples(:,5) = nan();
                             
                            %}
                             samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1); %horizontal position of the left eye in degrees
                             samples(:,3) = all_analogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees 
                             samples(:,4) = nan();
                             samples(:,5) = nan();
                             blinks = zeros(length(samples(:,1)),1);
                             recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);

                            % Runs the saccade detection
                            [saccades stats] = recording.FindSaccades();


                            % Plots a main sequence
                            enum = ClusterDetection.SaccadeDetector.GetEnum;
                            ampl = [ampl; saccades(:,enum.amplitude)];
                            veloc = [veloc; saccades(:,enum.peakVelocity)];
                            
                            fs     = abs(fft(samples(:,2))/length(samples(:,2))*2); 
                            %strTrial = sprintf('Trial%d',tr);
                            %filenb = sprintf('file%d', i);
                            Ses(tr,1:200,i) = fs(1:200);
                            

                                                 
                         end
                     end
                 end
                
             end
             
             
         figure();
         subplot(1,2,1)
        plot(ampl,veloc,'*')
        xlabel('Saccade amplitude (deg)');
        ylabel('Saccade peak velocity (deg/s)');
        xlim([0 2])
         title(session,'Interpreter', 'none') 
        subplot(1,2,2)
         plot(mean(Ses(:,:,i),1));
         %ylim([2 20]);

         xlabel('Amplitude')
         xlabel('Frequency band (Hz)')
                
             %}
           % catch
           %     cnt = cnt+1;
           %     disp(BRdatafile)
           % end

            
     

        end
    
end

%% compare file names from code from 12-15-2020 (get_good_singleunits_ns6_filenames.m) with code from 1-22-2021 (get_concat_filenames.m)
metadir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
oldfilenames = load(strcat(metadir, 'single_units_ns6_metadata.mat'));

ns6_filenames = oldfilenames.ns6FileName;
fullCells = find(~cellfun('isempty', ns6_filenames));
fullFilenames = ns6_filenames(fullCells);

for i =1:length(fullFilenames)
    filename = fullFilenames{i};
    sessions{i} = filename(1:8);
end
uniqueStrCell(sessions)

indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
newfilenames = load(strcat(indexdir,'concat_filenames.mat'));
fieldnames(newfilenames)

%% Analyze eye movement data with neural data
%{
monodatadir = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
neurfilename = [monodatadir 'su_peaks_03032020_corrected\all_units\clean_origin_sup_50'];
neuralDat = load(strcat(neurfilename));
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
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
eyeMovDir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\data\';
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
    
    for tr = 1:length(trialindex)
        clear trSaccLocs trPeakLocs
        if isfield(eyeMovDat.(xcluster),sprintf('t%d',trialindex(tr)))
            %get eye movement location (specify eye movement below 2 degrees..?)
            saccades = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).saccades;
            enum = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).enum;
            trSaccLocs = saccades(:,enum.startIndex)+200;
            %select eye movement locations within the peak spans (lets start
            %with +-125)
            trPeakLocs = condPeakLocs.(xcluster)(:,tr);%-200; %need to subtract 200 to peak location (or add 200 to msacc location) as the peak locations were estimated from trials that start 200 ms before stim onset (in get_clean_peaks_and_data.m)
            for pl =1:length(trPeakLocs)
                cnt =0;
                for sl=1:length(trSaccLocs)
                    if (trPeakLocs(pl)-125 <= trSaccLocs(sl)) && (trSaccLocs(sl)<= trPeakLocs(pl)+125)
                     cnt =cnt+1;   
                        selectSaccLocs(cnt,pl,tr) = trSaccLocs(sl); %dim1 =msacc location, dim2 = peak number, dim3=trial
                    end
                end
            end
            
        end
        
    end
    allSelectSaccLocs.(xcluster) = selectSaccLocs;
end


%2) Align neural data on microsaccade onset from  peak1 to peak 4
aligned_trials = struct();
for i =1:length(fieldnames(eyeMovDat))
    xcluster = xfilenames{i};
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
    
    
end


%%store all trials of a given peak, in a matrix, across all units
%clear aligned_trials
trials_dat = nan(251,6,41, length(fieldnames(aligned_trials)));%Dim1 =length data that we wanna look at (+-125ms for each msacc) on 4 peaks,Dim2: 41 = max number of trials found in a unit in the variable aligned_trials
clear i 
for i = 1:length(fieldnames(aligned_trials))
    xcluster = xfilenames{i};
    %if ~isnan(max_low_dist.(xcluster))
        for tr = 1:length(aligned_trials.(xcluster)(1,:,1,1))
            for pn = 1:4
                for ms =1:6
                    if ~all(isnan(aligned_trials.(xcluster)(:,tr,ms,pn)))
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
plot(trMean(:,i))
vline(125)
title(xfilenames{i},'Interpreter', 'none')
end

%compute change of firing rate 
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
%}

%% Example plot

xcluster = xfilenames{31};
        cluster = xcluster(2:end);
        underscore = strfind(cluster, '_');
        session =  cluster(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\concat2_bhv_selected_units\',cluster,'\');

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
            end
            
        end
            samplerate = 1000;

             %penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
             %STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
             trialindex = condSelectedTrialsIdx.(xcluster);
             
             ampl = [];
             veloc = [];
             ntr =0; %number of trials in which microsaccades were detected
            try
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
                             samples(:,2) = all_analogData{trialindex(tr)}.EyeSignal(:,1); %horizontal position of the left eye in degrees
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
            catch
                cnt = cnt+1;
                disp(xBRdatafile)
            end

          
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
                    u = (cumsum(u1)-cumsum(u2))/2;
                    plot(samples(:,1), u,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
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
                    uu = (cumsum(u3)-cumsum(u4))/2;
                    plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-6))).samples(:,1), uu,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
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
                    uuu = (cumsum(u5)-cumsum(u6))/2;
                    plot(eyeMovData.(xcluster).(sprintf('t%d',trialindex(tr-4))).samples(:,1), uuu,'k','linewidth',1)
                    xlim([-200 1500])
                    ylim([-20 20])
                    set(gca,'box','off');
                    set(gca, 'linewidth',1)
                    xlabel('Time from stimulus onset (ms)');
                    ylabel('Eye Position (deg)');

                    legend({'Left Horiz', 'Left Vert', 'Microsaccades'},'Location','bestoutside')
                end
                saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.png'));
                saveas(gcf,strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\main_sequence_eyeMovExamples',xcluster,'.svg'));

            