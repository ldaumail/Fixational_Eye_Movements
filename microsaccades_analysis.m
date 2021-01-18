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


metadir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']);
metafilename = load(strcat(metadir, 'single_units_ns6_metadata.mat'));

for i = 1:length(metafilename.ns6FileName)
    if ~isempty(metafilename.ns6FileName{i,1})
        underscore = strfind(metafilename.ns6FileName{i,1}, '_');
        session = metafilename.ns6FileName{i,1};
        session =  session(1:underscore(2)-1);
        directory = strcat('C:\Users\daumail\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\bhv_selected_units\',session,'\');

        BRdatafile = metafilename.ns6FileName{i,1};
        filename   = [directory BRdatafile]; 
        if exist(strcat(filename, '.bhv'),'file') 
            eye_info = concatBHV(strcat(filename,'.bhv'));

            samplerate = 1000;



             penetration_name = erase(selected_trials_idx.logicals(i).penetration, 'matDE50_NDE0');
             STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
             trialindex = selected_trials_idx.logicals(i).idx;
             %cnt =0;
             ampl = [];
             veloc = [];

             for tr = 1:length(trialindex)
                 %try
                 if trialindex(tr) <= length(eye_info.CodeNumbers)
                     codes                 = eye_info.CodeNumbers{trialindex(tr)};
                     times                 = eye_info.CodeTimes{trialindex(tr)};


                     %timestamps = zeros(length( eye_info.AnalogData{1,1}.EyeSignal(:,1)),1);
                     %timestamps(trialindex) = trialindex;
                     if nnz(find( codes == 23))
                         samples = [];
                         samples(:,1) = (-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,1)) - times(codes == 23));
                         %(-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{1,1}.EyeSignal(:,1)) - times(codes == 23)); %timestamps of the recording in miliseconds
                         %1:length(eye_info.AnalogData{1,1}.EyeSignal(:,1)); %timestamps of the recording in miliseconds
                         if ~isempty(samples) 
                             samples(:,2) = eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,1); %horizontal position of the left eye in degrees
                             samples(:,3) = eye_info.AnalogData{trialindex(tr)}.EyeSignal(:,2); %vertical position of the left eye in degrees 
                             samples(:,4) = nan();
                             samples(:,5) = nan();
                             blinks = zeros(length(samples(:,1)),1);
                             recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, blinks, samplerate);

                            % Runs the saccade detection
                            [saccades stats] = recording.FindSaccades();
                   % catch
                   %  cnt = cnt+1;
                   % end


                            % Plots a main sequence
                            enum = ClusterDetection.SaccadeDetector.GetEnum;
                            ampl = [ampl; saccades(:,enum.amplitude)];
                            veloc = [veloc; saccades(:,enum.peakVelocity)];
                         end
                     end
                 end
             end

        figure();
        subplot(2,2,1)
        plot(ampl,veloc,'*')
        xlabel('Saccade amplitude (deg)');
        ylabel('Saccade peak velocity (deg/s)');
        xlim([0 2])
        title(session,'Interpreter', 'none')
        % Plots the traces with the labeled microsaccades
                if ~isempty(samples)
                    subplot(2,2,[3:4])
                    plot(samples(:,1), samples(:,2:end));
                    hold
                    yl = get(gca,'ylim');
                    u1= zeros(size(samples(:,1)))+yl(1);
                    u2= zeros(size(samples(:,1)))+yl(1);
                    u1((saccades(:,enum.startIndex))) = yl(2);
                    u2(saccades(:,enum.endIndex)) = yl(2);
                    u = cumsum(u1)-cumsum(u2);
                    plot(samples(:,1), u,'k')

                    xlabel('Time (ms)');
                    ylabel('Eye Position (deg)');

                    legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})
                end
        end
    end
end


figure
subplot(2,2,1)
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');

%{
% Plots the traces with the labeled microsaccades
subplot(2,2,[3:4])
plot(samples(:,1), samples(:,2:end));
hold
yl = get(gca,'ylim');
u1= zeros(size(samples(:,1)))+yl(1);
u2= zeros(size(samples(:,1)))+yl(1);
u1((saccades(:,enum.startIndex))) = yl(2);
u2(saccades(:,enum.endIndex)) = yl(2);
u = cumsum(u1)-cumsum(u2);
plot(samples(:,1), u,'k')

xlabel('Time (ms)');
ylabel('Eye Position (deg)');

legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})

end
%}