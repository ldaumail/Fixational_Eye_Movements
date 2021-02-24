eyeMovDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\data\';
eyeMovDat = load( [eyeMovDir, 'all_eye_movement_data']); %eye movement data
indexdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'stim_selected_trials_idx']); %trial indices


monodatadir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
neurfilename = [monodatadir 'su_peaks_03032020_corrected\all_units\clean_origin_sup_50'];
neuralDat = load(strcat(neurfilename));
condSelectedTrialsIdx = struct();
for n =1:length(neuralDat.clean_origin_data)
    if ~isempty(neuralDat.clean_origin_data(n).unit)
        condSelectedTrialsIdx.(strcat('x',erase(selected_trials_idx.logicals(n).penetration, 'matDE50_NDE0'))) = selected_trials_idx.logicals(n).idx;
    end
end

amps = [];
xfilenames = fieldnames(eyeMovDat);

for i =24:length(fieldnames(eyeMovDat))
%for i =1:length(fieldnames(eyeMovDat))
    xcluster = xfilenames{i};
    trialindex = condSelectedTrialsIdx.(xcluster);
    for tr = 1:length(trialindex)
        if isfield(eyeMovDat.(xcluster),sprintf('t%d',trialindex(tr)))
            
            saccades = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).saccades;
            enum = eyeMovDat.(xcluster).(sprintf('t%d',trialindex(tr))).enum;
            amplitudes = saccades(:,enum.leftAmplitude);
            %amps = [amps; amplitudes];
            rightAmps = amplitudes(amplitudes <= 10);
            amps = [amps; rightAmps];
            
        end
    end
end

figure();
histogram(amps, 'BinWidth', 0.01, 'DisplayStyle', 'stairs')
xlim([0 5])
xlabel('Saccade amplitude (deg)')
ylabel('Number of microsaccades')
title('Microsaccade amplitude distribution of all eye movements data')

saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eyeMovements_hist_secondrun.png'));


figure();
%[N,edges] = histcounts(amps, 'Normalization','pdf');
[N,edges] = histcounts(amps);
edges = edges(2:end) - (edges(2)-edges(1))/2;
semilogx(edges, N);
xlim([0.01 10])
New_XTickLabel = get(gca,'xtick');
set(gca,'XTickLabel',New_XTickLabel);
set(gca, 'Box', 'off')
xlabel('Saccade amplitude (deg)')
ylabel('Number of microsaccades')
title('Microsaccade amplitude distribution of 2018 and over eye movements data')
saveas(gcf,strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data\single_units\microsaccades_adaptation_analysis\plots\eyeMovements_distribution_2018over.png'));


