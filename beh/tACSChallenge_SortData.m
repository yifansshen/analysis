function trials_sorted = tACSChallenge_SortData(data_path, subj, conditions)
%% script originally written by Florian Kasten, University of Oldenburg
%% modified by Benedikt Zoefel, CNRS Toulouse, in October 2021, April 2022 (clean tACS signal), and June 2022 (correction for imperfect stimulation frequency)
%% modified by Florian Kasten in March 2022 (correction of target onsets)
%% modified by BZ in June 2025 (check for number of blocks)

%% data_path refers to the folder the data is located in (each subject in separate folder).
%% example for condition labels: conditions = {'*Montage A*','*Montage B*','*Montage C*'};
%% subj is subject number in string format (e.g., '01'). data must be located in a folder with the same name.

trials_sorted = cell(length(conditions),1);
phase_bins = 0:pi/4:2*pi;

for c = 1:length(conditions)
    
    curr_files = dir([data_path, subj, filesep, conditions{c}]);
    trialcounter = 0;
    trials_sorted{c} = [];
    
    if length(curr_files) ~= 3 && strcmp('*_Sh_*',conditions{c})==0
        error('number of blocks for at least one condition is not three as expected');
    end
    
    for b = 1:length(curr_files)
        
        %% import data
        data = tACSChallenge_ImportData([curr_files(b).folder, filesep, curr_files(b).name]);
        
        %% get onsets of button presses
        respOnsets = diff(data.R_Button);
        respOnsets(respOnsets < 0) = 0;
        %% this is when the subject pressed the button (in sample points)
        RespLat = find(respOnsets>0);
        
        %% leverage the inactive central LED to remove intensity offset from all LED channels
        data.LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2));
        %% merge LED signals into one
        LED = max(data.LEDs, [], 2);
        %% get onsets of target LED
        LED = diff(LED);
        LED(LED < 0) = 0;
        %% this is when LEDs were on (in sample points)
        LEDLat=find(LED>0);
        
        %% tACS signal
        tACS = data.tACS; tACS = tACS-mean(tACS);
        % clean the signal - we will fit a sine wave, starting a few
        % seconds before the first target
        tACS_to_fit = tACS(LEDLat(1)-1000:end);
        
        % sometimes the tACS signal is not precisely 10 Hz. we therefore
        % estimate its frequency
        % needs to be in the range 9-11 Hz
        [filt_b,filt_a] = butter(2,[9/(data.Fs/2) 11/(data.Fs/2)],'bandpass');
        tACS_filtered = filtfilt(filt_b,filt_a,tACS_to_fit);
        tACS_filtered_h = hilbert(tACS_filtered);
        % instantaneous frequency
        ifq = instfreq(tACS_filtered,data.Fs);
              
        % now reconstruct tACS sine
        t_fit = 0:1/round(data.Fs):length(tACS_to_fit)/round(data.Fs)-1/round(data.Fs);
        % make a sine with the same average frequency and the same phase
        % (skip the first and last few values to avoid edge effects)
        tACS_sine = cos(2*pi*t_fit*mean(ifq(4:end-3))+angle(tACS_filtered_h(round(data.Fs)*10+1)));       
        tACS_phase = angle(hilbert(tACS_sine));
        
        %% go through all trials
        for i = 1:length(LEDLat) 
            
            trialcounter = trialcounter+1;
            curr_t = LEDLat(i);
            trials_sorted{c}(trialcounter,1) = curr_t;
            
            %% target is considered detected if button press occurred within 2 s
            if any(RespLat > curr_t & RespLat < curr_t + 2000)
                
                trials_sorted{c}(trialcounter,2) = 1; % hit
                RT = RespLat(RespLat > curr_t & RespLat < curr_t + 2000);
                trials_sorted{c}(trialcounter,3) = min(RT) - curr_t; % response time
            end
                t_tACS = curr_t-LEDLat(1)+1000; % this is because the tACS signal starts later (see signal cleaning above)
                trials_sorted{c}(trialcounter,4) = tACS_phase(t_tACS); % tACS phase at the time of target presentation
                [~,trials_sorted{c}(trialcounter,5)] = min(abs(tACS_phase(t_tACS)+pi-phase_bins)); % number of nearest phase bin (for visualisation in later scripts)
            
        end
    end
end


