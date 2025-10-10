function data = tACSChallenge_AnalyzeStaircase(filename, sj)


%% script adapted by Giuseppe Di Dona from tACSChallenge_ImportData function written by Florian Kasten

% The function imports data from logfiles of the staircase procedures and
% calculates the threshold brightness (base brightness + target brightness) 
% which can be then used in the experiment.

% This function will produce a picture of the staircase results including a
% text label with the threshold brightness calculated as the (target + base
% brightness) averaged along the last 10 trials. A .csv table will be also 
% created with the data from all trials. The last value of the column
% "ThresholdBrightness" is the value that will be used in the experiment as 
% thresholded brightness

% Note that the analysis will ignore trials in which the change in luminance results
% in no change as the target LED intensity is the same as base brightness



dataLines = [1, Inf];

%%
opts = delimitedTextImportOptions("NumVariables", 14);
% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Var1", "Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.SelectedVariableNames = ["Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
tmp_matrix = readtable(filename, opts);

%% Convert to output type
tmp_matrix = table2array(tmp_matrix);
tmp_matrix = tmp_matrix(2:end,:);

%% Extract Data
% Delta T
data.raw_dt = diff(tmp_matrix(:,1));

% handle reset of processor clock in teensy during recording
% check if any dt is negative (step backwards in time)
if sum(data.raw_dt<0) > 0

    % replace with the mean of all other dt
    data.raw_dt(data.raw_dt < 0) =  mean(data.raw_dt(data.raw_dt > 0),'omitnan');

end

data.dt = mean(data.raw_dt,'omitnan')/10^6; %uS to Seconds
% Sampling Rate
data.Fs = 1/data.dt;
% tACS
data.tACS = tmp_matrix(:,2);
% BNC (V)
data.BNC_In = tmp_matrix(:,4);
data.BNC_Out = tmp_matrix(:,5);
% Buttons
data.L_Button = tmp_matrix(:,6);
data.R_Button = tmp_matrix(:,7);
% LED
data.LEDs = tmp_matrix(:,8:13);


%% Get onsets of button presses

% Fuse R_Button and left Button

LR_Button = max([data.R_Button data.L_Button], [], 2); %takes the max resp over all L/R button so it keeps the value of the active button which is always superior to the other (one is 0 and the other is 1)

respOnsets = [ 0 ; diff(LR_Button)]; % codes onsets (+1) and offsets (-1), added a 0 in the beginning to correct indexing

respOnsets(respOnsets < 0) = 0; % removes offset going from -1 to 0, only onsets remain (+1)

%% this is when the subject pressed the button (in sample points)
RespLat = find(respOnsets>0); %finds the index of response onsets (+1)

%% leverage the inactive central LED to remove intensity offset from all LED channels

data.LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2)); %subtracts the value of the central led (1st column) transformed into a matrix from all leds

%% merge LED signals into one
LED = max(data.LEDs, [], 2); 

%takes the max led over all leds so it keeps the value of the active LED which is always superior to the other which are alwyas 0.2 (now zero for the previous step)

%% get onsets of target LED
LEDOnsets = [ 0; diff(LED)]; %same procedure used for coding the onsets of the responses
LEDOnsets(LEDOnsets < 0) = 0;


%% this is when LEDs were on (in sample points)
LEDLat=find(LEDOnsets>0); %find INDEX when LED were on

%% this is the brightness of the LEDS in sample points

LEDBright = LED(LEDLat); %this is the TARGET BRIGHTNESS which will need to be addedd to base brightness not the FINAL

BASEBright = 0.2;

% create matrix of behavioural responses
trialcounter = 0;

for i = 1:length(LEDLat)

    trialcounter = trialcounter+1;
    curr_t = LEDLat(i);
    trials_stair(trialcounter,1) = curr_t;

    if any(RespLat > curr_t & RespLat < curr_t + 1200)

        trials_stair(trialcounter,2) = 1; % hit
        RT = RespLat(RespLat > curr_t & RespLat < curr_t + 1200); %% target is considered detected if button press occurred within 1.2 s
        trials_stair(trialcounter,3) = min(RT) - curr_t; % response time

    end

    trials_stair(trialcounter,4) = LEDBright(trialcounter); % Target brightness (delta)
    trials_stair(trialcounter,5) = LEDBright(trialcounter) + BASEBright; % final brighness = Target + base
    trials_stair(trialcounter,6) = trialcounter; % trial number

  
end

%% Calculate the thresholded brightness based on the last 10 trials
% for each trial. The last values of column 7 is the thresholded brightness to be used

trials_stair(:,7) = zeros(size(trials_stair,1),1);


for i = 10:size(trials_stair,1)
    
    trials_stair(i,7) = mean(trials_stair(i-9 : i, 5));

end

threshold_bright = trials_stair(size(trials_stair,1),7); % threhsolded brightness

%% Save the table in .csv

tab = array2table(trials_stair,'VariableNames',{'Latency', 'Response', 'ISI','TargetBrightness_delta', 'FinalBrightness', 'TrialNumber','ThresholdBrightness'});
formatDescr = 'dd_mm_yyyy_HH_MM_SS';
datestring = datestr(now, formatDescr);
mkdir('StaircaseResults')
writetable(tab,['StaircaseResults/Staircase_SJ_' sj '_DATE_' datestring '.csv'])

%% Plot the result

markers = {'v', '^',};
colors = [1 0 0; 0 0.6 0]; %Red  Green
markers_idx = trials_stair(:,2) +1;  % This will determine the marker type (1 is X for miss and 2 is . for hit)


figure; 
plot(trials_stair(:,5))
legend;
xlabel('Trials');
ylabel('Final Brightness (Target + Base 0.2 )');
title(['Staircase Subject ' sj]);
hold on;

for i = 1:length(markers)
    idx = markers_idx == i;
    scatter( trials_stair(idx,6), trials_stair(idx,5), 30, 'filled', markers{i}, ...
        'MarkerEdgeColor',colors(i,:) ,'MarkerFaceColor', colors(i,:));
end

hold on;
plot([size(trials_stair,1)-9 size(trials_stair,1)],[threshold_bright threshold_bright], 'b-', 'LineWidth', 2);  % 'k-' is a black solid line
text(size(trials_stair,1)-15, threshold_bright+0.05, ['Threshold Brightness = ' num2str(threshold_bright)],'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');
legend({'','Miss','Hit',''})

% Save the result
saveas(gcf,['StaircaseResults/Staircase_SJ_' sj '_DATE_' datestring '.png'])


