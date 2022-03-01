% To avoid including a large data file in the repository, the data
% loaded here have been pre-binned using 10ms windows, i.e., each entry in
% the spikes matrices indicates the number of spikes recorded in a 10ms
% window
load mat_sample/sample_data.mat

%% Example computation of a cross-correlation map
% Figs. 3 and 4

% The units of the arguments are with respect to the binning window used
% to bin spikes.

argIn.BinWidth = 1;     % 10ms
argIn.MaxDelay = 10;    % 100ms
argIn.TimeStep = 4;     % 40ms
argIn.WindowLength = 8; % 80ms

%argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

disp(argIn)
argOut = ComputeCorrMap(spikes, expCond, argIn);

%%
CANONICAL_PAIR_IDX = 1;
mapDim = size(argOut.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay)*10; % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim)*10; % Convert to ms

figure(1);

imagesc( delays, t, argOut.CorrMap(:,:,CANONICAL_PAIR_IDX)' )

ax = gca;
ax.YDir = 'Normal';

xlabel('Delay')
ylabel('Time')

%% Example computation of the interaction structure analysis
% Fig. 6

clear argIn

argIn.TimePeriods = [...
    (  0:20:40)' ( 20:20:60)'; ...
    (128:20:168)' (148:20:188)'] + 5;

% Can take up to 15min due to the 10-fold cross-validation
argOut = CovStabilityAcrossTimeAnalysis(spikes, expCond, argIn);

%%
RANK_TO_PLOT = 2;

normFactor = diag(argOut.CvR(:,:,RANK_TO_PLOT));
numTimePeriods = size(argIn.TimePeriods, 1);

figure(2);

imagesc(argOut.CvR(:,:,RANK_TO_PLOT)./repmat(normFactor', numTimePeriods, 1))

ax = gca;
ax.YDir = 'Normal';

axis square

xlabel('Time Used For Correlation')
ylabel('Time Used For Fitting')

