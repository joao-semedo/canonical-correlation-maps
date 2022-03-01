function argOut = ComputeCorrMap(spikes, expCond, argIn)
% 
% argOut = ComputeCorrMap(spikes, expCond, argIn) computes a
% cross-correlation map using Canonical Correlation Analysis (CCA).
% 
%   p1: number of neurons in population 1
%   p2: number of neurons in population 2
%   T:  number of time points in a trial
%   N:  number of trials
% 
% INPUTS:
% 
% spikes  - a two element cell. spikes{1} is a p1 x T x N array containing
% the spiking activity in neuronal population 1. spikes{2} is a p2 x T x N
% array containing the spiking activity in neuronal population 2.
% 
% expCond - vector containing the trial label for each trial (N x 1)
% 
% argIn   - structure containing additional input arguments.
% 	argIn =
%       BinWidth:      Size of the binning window. Units are derived from
%       spikes, i.e., if each time point in spikes{1} corresponds to 1 ms
%       of activity, BinWidth will be in ms.
% 
%       MaxDelay:      Max delay considered between the windows of activity
%       in each area. Same units as BinWidth.
% 
%       TimeStep:      Time step with which activity windows progress
%       through the trial. Same units as BinWidth.
% 
%       WindowLength:  Length of the activity windows. Same units as
%       BinWidth.
% 
%		AreaIdx:       (Optional) [] (default), 1 or 2. Used to compute
%       within-area correlation maps. Activating this option deactivates
%       the ShuffleType argument.
% 
%       CrossValidate: (Optional) False (default) or True. Used to compute
%       a cross-validated cross-correlation map. When using
%       cross-validation, only the map associated with the first pair is
%       computed.
% 
%       NumWorkers:    (Optional) Number of workers used for parallel
%       computation. 0 (default) runs computations serially. Any number
%       larger than 0 will engage parallel computation. Setting
%       NumWorkers = Inf will use as many workers as the number of
%       physical CPU cores present in the machine. This option can be
%       *very* memory intensive.
% 
%       NumDraws:      (Optional) Number of random splits when computing
%       within-area correlation maps. Defaults to 1.
% 
%       NumShuffles:   (Optional) Number of shuffles. Defaults to 0 (i.e.,
%       by default, no shuffle control is performed).
% 
%       ShuffleType:   (Optional) Type of shuffle control. 1 - Shuffle
%       across-area correspondence across time and trials. 2 - Shuffle
%       across-area correspondence across time but not trials. 3 - Shuffle
%       across-area correspondence across trials but not time. 4 - Shuffle
%       across-area correspondence across trials but not time, and only
%       within each trial condition. 5 - Jitter shuffle.
% 
%       TimePoints:    (Optional) 1D array. Only compute cross-correlation
%       map at specific time points in the trial. By default, the
%       cross-correlation map is computed for all timepoints, i.e.,
%       t = 1:argIn.TimeStep_(T - WindowLength-1). When TimePoins is set,
%       t = argIn.TimePoints. Same units as BinWidth.
% 
% OUTPUTS:
% 
% argOut - structure containing outputs
% 	argOut =
%       CorrMap:         a 3D array. First dimension is time for the
%       activity window in area 1, which takes values
%       t1 = 1:argIn.TimeStep_(T - WindowLength-1). Second dimension is
%       relative time for the activity window in area 2, which takes values
%       t1-argIn.MaxDelay:t1+argIn.MaxDelay. Third dimension corresponds to
%       the number of canonical pairs used, where the first index contains
%       the map associated with the first canonical pair, etc. Each entry
%       contains the canonical correlation between the two activity
%       windows.
% 
%       FrMap:           a 3D array with the same dimensions and structure
%       as CorrMap. Each entry contains the geometric mean of the average
%       firing rate in each population.
% 
%       ShuffledCorrMap: a 3D array with the same dimensions and structure
%       as CorrMap, except for the third dimension, which corresponds to
%       each shuffle.
% 
%       UnitsRmvd:       a two element cell. UnitsRmvd{1} contains the
%       indexes of the units in population 1 that were excluded from the
%       analysis. Similarly for UnitsRmvd{2}.

argOut = struct(           ...
            'CorrMap', [], ...
              'FrMap', [], ...
    'ShuffledCorrMap', [], ...
          'UnitsRmvd', []  ...
    );

C_CV_NUM_FOLDS = 10;
C_MIN_FIRING_RATE = 0.5;

C_INTERNAL_BIN_WIDTH = 1; % This is for debugging purposes only, change at your own risk!
argIn.MaxDelay = argIn.MaxDelay/C_INTERNAL_BIN_WIDTH;
argIn.TimeStep = argIn.TimeStep/C_INTERNAL_BIN_WIDTH;
argIn.WindowLength = argIn.WindowLength/C_INTERNAL_BIN_WIDTH;

pop1 = 1;
pop2 = 2;

if ~isfield(argIn, 'AreaIdx')
    argIn.AreaIdx = [];
elseif (argIn.AreaIdx == pop1) || (argIn.AreaIdx == pop2)
    argIn.NumShuffles = 0;
end
if ~isfield(argIn, 'CrossValidate')
    argIn.CrossValidate = false;
end
if ~isfield(argIn, 'NumWorkers')
    argIn.NumWorkers = 0;
end
if ~isfield(argIn, 'NumDraws')
    argIn.NumDraws = 1;
end
if ~isfield(argIn, 'NumShuffles')
    argIn.NumShuffles = 0;
end
if ~isfield(argIn, 'ShuffleType')
    argIn.ShuffleType = 1;
end

if argIn.CrossValidate
    C_MAX_NUM_CANON_PAIRS = 1;
else
    C_MAX_NUM_CANON_PAIRS = 2;
end

[~, trialLength, ~] = size(spikes{pop1});

mapDim = trialLength - (argIn.WindowLength-1);

if isfield(argIn, 'TimePoints') && ~isempty(argIn.TimePoints)
    timePoints = argIn.TimePoints;
else
    argIn.TimePoints = [];
    timePoints = 1:argIn.TimeStep:mapDim;
end
numTimePoints = numel(timePoints);

[spikes, argOut.UnitsRmvd] ...
    = RmvSilentUnits(spikes, timePoints, ...
        argIn.BinWidth, argIn.MaxDelay, argIn.WindowLength);

firingRates = ComputeFiringRates(spikes, C_INTERNAL_BIN_WIDTH);
spikes{pop1}(firingRates{pop1} < C_MIN_FIRING_RATE,:,:) = [];
spikes{pop2}(firingRates{pop2} < C_MIN_FIRING_RATE,:,:) = [];

[numUnits(pop1), ~, numTrials] = size(spikes{pop1});
numUnits(pop2) = size(spikes{pop2}, 1);

if ~isempty(argIn.AreaIdx)
    if argIn.AreaIdx == pop1
        
        spikes{pop2} = spikes{pop1};
        
    elseif argIn.AreaIdx == pop2
        
        spikes{pop1} = spikes{pop2};
        
    end
end

numUnits(pop1) = size(spikes{pop1}, 1);
numUnits(pop2) = size(spikes{pop2}, 1);

[~, resids] = ComputePsths( ZScoreSpikes(spikes, expCond), expCond );

numDelays = 2*argIn.MaxDelay + 1;
N = argIn.WindowLength/(argIn.BinWidth/C_INTERNAL_BIN_WIDTH)*numTrials;
Y = zeros(          ...
    N,              ...
    numUnits(pop2), ...
    numDelays );

numCanonPairs = min(...
    [C_MAX_NUM_CANON_PAIRS, numUnits(pop1), numUnits(pop2)]);

argOut.CorrMap = zeros(     ...
    numDelays,              ...
    ceil(numTimePoints),    ...
    numCanonPairs );

argOut.FrMap = zeros(       ...
    numDelays,              ...
    ceil(numTimePoints) );

argOut.LowRankWarning = zeros(  ...
    numDelays,                  ...
    ceil(numTimePoints) );

argOut.ShuffledCorrMap = zeros( ...
    numDelays,                  ...
    ceil(numTimePoints),        ...
    argIn.NumShuffles );

iterTime = zeros(numTimePoints, 1);
msg = [];

numWorkers = StartParPool( min(numDelays, argIn.NumWorkers) );
for i = 1:numTimePoints
    
    tic
    
    t1 = timePoints(i);
    
    frX = mean( mean( Squash( BinTime(                          ...
        spikes{pop1}(:,t1:t1 + (argIn.WindowLength-1),:)		...
        , argIn.BinWidth ) ) ) );
    
    auxSpikeRastersX = spikes{pop1}(:,t1:t1 + (argIn.WindowLength-1),:);
    
    auxResidsX = BinTime(										...
        resids{pop1}(:,t1:t1 + (argIn.WindowLength-1),:)		...
        , argIn.BinWidth );
    T = size(auxResidsX, 2);
    
    X = Squash( auxResidsX )';
    
    t2 = max(t1 - argIn.MaxDelay, 1):min(t1 + argIn.MaxDelay, mapDim);
    
    % Prepare Y
    fr = zeros( numel(t2), 1 );
    for j = 1:numel(t2)
        t = t2(j);
        auxResidsY = BinTime( ...
            resids{pop2}(:,t:t + (argIn.WindowLength-1),:), ...
            argIn.BinWidth );
        Y(:,:,j) = Squash( auxResidsY )';
        
        frY = mean( mean( Squash( BinTime(               ...
            spikes{pop2}(:,t:t + (argIn.WindowLength-1),:)	...
            , argIn.BinWidth ) ) ) );
        
        fr(j) = geomean([frX frY]);
    end
    
    r = zeros( numel(t2), numCanonPairs );
    warnings = zeros( numel(t2), 1 );
    r_s = zeros( numel(t2), argIn.NumShuffles );
    parfor (j = 1:numel(t2), numWorkers)
        
        Yj = Y(:,:,j);
        
        if isempty(argIn.AreaIdx) %#ok<PFBNS>
            if argIn.CrossValidate
                c = cvpartition(N, 'kFold', C_CV_NUM_FOLDS);
                numPairs = min( size(X, 2), size(Yj, 2) );
                cvr = zeros(10, numPairs);
                for foldIdx = 1:C_CV_NUM_FOLDS
                    cvr(foldIdx,:) = CanonCorrFitAndPredict( ...
                        X(c.training(foldIdx),:), ...
                        Yj(c.training(foldIdx),:), ...
                        X(c.test(foldIdx),:), ...
                        Yj(c.test(foldIdx),:) ...
                        );
                end
                if find(isnan(cvr))
                    warnings(j) = true;
                end
                
                cvr = nanmean(cvr);
                
                r(j,:) = cvr(1:numCanonPairs);
            else
                [~, ~, aux, warnings(j)] = CanonCorr(X, Yj);
                
                r(j,:) = aux(1:numCanonPairs);
            end
            
            r_s_j = zeros(1, argIn.NumShuffles);
            switch argIn.ShuffleType
                
                case 1
                    for k = 1:argIn.NumShuffles
                        [~, ~, r_s_j_k] ...
                            = CanonCorr( X(randperm(N),:), Yj );
                        r_s_j(k) = r_s_j_k(1);
                    end
                    
                case 2
                    for k = 1:argIn.NumShuffles
                        auxX = auxResidsX;
                        for l = 1:size(auxX,3)
                            auxX(:,:,l) = auxX(:,randperm(T),l);
                        end
                        auxX = Squash( auxX )';
                        
                        [~, ~, r_s_j_k] ...
                            = CanonCorr( auxX, Yj );
                        r_s_j(k) = r_s_j_k(1);
                    end
                    
                case 3
                    for k = 1:argIn.NumShuffles
                        auxX = auxResidsX;
                        for l = 1:T
                            auxX(:,l,:) = auxX(:,l,randperm(numTrials));
                        end
                        auxX = Squash( auxX )';
                        
                        [~, ~, r_s_j_k] ...
                            = CanonCorr( auxX, Yj );
                        r_s_j(k) = r_s_j_k(1);
                    end
                    
                case 4
                    expCondId = unique(expCond);
                    numExpCond = numel(expCondId);
                    for k = 1:argIn.NumShuffles
                        auxX = auxResidsX;
                        for l = 1:T
                            for expCondIdx = 1:numExpCond
                                idxs = ...
                                    find(expCond == expCondId(expCondIdx));
                                randomize = randperm( numel(idxs) );
                                auxX( :,l,idxs ) ...
                                    = auxX( :,l,idxs(randomize) );
                            end
                        end
                        auxX = Squash( auxX )';
                        
                        [~, ~, r_s_j_k] ...
                            = CanonCorr( auxX, Yj );
                        r_s_j(k) = r_s_j_k(1);
                    end
                    
                case 5
                    expCondId = unique(expCond);
                    numExpCond = numel(expCondId);
                    for k = 1:argIn.NumShuffles
                        jSpikeRasters = zeros( size(auxSpikeRastersX) );
                        for expCondIdx = 1:numExpCond
                            idxs = find(expCond == expCondId(expCondIdx));
                            
                            jSpikeRasters( :,:,idxs ) = JitterSpikes...
                                ( auxSpikeRastersX(:,:,idxs), ...
                                argIn.JitterWindow );
                        end
                        [~, jResids] ...
                            = ComputePsths({jSpikeRasters}, expCond);
                        jResids{1} ...
                            = BinTime(jResids{1}, argIn.BinWidth);
                        auxX = Squash( jResids{1} )';
                        
                        [~, ~, r_s_j_k] ...
                            = CanonCorr( auxX, Yj );
                        r_s_j(k) = r_s_j_k(1);
                    end
                
            end
            r_s(j,:) = r_s_j;
        else
            
            auxIdxs = randperm( numUnits(pop1) ); %#ok<PFBNS>
            unitIdxsX = auxIdxs( 1:ceil(numUnits(pop1)/2) );
            unitIdxsY = auxIdxs( ceil(numUnits(pop1)/2) + 1:end );
            
            Xn = X(:,unitIdxsX);
            Yjn = Yj(:,unitIdxsY);
            
            r_j = zeros(argIn.NumDraws, numCanonPairs);
            for k = 1:argIn.NumDraws
                [~, ~, r_j_k] = CanonCorr( Xn, Yjn );
                r_j(k,:) = r_j_k(1:numCanonPairs);
            end
            r(j,:) = mean(r_j, 1);
            
        end
        
    end
    
    argOut.CorrMap                          ...
        (t2 - t1 + argIn.MaxDelay+1,i,:)    ...
        = permute(r, [1 3 2]);
    
    argOut.FrMap                            ...
        (t2 - t1 + argIn.MaxDelay+1,i,:)    ...
        = fr;
    
    argOut.LowRankWarning                   ...
        (t2 - t1 + argIn.MaxDelay+1,i,:)    ...
        = warnings;
    
    argOut.ShuffledCorrMap                  ...
        (t2 - t1 + argIn.MaxDelay+1,i,:)    ...
        = permute(r_s, [1 3 2]);
    
    iterTime(i) = toc;
    
    DeleteMsg(msg);
    msg = PrintProgress(i, numel(timePoints), iterTime);
    
end

fprintf('\n')

end



function [spikes, unitsRmvd] = RmvSilentUnits...
	(spikes, timePoints, binWidth, maxDelay, windowLength)

pop1 = 1;
pop2 = 2;

numPops = numel(spikes);

numUnits = zeros(numPops, 1);
for popIdx = 1:numPops
	numUnits(popIdx) = size(spikes{popIdx}, 1);
end

trialLength = size(spikes{pop1}, 2);

mapDim = trialLength - (windowLength-1);

unitsRmvd = cell(numPops, 1);
unitsRmvd{pop1} = zeros(numUnits(pop1), 1);
unitsRmvd{pop2} = zeros(numUnits(pop2), 1);
for t1 = timePoints
    
    X = Squash( BinTime( 						...
        spikes{pop1}(:,t1:t1+windowLength-1,:)	...
        , binWidth ) );
    
    unitsRmvd{pop1}(std(X, 0, 2) == 0) = 1;
    
    [~,I,~] = unique(X, 'rows');
    ixDupRows = setdiff(1:size(X,1), I);
    
    unitsRmvd{pop1}(ixDupRows) = 1;
    
    for t2 = max(t1 - maxDelay, 1):min(t1 + maxDelay, mapDim)
        
        Y = Squash( BinTime( 						...
            spikes{pop2}(:,t2:t2+windowLength-1,:)	...
            , binWidth ) );
        
        unitsRmvd{pop2}(std(Y, 0, 2) == 0) = 1;
        
        [~,I,~] = unique(Y, 'rows');
        ixDupRows = setdiff(1:size(Y,1), I);
        
        unitsRmvd{pop2}(ixDupRows) = 1;
        
    end
    
end

for popIdx = 1:numPops

	unitsRmvd{popIdx} = (unitsRmvd{popIdx} == 1);
	
	spikes{popIdx}(unitsRmvd{popIdx},:,:) = [];

end

end


