function argOut = CovStabilityAcrossTimeAnalysis...
    (spikes, expCond, argIn)
% 
% argOut = CovStabilityAcrossTimeAnalysis(spikes, expCond, argIn) computes 
% a interaction structure similarity.
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
% expCond - a vector containing the trial label for each trial (N x 1)
% 
% argIn   - structure containing additional input arguments.
% 	argIn =
%       TimePeriods:    a matrix indicating the time intervals to be
%       considered for each trial period. Each row corresponds to a period
%       with the first column indicating start time and the seconds column
%       indicating end time. Units are indexes into the second dimension of
%       spikes{1} and spikes {2}
% 
% OUTPUTS:
% 
% argOut - structure containing outputs
% 	argOut =
%       CvR: a 3D array containing the captured correlations. The first 
%       dimension corresponds to the time period the CCA model was fit to. 
%       The second dimensions corresponds to the time period the CCA model 
%       was tested on. The third dimension indicates the dimensionality of 
%       the CCA model.
% 
%       cvLogLike: a matrix containing the cross-validated log-likelihood 
%       of pCCA for each time period. Rows correspond to the different time
%       periods, columns to the different model dimensionalities.

pop1 = 1;
pop2 = 2;

if isfield(argIn, 'FixedDims')
	fixedDims = argIn.FixedDims;
else
	fixedDims = [];
end

spikes{pop1} = RmvSilentUnits(spikes{pop1}, argIn.TimePeriods);
spikes{pop2} = RmvSilentUnits(spikes{pop2}, argIn.TimePeriods);

numTimePeriods = size(argIn.TimePeriods, 1);
lenTimePeriods = range(argIn.TimePeriods(1,:));

numTrials = size(spikes{pop1}, 3);

p = size(spikes{pop1}, 1);
K = size(spikes{pop2}, 1);

[~, resids] = ComputePsths(spikes, expCond);

firingRates{pop1} = mean(spikes{pop1}, 3);
firingRates{pop2} = mean(spikes{pop2}, 3);

rank = 1:3;%min( min(p, K) - 1, C_MAX_RANK_BETA );

X = zeros(numTrials*lenTimePeriods, p, numTimePeriods);
for timeIdx = 1:numTimePeriods
    
    X(:,:,timeIdx) = Squash( resids{pop1}(:,...
        argIn.TimePeriods(timeIdx,1) + 1:argIn.TimePeriods(timeIdx,2),...
        :) )'; %#ok<*PFBNS>
    
end

% ===================================================================
%   Regress pop2 residuals on pop1 residuals
% ===================================================================

Y = zeros(numTrials*lenTimePeriods, K, numTimePeriods);
for timeIdx = 1:numTimePeriods
    
    Y(:,:,timeIdx) = Squash( resids{pop2}(:,...
        argIn.TimePeriods(timeIdx,1) + 1:argIn.TimePeriods(timeIdx,2),...
        :) )';
    
end

[cvR, cvLogLike, ccaVarX, ccaVarY, estParams] = CanonCorrStabilityAcrossTime...
	(X, Y, fixedDims, rank);

argOut = SaveArgOut(cvR, cvLogLike, ccaVarX, ccaVarY, estParams, rank, firingRates);

end



function spikes = RmvSilentUnits(spikes, timePeriods)

minThreshold = .001*range(timePeriods(1,:))*size(spikes, 3);

silentUnits = zeros(size(spikes, 1), 1);

numTimePeriods = size(timePeriods, 1);
for timeIdx = 1:numTimePeriods
    numSpikes = sum( Squash( spikes(:,...
        timePeriods(timeIdx,1) + 1:timePeriods(timeIdx,2),...
        :) ), 2 );
    
    silentUnits( numSpikes < minThreshold ) = 1;
end

silentUnits = (silentUnits == 1);

spikes(silentUnits,:,:) = [];

end



function [cvR, cvLogLike, ccaVarX, ccaVarY, estParams] = CanonCorrStabilityAcrossTime...
	(X, Y, fixedDims, rank)

C_NUM_FOLDS = 10; % Num of cross-validation folds

N = size(X, 1);
cvp = cvpartition(N, 'kFold', C_NUM_FOLDS);

numCases = size(X, 3);

numelRank = numel(rank);

cvR = zeros(numCases, numCases, numelRank, C_NUM_FOLDS);
cvLogLike = zeros(numCases, numelRank, C_NUM_FOLDS);
ccaVarX = zeros(numCases, min( size(X, 2), size(Y, 2) ), C_NUM_FOLDS);
ccaVarY = zeros(numCases, min( size(X, 2), size(Y, 2) ), C_NUM_FOLDS);
for fold = 1:C_NUM_FOLDS
    
    for i = 1:numCases
        
        Xtrain_i = X( cvp.training(fold), :, i );
        Ytrain_i = Y( cvp.training(fold), :, i );
        
        [Ai, Bi, r] = canoncorr(Xtrain_i,Ytrain_i);
        
        ccaVarX(i,:,fold) = VarExp(Xtrain_i, Ai);
        ccaVarY(i,:,fold) = VarExp(Ytrain_i, Bi);
        
        estParams_i = ComputeProbCanonCorrParams...
			(Xtrain_i, Ytrain_i, Ai, Bi, r, rank);
        
        for j = 1:numCases
			
            Xtrain_j = X( cvp.training(fold), :, j );
            Ytrain_j = Y( cvp.training(fold), :, j );
            
            [A, B, r] = canoncorr(Xtrain_j,Ytrain_j);
            
            estParams_j = ComputeProbCanonCorrParams...
                (Xtrain_j, Ytrain_j, A, B, r, rank);
			
            Xtest = X( cvp.test(fold), :, j );
            Ytest = Y( cvp.test(fold), :, j );
            
            cvR(i,j,:,fold) = ComputeSubspaceCorr(Xtest, Ytest, ...
				estParams_i, estParams_j, ...
				fixedDims, rank);
            
            if i == j
                cvLogLike(i,:,fold) ...
                    = ProbCanonCorrLogLike(Xtest, Ytest, estParams_i);
                
                ccaVarX(i,:,fold) = VarExp(Xtest, Ai);
                ccaVarY(i,:,fold) = VarExp(Ytest, Bi);
            end
            
        end
        
    end
    
end

cvR = mean(cvR, 4);
cvLogLike = mean(cvLogLike, 3);

estParams = cell(numCases, 1);
for i = 1:numCases
    
    [A, B, r] = canoncorr(X(:,:,i),Y(:,:,i));

    estParams{i} = ComputeProbCanonCorrParams...
        (X(:,:,i), Y(:,:,i), A, B, r, rank);
end

end



function subR = ComputeSubspaceCorr(Xtest, Ytest, ...
	estParams_i, estParams_j, fixedDims, rank)

SX = cov(Xtest);
SY = cov(Ytest);

subR = zeros(1, numel(rank));
for rankIdx = 1:numel(rank)
	
	if rank(rankIdx) == 0
		r = 0;
	else
		
		if isempty(fixedDims)
			W1 = estParams_i(rankIdx).W1;
			W2 = estParams_i(rankIdx).W2;
		elseif fixedDims == 1
			W1 = estParams_i(rankIdx).W1;
			W2 = estParams_j(rankIdx).W2;
		elseif fixedDims == 2
			W1 = estParams_j(rankIdx).W1;
			W2 = estParams_i(rankIdx).W2;
		end
		
		A = SX\W1;
		B = SY\W2;
		
% 		[V1, ~, ~] = svd(W1, 0);
% 		[V2, ~, ~] = svd(W2, 0);
% 		
% 		A = SX\V1;
% 		B = SY\V2;
		
		[~, ~, r] ...
			= canoncorr( Xtest*A, Ytest*B );
	end
	
	subR(rankIdx) = sum(r);
	
end

end


function var = VarExp(X, A)
numDims = size(A, 2);

V = pca(X);

[Q, ~] = qr(A);

var = zeros(1, numDims);
for dim = 1:numDims
    var(dim) = trace( cov(X*Q(:,1:dim)) )/trace( cov(X*V(:,1:dim)) );
end

end



function argOut = SaveArgOut(cvR, cvLogLike, ccaVarX, ccaVarY, estParams, rank, firingRates)

argOut.CvR = cvR;

argOut.CvLogLike = cvLogLike;

argOut.CCAVarX = ccaVarX;
argOut.CCAVarY = ccaVarY;

argOut.EstParams = estParams;

argOut.RankBeta = rank;

argOut.FiringRates = firingRates;

end % of SaveArgOut


