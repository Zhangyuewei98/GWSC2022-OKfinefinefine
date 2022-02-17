function llr = glrtqcsig_Mod(qcCoefs, params)
% function llr = glrtqcsig_Mod(dataVec, ,qcCoefs, psdPosFreq, sampFreq, nSamples)
% Calculate the GLRT for a quadratic chirp signal with unknown amplitude
% LLR = GLRTQCSIG_MOD(X, P)
% Calculate the GLRT for a quadratic chirp signal data whose parameter 
% vector [a1, a2, a3] = X with a parameter vector is P. 
% The fields of P are:
% 'dataY': The data vector (a uniformly sampled time series).
% 'dataX': The time stamps of the data samples.
% 'dataXSq': dataX.^2
% 'dataXCb': dataX.^3
% 'sampFreq': sampling frequency
% 'psdPosFreq': power spectral density
% 'rmin', 'rmax': The minimum and maximum values of the three parameters
%                 a1, a2, a3 in the candidate signal:
%                 a1*dataX+a2*dataXSq+a3*dataXCb
% rmin and rmax are used
%to convert X(i,j) internally before computing the fitness:
%X(:,j) -> X(:,j)*(rmax(j)-rmin(j))+rmin(j). 

% Yuewei Zhang, January 2022

nSamples = length(params.dataX);
sampFreq = params.sampFreq;
psdPosFreq = params.psdPosFreq;
timeVec = (0:(nSamples-1))/sampFreq;

dataVec = params.dataY;

% the value used for 'A' does not matter because we are going to normalize 
% the signal anyway.
A = 1; 
sigVec = crcbgenqcsig(timeVec,A,qcCoefs);
%We do not need the normalization factor, just the  template vector
[templateVec,~] = normsig4psd(sigVec,sampFreq,psdPosFreq,1);
% Calculate inner product of data with template
llr = innerprodpsd(dataVec,templateVec,sampFreq,psdPosFreq);
%GLRT is its square
llr = llr^2;
llr = -1 * llr;    % transform maximizing llr to minimizing -llr 
% disp("GLRT value = " + llr);