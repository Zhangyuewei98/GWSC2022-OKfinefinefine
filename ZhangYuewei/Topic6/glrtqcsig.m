function llr = glrtqcsig(dataVec,qcCoefs, psdPosFreq, sampFreq, nSamples)
% Calculate the GLRT for a quadratic chirp signal with unknown amplitude
% LLR = GLRTQCSIG(X,C,F,S,N)
% Calculate the GLRT for a data vector X with a quadratic chirp signal 
% whose parameter vector is C = [a1,a2,a3]. The frequency vector is F, the
% sampling rate is S and the number of sampling points is N 

% Yuewei Zhang, January 2022

timeVec = (0:(nSamples-1))/sampFreq;
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
% disp("GLRT value = " + llr);