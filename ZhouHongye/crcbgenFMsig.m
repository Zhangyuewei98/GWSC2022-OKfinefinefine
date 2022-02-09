function sigVec = crcbgenFMsig(dataX,snr,FMCoefs)
% Generate a  Frequency modulated (FM) sinusoid signal
% S = CRCBGENFMSIG(X,SNR,C)
% Generates a Frequency modulated (FM) sinusoid signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% four coefficients [A,b,f0,f1] that parametrize the signal:


%Zhou Hongye, Feb. 2022
sigVec = FMCoefs(1)*sin(2*pi*FMCoefs(3)*dataX+FMCoefs(2)*cos(2*pi*FMCoefs(4)*dataX));
sigVec = snr*sigVec/norm(sigVec);