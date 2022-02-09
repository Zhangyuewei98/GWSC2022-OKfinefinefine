function sigVec = crcbgenAMsig(dataX,snr,AMCoefs)
% Generate a  Frequency modulated (FM) sinusoid signal
% S = CRCBGENAMSIG(X,SNR,C)
% Generates a Frequency modulated (FM) sinusoid signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% four coefficients [A,phai0,f0,f1] that parametrize the signal:


%Zhou Hongye, Feb. 2022
sigVec = FMCoefs(1)*cos(2*pi*AMCoefs(4)*dataX)*sin(AMCoefs(4)*dataX+AMCoefs(2));
sigVec = snr*sigVec/norm(sigVec);