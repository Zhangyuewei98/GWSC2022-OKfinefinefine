function sigVec = crcbgenMIX3sig(dataX,snr,A_list,f_list,phai_list)
% Generate a  Frequency modulated (FM) sinusoid signal
% S = CRCBGENFMSIG(X,SNR,C)
% Generates a Frequency modulated (FM) sinusoid signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% four coefficients [A,b,f0,f1] that parametrize the signal:


%Zhou Hongye, Feb. 2022
sig1 = A_list(1)*sin(2*pi*f_list(1)*dataX+phai_list(1));
sig2 = A_list(2)*sin(2*pi*f_list(2)*dataX+phai_list(2));
sig3 = A_list(3)*sin(2*pi*f_list(3)*dataX+phai_list(3));
sigVec = sig1+sig2+sig3;
sigVec = snr*sigVec/norm(sigVec);