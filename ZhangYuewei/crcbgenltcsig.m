function sigVec = crcbgenltcsig(dataX,snr,sinCoefs)
% Generate a linear transient chirp signal
% S = CRCBGENSINSIG(X,SNR,C)
% Generates a linear transient chirp signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% five coefficients [ta,f0,f1,phi0,L] that parametrize the phase of the signal:
% 2*pi*(f0*(t-ta)+f1*(t-ta)^2)+phi0 between the time interval:[ta,ta+L], while
% the signal value outside this time interval is 0.

%Yuewei Zhang, February 2022

phaseVec = 2*pi*(sinCoefs(2)*(dataX-sinCoefs(1))+sinCoefs(3)*(dataX-sinCoefs(1)).^2)+sinCoefs(4);
sigVec = sin(phaseVec);
sigVec = snr*sigVec/norm(sigVec);
for i =1:length(dataX)
    if dataX(i) < sinCoefs(1) || dataX(i) > (sinCoefs(1)+sinCoefs(5))
        sigVec(i) = 0;
    end
end


