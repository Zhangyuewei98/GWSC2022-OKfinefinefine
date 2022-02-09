%% Band pass filter demo
addpath ../SIGNALS;
sampFreq = 1024;
nSamples = 2048;

timeVec = (0:(nSamples-1))/sampFreq;

%% MIX3 signal
% Signal parameters
A1 = 10;
A2 = 5;
A3 = 2.5;
f1 = 100;
f2 = 200;
f3 = 300;
phai1 = 0;
phai2 = pi/6;
phai3 = pi/4;
mix_snr = 10;
% Signal length
sigLen = (nSamples-1)/sampFreq;

% Generate signal
sigVec = crcbgenMIX3sig(timeVec,mix_snr,[A1,A2,A3],[f1,f2,f3],[phai1,phai2,phai3]);

%% Remove frequencies 
% Design low pass filter
filtOrdr = 30;
%Change the frequency band you need

%w_low =90;
%w_high = 110;
%w_low =190;
%w_high = 210;
w_low =290;
w_high = 310;
wn = [w_low/sampFreq w_high/sampFreq];
b = fir1(filtOrdr,wn);
% Apply filter
filtSig = fftfilt(b,sigVec);

%% Plots
figure;
hold on;
plot(timeVec,sigVec);
plot(timeVec,filtSig);