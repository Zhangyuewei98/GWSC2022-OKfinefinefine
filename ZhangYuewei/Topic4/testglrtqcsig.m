clc;clear;close all;

%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Supply PSD values
% This is the noise psd we will use.
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% Generate  data realization
% Noise + SNR=10 signal. 
a1=9.5;
a2=2.8;
a3=3.2;
A=10;
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
sig4data = crcbgenqcsig(timeVec,A,[a1,a2,a3]);
% Signal normalized to SNR=10
[sig4data,~]=normsig4psd(sig4data,sampFreq,psdPosFreq,10);
dataVec = noiseVec+sig4data;

figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sig4data);
xlabel('Time (sec)')
ylabel('Data');
title('Data realization for calculation of LR');

figure;
kNyq = floor(nSamples/2)+1;
dataLen = nSamples/sampFreq;
posFreq = (0:(kNyq-1))*(1/dataLen);
datFFT = abs(fft(dataVec));
sigFFT = abs(fft(sig4data));
plot(posFreq,datFFT(1:kNyq));
hold on;
plot(posFreq,sigFFT(1:kNyq));
xlabel('Frequency (Hz)');
ylabel('Periodogram');

figure;
[S,F,T] = spectrogram(dataVec,64,60,[],sampFreq);
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)')
ylabel('Frequency (Hz)');

%% Compute GLRT
llr = glrtqcsig(dataVec,[a1,a2,a3],psdPosFreq,sampFreq,nSamples); 
% set llr value as observed value
num = 0;
N = 500;
for i=1:N
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    testllr = glrtqcsig(noiseVec,[a1,a2,a3],psdPosFreq,sampFreq,nSamples);
    if testllr >= llr
        num = num + 1;
    end
end
alpha = num / N;
disp("alpha = " + alpha);