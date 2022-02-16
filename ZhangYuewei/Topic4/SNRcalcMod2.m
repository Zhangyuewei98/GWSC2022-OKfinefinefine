clc;clear;close all;
%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

%%
% Path to folder containing signal and noise generation codes
% addpath ../SIGNALS
% addpath ../NOISE

%%
% This is the target SNR for the LR
snr = 10;

%%
% Data generation parameters
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;


%%
% Generate the signal that is to be normalized
ta = 1;
f0 = 1;
f1 = 3;
phi0 = 0.5 * pi;
L = 2;
% Amplitude value does not matter as it will be changed in the normalization
A = 10; 
sigVec = crcbgenltcsig(timeVec,1,[ta,f0,f1,phi0,L]);

%%
% We will use the noise PSD of initial LIGO design sensitivity PSD.
% noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
LIGOPSD = load('iLIGOSensitivity.txt','-ascii');
% Modify the LIGO PSD
LIGOPSD_Mod = zeros(size(LIGOPSD));
LIGOPSD_Mod(:,1) = LIGOPSD(:,1);
for i=1:length(LIGOPSD(:,1))
    if (LIGOPSD(i,1)<=50) && (LIGOPSD(i+1,1)>50)
        LIGOPSD_Mod(1:i,2) = LIGOPSD(i+1,2)*ones(i,1);
    elseif (LIGOPSD(i,1)<700) && (LIGOPSD(i+1,1)>=700)
        LIGOPSD_Mod(i+1:end,2) = LIGOPSD(i,2)*ones(length(LIGOPSD(:,1))-i,1);
    else
        LIGOPSD_Mod(i,2) = LIGOPSD(i,2);
    end
end

%%
% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies. 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
% psdPosFreq = noisePSD(posFreq);
% figure;
% plot(posFreq,psdPosFreq);
% axis([0,posFreq(end),0,max(psdPosFreq)]);
% xlabel('Frequency (Hz)');
% ylabel('PSD ((data unit)^2/Hz)');
% title('Noise Power Spectral Density')

%% Interpolate LIGOPSD according to posFreq
LIGOpsdPosFreq = interp1(LIGOPSD_Mod(:,1),LIGOPSD_Mod(:,2),posFreq,'linear','extrap');
figure;
plot(posFreq,LIGOpsdPosFreq);
axis([0,posFreq(end),0,max(LIGOpsdPosFreq)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');
title('LIGO Noise Power Spectral Density')

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,LIGOpsdPosFreq);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),LIGOpsdPosFreq(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,LIGOpsdPosFreq);
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),LIGOpsdPosFreq(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,LIGOpsdPosFreq);
end
%%
% Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);
disp("Targeted SNR = " + snr);
disp("Estimated SNR =" + estSNR);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
legend('H_0','H_1');
title(['Estimated SNR = ',num2str(estSNR)]);

%%
% A signal realization
figure;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Noise');
title('Signal Realization')
%%
% A noise realization
figure;
plot(timeVec,noiseVec);
xlabel('Time (sec)');
ylabel('Noise');
title('Noise Realization')
%%
% A data realization
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Data');
title('Data Realization')
%%
% periodgram of the signal
figure;
fftSig = fft(sigVec);
fftSig = fftSig(1:kNyq);
plot(posFreq,abs(fftSig));
xlim([0,30]);
title('Signal Periodogram');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
%%
% periodgram of the data
figure;
fftData = fft(dataVec);
fftData = fftData(1:kNyq);
plot(posFreq,abs(fftData));
xlim([0,30]);
title('Data Periodogram');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
%%
% spectrogram of the signal
figure;
winLen = 0.2;%sec
ovrlp = 0.1;%sec
%1.Convert to integer number of samples(fs=NyquistFrequency) 
winLenSmpls = floor(winLen*sampFreq);
ovrlpSmpls = floor(ovrlp*sampFreq);
[S,F,T]=spectrogram(sigVec,winLenSmpls,ovrlpSmpls,[],sampFreq);
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Signal Spectrogram');
%%
% spectrogram of the data
figure;
winLen = 0.2;%sec
ovrlp = 0.1;%sec
%1.Convert to integer number of samples(fs=NyquistFrequency) 
winLenSmpls = floor(winLen*sampFreq);
ovrlpSmpls = floor(ovrlp*sampFreq);
[S,F,T]=spectrogram(dataVec,winLenSmpls,ovrlpSmpls,[],sampFreq);
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Data Spectrogram');

