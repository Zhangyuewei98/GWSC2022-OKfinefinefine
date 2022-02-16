clc;clear;close all
%% Plot the sinusoidal signal
% Signal parameters
A = 10;
ta = 1;
f0 = 1;
f1 = 3;
phi0 = 0.5 * pi;
L = 2;

% Stop time of the whole observing time
stoptime = ta + L + 1.0;

% Instantaneous frequency after 1 sec is 
maxFreq = 2*pi*(f0+2*f1*(ta+L)-2*f1*ta);  % f(t) = 2*pi*(f0+2*f1*t-2*f1*ta). maxFreq = f(t=ta+L)
disp('maximum value of instantaneous frequency(Hz):')
disp(maxFreq/2/pi)

% Nyquist sampling frequency
NyquistFreq = 2 * maxFreq;
disp('Nyquist sampling frequency(Hz):')
disp(NyquistFreq/2/pi)

% Set sample frequency
samplFreq1= NyquistFreq;
samplFreq2 = 5*NyquistFreq;
samplFreq3 = 0.5*NyquistFreq;
% Calculate sanpling interval delta
samplIntrvl1 = 1/samplFreq1;
samplIntrvl2 = 1/samplFreq2;
samplIntrvl3 = 1/samplFreq3;

% Time samples
timeVec1 = 0:samplIntrvl1:stoptime;
timeVec2 = 0:samplIntrvl2:stoptime;
timeVec3 = 0:samplIntrvl3:stoptime;
% Number of samples
nSamples1 = length(timeVec1);
nSamples2 = length(timeVec2);
nSamples3 = length(timeVec3);

% Generate the signal
sigVec1 = crcbgenltcsig(timeVec1,A,[ta,f0,f1,phi0,L]);
sigVec2 = crcbgenltcsig(timeVec2,A,[ta,f0,f1,phi0,L]);
sigVec3 = crcbgenltcsig(timeVec3,A,[ta,f0,f1,phi0,L]);

%Plot the signal
p1 = figure;
plot(timeVec1,sigVec1,'Marker','.','MarkerSize',8);
title('Linear transient chirp signal(fs=NyquistFrequency)');
xlabel('Time');
ylabel('Amplitude');
print(p1,'-djpeg','Linear transient chirp signal(fs=NyquistFrequency).jpg');
% another way to save picture:saveas(gcf, 'Linear transient chirp signal(fs=NyquistFrequency).jpg') 
p2 = figure;
plot(timeVec2,sigVec2,'Marker','.','MarkerSize',8);
title('Linear transient chirp signal(fs=5*NyquistFrequency)');
xlabel('Time');
ylabel('Amplitude');
print(p2,'-djpeg','Linear transient chirp signal(fs=5NyquistFrequency).jpg');
p3 = figure;
plot(timeVec3,sigVec3,'Marker','.','MarkerSize',8);
title('Linear transient chirp signal(fs=0.5*NyquistFrequency)');
xlabel('Time');
ylabel('Amplitude');
print(p3,'-djpeg','Linear transient chirp signal(fs=0.5NyquistFrequency).jpg');

%Plot the periodogram
%--------------
%Length of data 
dataLen1 = timeVec1(end)-timeVec1(1)+samplIntrvl1;
dataLen2 = timeVec2(end)-timeVec2(1)+samplIntrvl2;
dataLen3 = timeVec3(end)-timeVec3(1)+samplIntrvl3;
%DFT sample corresponding to Nyquist frequency
kNyq1 = floor(nSamples1/2)+1;
kNyq2 = floor(nSamples2/2)+1;
kNyq3 = floor(nSamples3/2)+1;
% Positive Fourier frequencies
posFreq1 = (0:(kNyq1-1))*(1/dataLen1);
posFreq2 = (0:(kNyq2-1))*(1/dataLen2);
posFreq3 = (0:(kNyq3-1))*(1/dataLen3);
% FFT of signal
fftSig1 = fft(sigVec1);
fftSig2 = fft(sigVec2);
fftSig3 = fft(sigVec3);
% Discard negative frequencies
fftSig1 = fftSig1(1:kNyq1);
fftSig2 = fftSig2(1:kNyq2);
fftSig3 = fftSig3(1:kNyq3);

%Plot periodogram
p4 = figure;
plot(posFreq1,abs(fftSig1));
xlim([0,30]);
title('Periodogram(fs=NyquistFrequency)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
print(p4,'-djpeg','Periodogram(fs=NyquistFrequency).jpg');
p5 = figure;
plot(posFreq2,abs(fftSig2));
xlim([0,30]);
title('Periodogram(fs=5*NyquistFrequency)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
print(p5,'-djpeg','Periodogram(fs=5NyquistFrequency).jpg');
p6 = figure;
plot(posFreq3,abs(fftSig3));
xlim([0,30]);
title('Periodogram(fs=0.5*NyquistFrequency)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
print(p6,'-djpeg','Periodogram(fs=0.5NyquistFrequency).jpg');

%Plot a spectrogram
%----------------
winLen = 0.2;%sec
ovrlp = 0.1;%sec
%1.Convert to integer number of samples(fs=NyquistFrequency) 
winLenSmpls1 = floor(winLen*samplFreq1);
ovrlpSmpls1 = floor(ovrlp*samplFreq1);
[S1,F1,T1]=spectrogram(sigVec1,winLenSmpls1,ovrlpSmpls1,[],samplFreq1);
p7 = figure;
imagesc(T1,F1,abs(S1)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Spectrogram(fs=NyquistFrequency)');
print(p7,'-djpeg','Spectrogram(fs=NyquistFrequency).jpg');
%2.Convert to integer number of samples(fs=5*NyquistFrequency) 
winLenSmpls2 = floor(winLen*samplFreq2);
ovrlpSmpls2 = floor(ovrlp*samplFreq2);
[S2,F2,T2]=spectrogram(sigVec2,winLenSmpls2,ovrlpSmpls2,[],samplFreq2);
p8 = figure;
imagesc(T2,F2,abs(S2)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Spectrogram(fs=5*NyquistFrequency)');
print(p8,'-djpeg','Spectrogram(fs=5NyquistFrequency).jpg');
%3.Convert to integer number of samples(fs=0.5*NyquistFrequency) 
winLenSmpls3 = floor(winLen*samplFreq3);
ovrlpSmpls3 = floor(ovrlp*samplFreq3);
[S3,F3,T3]=spectrogram(sigVec3,winLenSmpls3,ovrlpSmpls3,[],samplFreq3);
p9 = figure;
imagesc(T3,F3,abs(S3)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Spectrogram(fs=0.5*NyquistFrequency)');
print(p9,'-djpeg','Spectrogram(fs=0.5NyquistFrequency).jpg');