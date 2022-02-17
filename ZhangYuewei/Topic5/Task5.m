clc;clear;close all;

%% Parameters for data realization
% Number of samples and sampling frequency.
nRuns = 8;
snr = 10;
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%%
% Generate the signal that is to be normalized
a1=3;
a2=3;
a3=8;
% Amplitude value does not matter as it will be changed in the normalization
A = 10; 
sigVec = crcbgenqcsig(timeVec,1,[a1,a2,a3]);

%% Set rmin and rmax
rmin = [1, 2, 5];
rmax = [5, 6, 10];

%%
% We will use the noise PSD used in colGaussNoiseDemo.m but add a constant
% to remove the parts that are zero. (Exercise: Prove that if the noise PSD
% is zero at some frequencies but the signal added to the noise is not,
% then one can create a detection statistic with infinite SNR.)
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

%%
% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies. 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
figure;
plot(posFreq,psdPosFreq);
axis([0,posFreq(end),0,max(psdPosFreq)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,psdPosFreq);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Calculate noise vector
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);

%% Calculate data vector
dataVec = noiseVec + sigVec;
dataVecSq = dataVec.^2;
dataVecCb = dataVec.^3;


inParams = struct('dataY', dataVec, ...
                    'dataX', timeVec, ...
                    'dataXSq', dataVecSq, ...
                    'dataXCb', dataVecCb, ...
                    'sampFreq', sampFreq, ...
                    'psdPosFreq', psdPosFreq, ...
                    'rmin',rmin,...
                     'rmax',rmax ...
                  );


outResult = crcbqcpso_Mod(inParams, [], nRuns);
a1PSO = outResult.bestQcCoefs(1);
a2PSO = outResult.bestQcCoefs(2);
a3PSO = outResult.bestQcCoefs(3);
disp("Actual parameter:")
disp("a1 = " + a1);
disp("a2 = " + a2);
disp("a3 = " + a3);
disp("PSO parameter:")
disp("a1 = " + a1PSO);
disp("a2 = " + a2PSO);
disp("a3 = " + a3PSO);