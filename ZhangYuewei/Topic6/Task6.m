clc;clear;close all;

%% Load noise and analysis data
noise = load("TrainingDataTF.mat");
noiseVec = noise.trainData;
data = load("analysisDataTF.mat");
dataVec = data.dataVec;

%% Parameters for data realization
% Number of samples and sampling frequency.
nRuns = 8;
nSamples = length(dataVec);
sampFreq = noise.sampFreq;
timeVec = (0:(nSamples-1))/sampFreq;

%% Estimate power spectral density
[psdPosFreq, posFreq] = pwelch(noiseVec,nSamples,[],[],sampFreq);
psdPosFreq = psdPosFreq';
posFreq = posFreq';
figure;
plot(posFreq,psdPosFreq);
axis([0,posFreq(end),0,max(psdPosFreq)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');

%% Set rmin and rmax
rmin = [40, 1, 1];
rmax = [100, 50, 15];

%% Calculate data vector
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
llrmax = outResult.bestFitness * (-1);

%% Detection(Calculate significance)
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
num = 0;
sigVec = crcbgenqcsig(timeVec,1,[a1PSO,a2PSO,a3PSO]);
[templateVec,~] = normsig4psd(sigVec,sampFreq,psdPosFreq,1);
llrH0_Vec = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    llrH0 = innerprodpsd(noiseVec,templateVec,sampFreq,psdPosFreq);
    llrH0_Vec(lp) = llrH0;
    if llrH0 >= llrmax
        num = num + 1;
    end
end
nH1Data = 1000;
llrH1_Vec = zeros(1,nH1Data);
for lp = 1:nH1Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    llrH1_Vec(lp) = innerprodpsd(dataVec,templateVec,sampFreq,psdPosFreq);
end

%% Calculate significance
significance = num / nH0Data;
if significance < 0.05
    disp("There is a signal in the analysis data.");
end
disp("significance = " + significance);

%% Calculate SNR
estSNR = (mean(llrH1_Vec)-mean(llrH0_Vec))/std(llrH0_Vec);
ActualParams = load("keyFileTF.mat");
snr = ActualParams.snr;
disp("Targeted SNR = " + snr);
disp("Estimated SNR =" + estSNR);

%% Load Actual parameters [a1,a2,a3]
a1 = ActualParams.a1;
a2 = ActualParams.a2;
a3 = ActualParams.a3;
disp("Targeted parameter:")
disp("a1 = " + a1);
disp("a2 = " + a2);
disp("a3 = " + a3);
disp("PSO estimated parameter:")
disp("a1 = " + a1PSO);
disp("a2 = " + a2PSO);
disp("a3 = " + a3PSO);