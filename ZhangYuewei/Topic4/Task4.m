clc;clear;close all;

%% Input 3 data files
dataVec1 = load("data1.txt","-ascii")';
dataVec2 = load("data2.txt","-ascii")';
dataVec3 = load("data3.txt","-ascii")';

%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = length(dataVec1);
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Supply PSD values
% This is the noise psd we will use.
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% 
% Signal parameter
a1=10;
a2=3;
a3=3;

%% Calculate GLRT value
llr1 = glrtqcsig(dataVec1,[a1,a2,a3],psdPosFreq,sampFreq,nSamples);
disp("GLRT value for data1 =" + llr1);
llr2 = glrtqcsig(dataVec2,[a1,a2,a3],psdPosFreq,sampFreq,nSamples);
disp("GLRT value for data1 =" + llr2);
llr3 = glrtqcsig(dataVec3,[a1,a2,a3],psdPosFreq,sampFreq,nSamples);
disp("GLRT value for data1 =" + llr3);

%% Estimate the significance of data1
M = 10000;
num1 = 0; 
num2 = 0;
num3 = 0;
for i=1:M
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    testllr = glrtqcsig(noiseVec,[a1,a2,a3],psdPosFreq,sampFreq,nSamples);
    if testllr >= llr1
        num1 = num1 + 1;
    end
    if testllr >= llr2
        num2 = num2 + 1;
    end
    if testllr >= llr3
        num3 = num3 + 1;
    end
end
alpha1 = num1 / M;
alpha2 = num2 / M;
alpha3 = num3 / M;
disp("Significance value for data1 = " + alpha1);
disp("Significance value for data2 = " + alpha2);
disp("Significance value for data3 = " + alpha3);