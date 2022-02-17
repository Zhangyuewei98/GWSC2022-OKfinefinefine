% %Demo for colored Gaussian noise generation
% %Sampling frequency for noise realization
% sampFreq = 1024; %Hz
% %Number of samples to generate
% nSamples = 16384;
% %Time samples
% timeVec = (0:(nSamples-1))/sampFreq;
% 
% %Target PSD given by the inline function handle
% targetPSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000;
% 
% %Plot PSD
% freqVec = 0:0.1:512;
% psdVec = targetPSD(freqVec);
% %plot(freqVec,psdVec);
% 
% %
% % Design FIR filter with T(f)= square root of target PSD
%  sqrtPSD = sqrt(psdVec);
%  fltrOrdr = 500;
%  b = fir2(fltrOrdr,freqVec/(sampFreq/2),sqrt(psdVec));
%  
% % %%
% % % Generate a WGN realization and pass it through the designed filter
% % % (Comment out the line below if new realizations of WGN are needed in each run of this script)
% % rng('default'); %默认种子
% % in_Noise = randn(1,nSamples);
% % out_Noise = fftfilt(b,inNoise);
% % figure;
% % plot(timeVec,out_Noise);
% % xlabel('time(t)');
% % ylabel('WGN');
% % Generate noise realization
% outNoise = statgaussnoisegen(nSamples,[freqVec(:),psdVec(:)],fltrOrdr,sampFreq);
% %
% % Estimate the PSD
%  figure;
%  pwelch(outNoise, 512,[],[],sampFreq);
%  hold on;
%  pwelch(inNoise, 512, [], [], sampFreq);
% %Pwelch plots in dB (= 10*log10(x)); plot on a linear scale
% [pxx,f]=pwelch(outNoise, 256,[],[],sampFreq);
% figure;
% plot(f,pxx);
% xlabel('Frequency (Hz)');
% ylabel('PSD');
% % Plot the colored noise realization
% figure;
% plot(timeVec,outNoise);
% xlabel('time (s)');
% ylabel('outNoise Magnitude');
% 
y = load('iLIGOSensitivity.txt','-ascii');

sampFreq = 20000; %Hz
nSamples = 16384;
%targetPSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000;
freqVec = y(:,1);
psdVec = y(:,2);
freqVec = [0,freqVec',sampFreq/2]';
psdVec = [0,psdVec',psdVec(97)*2]';
for i = 1:99
    if freqVec(i) < 50
        psdVec(i) = psdVec(42);
    end
    if freqVec(i) > 700
        psdVec(i) = psdVec(70);
    end
end
%psdVec(99) = sampFreq/2;
%psdVec = abs(log10(psdVec));
f_matrix = [freqVec,psdVec];

timeVec = (0:(nSamples-1))/sampFreq;
filt_ord = 500;
outNoise = colorGN(nSamples,f_matrix,filt_ord,sampFreq);
% plot(timeVec,outNoise);
% xlabel('time (s)');
% ylabel('outNoise Magnitude');

% Estimate the PSD
 figure;
 pwelch(outNoise, 2000,[],[],sampFreq);
 hold on;
 pwelch(inNoise, 2000, [], [], sampFreq);
%Pwelch plots in dB (= 10*log10(x)); plot on a linear scale
[pxx,f]=pwelch(outNoise, 2000,[],[],sampFreq);
figure;
plot(f,pxx);
xlabel('Frequency (Hz)');
ylabel('PSD');

function sigVec = colorGN(smp_num,f_matrix,filt_ord,fs)
%Sampling frequency for noise realization
sampFreq = fs; %Hz
%Number of samples to generate
nSamples = smp_num;
freqVec = f_matrix(:,1)';
psdVec = f_matrix(:,2)'.*f_matrix(:,2)';
 
% Generate noise realization
 sigVec = statgaussnoisegen(nSamples,[freqVec(:),psdVec(:)],filt_ord,sampFreq);
 end