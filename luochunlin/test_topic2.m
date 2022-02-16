%% lab topic2
% signal parameters
A = 25;
B = 10;
f0 = 1e-5;
phi0 = pi/3;
% detector parameter
theta = [1/6,1/5,1/4]*pi;
phi = [1/3,1/2,2/3]*pi;
psi = theta;

%f_max at t=1
f_max = f0;
f_sampl = 5 * f0;
%sampling dt
dt = 1/f_sampl;

%time vector
t = 0:dt:365*24*3600;
%number of samples
n_sampl = length(t);

%generate signal
hp = gen_sin_sig(t,f0,A,0);
hc = gen_sin_sig(t,f0,B,phi0);

figure(1);
plot(t,hp,t,hc);
legend('h+','hx')
title('hplus and hcross of the origin signal')
%fp,fc
[fp,fc]=det_frame_rot_fpfc(theta,phi,psi);
%strain
for i = 1:length(fp)
    s(i,:) = hp * fp(i) + hc * fc(i);
end

%plot
figure(2);
plot(t,s(1,:),t,s(2,:),t,s(3,:))
legend('s1','s2','s3')
title('land for strain1/strain2/strain3')

%% LISA response
%LISA fp fc
[Fp,Fc] = LISA_FpFc(theta,phi,psi,t);
strain1 = Fp(1,:,1).*hp+Fc(1,:,1).*hc;
strain2 = Fp(1,:,2).*hp+Fc(1,:,2).*hc;
strain3 = Fp(3,:,1).*hp+Fc(3,:,1).*hc;      %different location compare to strain1

figure(3)
plot(t,strain1)
legend('tdi-s1')
title('TDI strain1')

n_sampl = length(t);
kNyq = floor(n_sampl/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/t(end));
fft_hp = fft(hp);
fft_hc = fft(hc);
fft_strain1 = fft(strain1);
fft_strain2 = fft(strain2);
fft_strain3 = fft(strain3);
fft_hp = fft_hp(1:kNyq);
fft_hc = fft_hc(1:kNyq);
fft_strain1 = fft_strain1(1:kNyq);
fft_strain2 = fft_strain2(1:kNyq);
fft_strain3 = fft_strain3(1:kNyq);
%plot
figure(4);
plot(posFreq,abs(fft_hp),posFreq,abs(fft_strain1),posFreq,abs(fft_strain3))
legend('ffth+','fft-s1','fft-s3')
title('FFT of the signal')

