function sig = gen_sin_sig(t,f0,A,phi0)

sig = A * sin(2*pi * f0 * t + phi0);
end