function [ F0, Voc, time_marks ] = irapt(sig, fs, est_type, sig_type)
%IRAPT -- implementation of instanteneous RAPT algorithm.
%  sig -- input signal
%  fs  -- sampling frequency
%  time_marks -- corresponding time marks (in seconds)
%  est_type --  estimator type ('irapt1' or 'irapt2')
%  sig_type -- type of input signal ('speech' or 'sustain phonation')
%  F0  -- f0 estimations

if(fs~=44100)
    sig=resample(sig,44100,fs);
    fs=44100;
end

% Stage 1
if (strcmpi(sig_type,'speech'))
    [F0,Voc,Cfg] = irapt1(sig');
else
    [F0,Voc,Cfg] = irapt_sus(sig');
end

% Stage 2
if (strcmpi(est_type,'irapt2'))
    F0 = irapt2(Cfg,sig',F0);   %IRAPT2
end

time_marks = (0:length(F0)-1)*Cfg.step_smp/fs;

end

