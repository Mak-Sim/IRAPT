function [ F0, time_marks ] = irapt(Sig, Fs, type)
%IRAPT -- implementation of instanteneous RAPT algorithm.
%  Sig -- input signal
%  Fs  -- sampling frequency
%  F0  -- f0 estimations
%  time_marks -- corresponding time marks (in seconds)
%  type -- 'irapt1' or 'irapt2'

if(Fs~=44100)
    Sig=resample(Sig,44100,Fs);
    Fs=44100;
end

[F0,~,Cfg] = irapt1(Sig');
if (strcmpi(type,'ipart2'))
    F0 = irapt2(Cfg,Sig',F0);   %IRAPT2
end

time_marks = (0:length(F0)-1)*Cfg.step_smp/Fs;

end

