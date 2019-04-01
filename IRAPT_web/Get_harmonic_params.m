function [ Amps, Frc, Phs, Zo ] = Get_harmonic_params(x,E_mx,M,S,Zi,ind_array)
%x - samples packet of length M (column vector)
%E_mx - polyphase matrix
%M - number of channels
%S - subsampling factor
%Zi - input filters states (column vector)
%Amps, Frc, Phs  - output params
%Zo - output filters states (column vector)

[fb_out1, Zi]=F0_Bank_get_samples(x(1:S-1),E_mx,M,S-1,Zi);
[fb_out2, Zi]=F0_Bank_get_samples(x(S),E_mx,M,1,Zi);

y1 = sum(fb_out1(ind_array),1).';
y2 = sum(fb_out2(ind_array),1).';  

Zo = Zi;

Amps=abs(y2)*2;
Phs=angle(y2);
Frc=diff(My_unwrap([angle(y1) Phs]),1,2);

end



