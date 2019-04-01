function [y,Zo]=F0_Bank_get_samples(x,E_mx,M,S,Zi)
%x - samples packet of length M in the source order (column vector)
%E_mx - polyphase matrix
%M - number of channels
%S - subsampling factor
%Zi - input filters states (column vector)
%y  - subband signals (1 sample per channel)
%Zo - output filters states (column vector)


%Making delay chain matrix
Zi = [M*x(S:-1:1); Zi(1:end-S)];

%Polyphase filtration
D = reshape(Zi,M,length(Zi)/M);
P = sum(D.*E_mx,2);

%DFT-modulation
y = ifft(P);
% y = (1/M)*fft(P);

%Updating output memory 
Zo=Zi;