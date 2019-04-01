function [y, Zo] = dft_fb(x,E_mx,M,oflag,Zi)
%DFT_FB Implements filtering by dft-modulated filter bank
%   x - input signal
%   E_mx - polyphase matrix
%   M - number of channels
%   S - subsampling factor
%   oflag - flag: 0 - no output, 1 - calculate output;
%   Zi - input filters states (column vector)
%   y  - subband signals (1 sample per channel)
%   Zo - output filters states (column vector)

Npt = length(x);
y = zeros(M,1);

% Making delay chain matrix
Zi = [x(Npt:-1:1)*M; Zi(1:end-Npt)];

if oflag == 1
    %Polyphase filtration
    D = reshape(Zi,M,length(Zi)/M);
    P = sum(D.*E_mx,2);    
    %DFT-modulation
%     y = (1/M)*fft(P);
    y = ifft(P);
end

%Updating output memory 
Zo=Zi;
end

