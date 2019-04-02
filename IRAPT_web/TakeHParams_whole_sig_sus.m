function [Amps, Frc, Phs] = TakeHParams_whole_sig_sus(Sig, Cfg)
%TakeHParams_whole_sig -- takes harmonic parameters from signal.

M = 360;                % Number of channel
K = 4;                  % Number of merging subchannel
S = Cfg.step_sub_smp;   % Decimation coefficient
% m = 1;                  % Polyphase order
% N = m*M+1;            % Порядок фильтра-прототипа


% Original IRAPT filter-prototype
% % N=317;
% % h=fir1(N-1, 1/M, hamming(N), 'noscale');
% % h=[zeros(1,M/2+22) h zeros(1,M/2+22)];
% % N = length(h);

% Modified filter prototype for sustain phonation
N=720;
h=fir1(N, 1/M, hamming(N+1), 'noscale');
N = length(h);

% Filter bank Delay chain
Zi=zeros(N, 1);  
% Polyphase matrix
if mod(N,M)~=0 
    h = [h zeros(1,M-mod(N,M))];
    Zi=zeros(N+(M-mod(N,M)), 1);  
end

E_mx=reshape(h,M,ceil(N/M));

N_frames=ceil(length(Sig)/S);

% Zero padding of signal
Sig = [Sig zeros(1,ceil((N-1)/2))];

% Memory reserving
Amps = zeros(N_frames,(M/2)-(K-1));
Frc  = zeros(N_frames,(M/2)-(K-1));
Phs  = zeros(N_frames,(M/2)-(K-1));

offset = (N-1)/2 + 1 - S;
Frame=Sig(1:offset);
[~, Zi]=F0_Bank_get_samples(Frame', E_mx ,M, length(Frame), Zi);
Index_part1=repmat((1:K)',1,(M/2)-(K-1));
Index_part2=repmat((0:(M/2)-(K-1)-1),K,1);
ind_array = Index_part1+Index_part2;
ind_N=1;
for N=1:N_frames
    Frame=Sig((N-1)*S+offset+1:N*S+offset);
    [Amps(ind_N,:), Frc(ind_N,:), Phs(ind_N,:), Zi] = Get_harmonic_params(Frame', E_mx, M, S, Zi,ind_array);
    ind_N = ind_N + 1;
end

Frc = Frc/(2*pi)*Cfg.fs_f0;
end

