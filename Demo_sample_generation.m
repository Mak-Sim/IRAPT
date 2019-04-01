%% Pattern generation

Skeleton=[0 0.05 0.1 0.25 0.5 0.75 0.8 0.9]*44100;
Skeleton=[Skeleton; 220 220 220 380 130 220 220 220];
Skeleton=[Skeleton; 220 220 220 400 80 220 220 220];
Pattern_smooth=interp1(Skeleton(1,:),Skeleton(2,:),0:Skeleton(1,end),'spline');
Pattern_sharp=interp1(Skeleton(1,:),Skeleton(3,:),0:Skeleton(1,end),'linear');
Pattern_prokladka=zeros(1,5000)+Skeleton(2,1);
Pattern=[Pattern_smooth zeros(1,2000)+Skeleton(2,1) Pattern_sharp];
Frc=Pattern_prokladka;

for Step=2:2:8
    Frc=[Frc Pattern(1:Step:end)];
    Frc=[Frc Pattern_prokladka];
end

%% Signal generation
Fs=44100;
N_Harm=10;
Ln=length(Frc);

Phs=cumsum(Frc)/Fs*2*pi;
Phs=Phs(:);
rand('state',0);
Phs_relative=randn(N_Harm,1);

Amps=(N_Harm:-1:1);
Amps=Amps/sum(Amps);

Sig=zeros(Ln,1);
for N=1:N_Harm
    Sig=Sig+Amps(N)*cos(Phs*N+Phs_relative(N));
end
% % plot(Sig);
wavwrite(Sig,Fs,16,'web_src/Demo');
True_frc=Frc(1:224:end)';

