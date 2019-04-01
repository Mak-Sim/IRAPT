function [Amp_out, Frc_out, Phs_out] = warpingAnalyser(Sig, FN, Time_step, fs, ...
                            Phase_stages, hashTableName, actualChannels,numberCycles)
%WARPINGANALYSER -- (super) special speech analyser
% Sig -- signal
% FN  -- frequency contour
% fs  -- sampling frequency
% hashTableName -- Sinc-functions table with different shifts

half_size=100;              % Interpolation window size
Ln = length(Sig);
Ph_cum=0;
Frc=FN;
Inds_temp=1:Time_step:Ln;
Frc=MyFit(Inds_temp,Frc,1:Ln);

%% Filter bank preparation
nChannel = Phase_stages;    % Number of channel of filter bank
N = Phase_stages*2*numberCycles+1;     % (always odd)

h=fir1(N-1,1/nChannel/4,kaiser(N,9),'noscale');

% Initial delay-chain values
Zi=zeros(N, 1);  

% Making polyphase matrix
if mod(N, nChannel)~=0 
    h = [h zeros(1,nChannel-mod(N,nChannel))];
    Zi=zeros(N+(nChannel-mod(N,nChannel)), 1);  
end

E_mx=reshape(h,nChannel,ceil(N/nChannel));
%% Analyser parameters
N_chanels=length(actualChannels);
Phase_step = 2*pi/Phase_stages;
Current_sample=Time_step;
Latency=floor(N/2);  %Half impulse response of the filter bank (in phase samples) -1
Latency_timer=[];
Pre_mark=0;
Post_mark_integer=0;
Time_marks=[];
Flag_get_params_event=false;

Ln_interp=length(Inds_temp);

Amp_out = zeros(N_chanels, Ln_interp);
Phs_out = zeros(N_chanels, Ln_interp);
Frc_out = zeros(N_chanels, Ln_interp);

Frame_ind=1;

Test_moments1=[];
Test_moments2=[];
Test_ind=1;

load(hashTableName);
N_sinc_hash=size(Sinc_hash,1)-1;

for N=1:Ln
    Ph_last=Ph_cum;
    Ph_cum=Ph_cum+Frc(N)/fs*2*pi;
    P_points=floor(Ph_cum/Phase_step);
    
    if(Current_sample==Time_step)
        Latency_timer=[Latency_timer Latency];
        Current_sample=0;
        Flag_get_params_event=true;
        Post_mark_integer=0;
    end
    
    if(P_points>0)
        Moments=MyFit([Ph_last Ph_cum],[0 1],Phase_step:Phase_step:Phase_step*P_points);
        
        Left_ind=N-half_size;
        Right_ind=Left_ind+half_size*2+1;
        
        Left_zeros=[];
        Right_zeros=[];
        
        if(Left_ind<1)
            Left_zeros=zeros(1,-Left_ind+1);
            Left_ind=1;
        end
        
        if(Right_ind>Ln)
            Right_zeros=zeros(1,Right_ind-Ln);
            Right_ind=Ln;
        end
        
        Frame=[Left_zeros Sig(Left_ind:Right_ind) Right_zeros];
        
        if(Flag_get_params_event)
            Time_marks=[Time_marks [Pre_mark;Moments(1)+Post_mark_integer]];
            Test_moments1=[Test_moments1 Moments(1)+Post_mark_integer-Pre_mark];
            Test_ind=Test_ind+1;
            
            Flag_get_params_event=false;
        end
        
        for M=1:length(Moments)
            sigNew = sum(Sinc_hash(round(Moments(M)*N_sinc_hash)+1,:).*Frame(2:end-1));
            if(isempty(Latency_timer))
                % Send sample away
                [~,Zi]=dft_fb(sigNew, E_mx, nChannel, 0, Zi);  
            elseif(Latency_timer(1)>1)
                % The same
                [~,Zi]=dft_fb(sigNew, E_mx, nChannel, 0, Zi);  
                Latency_timer=Latency_timer-1;
            elseif(Latency_timer(1)==1)
                % Send sample away and get the first pack in return
                [Pre_params,Zi]=dft_fb(sigNew, E_mx, nChannel, 1, Zi);  
                Latency_timer=Latency_timer-1;
            else                
                % Send sample away and get second pack in return and calculate
                % interpolation
                [Post_params, Zi] = dft_fb(sigNew, E_mx, nChannel, 1, Zi);
                
                % Harmonic params calculation
                % Amplitude
                k_order=0;
                for k=actualChannels
                    k_order=k_order+1;
                    Amp_out(k_order,Frame_ind)=2*MyFit(Time_marks(:,1)',...
                        [abs(Pre_params(k)) abs(Post_params(k))],0);
                    % Phase
                    Phs_out(k_order,Frame_ind) = MyFit(Time_marks(:,1)',...
                        My_unwrap([angle(Pre_params(k)) angle(Post_params(k))]),0);
                    % Frequency
                    Frc_out(k_order,Frame_ind)=(diff(My_unwrap([angle(Pre_params(k)) ...
                        angle(Post_params(k))]))')/(Time_marks(2,1)-Time_marks(1,1))/2/pi*fs;
                    
                end
                
                Test_moments2=[Test_moments2 Time_marks(2,1)-Time_marks(1,1)];
                Frame_ind=Frame_ind+1;
                
                Latency_timer=Latency_timer-1;
                Latency_timer(1)=[];
                Time_marks(:,1)=[];
            end                                        

        end
        Pre_mark=Moments(end)-1;
        
        Ph_cum=Ph_cum-P_points*Phase_step;
    else
        if(~Flag_get_params_event)
            Pre_mark=Pre_mark-1;
        else
            Post_mark_integer=Post_mark_integer+1;
        end
    end
        
    Current_sample=Current_sample+1;
end

% Tail processing
N_tail=max(Latency_timer);

for N=1:N_tail
    if(isempty(Latency_timer))
        % Send sample away
        [~,Zi]=dft_fb(0, E_mx, nChannel, 0, Zi);
    elseif(Latency_timer(1)>1)
        % The same
        [~,Zi]=dft_fb(0, E_mx, nChannel, 0, Zi);
        Latency_timer=Latency_timer-1;
    elseif(Latency_timer(1)==1)
        % Send sample away and get the first pack in return
        [Pre_params,Zi]=dft_fb(0, E_mx, nChannel, 1, Zi);
        Latency_timer=Latency_timer-1;
    else
        % Send sample away and get second pack in return and calculate
        % interpolation
        [Post_params,Zi]=dft_fb(0, E_mx, nChannel, 1, Zi);
        
        % Harmonic params calculation
        % Amplitude
        k_order=0;
        for k=actualChannels
            k_order=k_order+1;
            Amp_out(k_order,Frame_ind)=2*MyFit(Time_marks(:,1)',...
                [abs(Pre_params(k)) abs(Post_params(k))],0);
            % Phase
            Phs_out(k_order,Frame_ind) = MyFit(Time_marks(:,1)',...
                My_unwrap([angle(Pre_params(k)) angle(Post_params(k))]),0);
            % Frequency
            Frc_out(k_order,Frame_ind)=(diff(My_unwrap([angle(Pre_params(k)) ...
                angle(Post_params(k))]))')/(Time_marks(2,1)-Time_marks(1,1))/2/pi*fs;
            
        end
        
        Frame_ind=Frame_ind+1;
        
        Latency_timer=Latency_timer-1;
        Latency_timer(1)=[];
        Time_marks(:,1)=[];
    end
end
end

