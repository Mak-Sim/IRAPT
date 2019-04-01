function [F0,Voc_value,Cfg]=irapt1(Signal)

%Global params
Cfg.fs=44100;
Cfg.fs_f0_target=6000;
Cfg.src_sub_ratio=round(Cfg.fs/Cfg.fs_f0_target);
Cfg.fs_f0=Cfg.fs/Cfg.src_sub_ratio;
Cfg.FD=35;                                               %Half band of the analysis filter
Cfg.max_harmonic_freq=14000;                              %Maximum harmonic frequency
Cfg.max_harmonic_number=100;                              %Maximum harmonic number

%Offset params
Cfg.step_sec=0.005;
Cfg.step_sub_smp=round(Cfg.step_sec*Cfg.fs_f0);
Cfg.step_smp=Cfg.step_sub_smp*Cfg.src_sub_ratio;

%Frame params
Cfg.frame_sec=0.05;
Cfg.frame_sub_smp=round(Cfg.frame_sec*Cfg.fs_f0/2)*2+1;  %F0 filter size size in samples (initial estimation)
Cfg.frame_smp=round(Cfg.frame_sec*Cfg.fs/2)*2+1;         %Harmonic analysis filter size in samples

%Chunk F0 params (initial estimation)
Cfg.chunk_f0_sec=0.3;
Cfg.chunk_f0_size=round(Cfg.chunk_f0_sec/Cfg.step_sec);
Cfg.chunk_f0_freqs=50:450;                               %Freq grid for DP
Cfg.f0_limits = [50 450];
Cfg.f0_max_step=23;                                      %Max F0 step in values of Cfg.chunk_f0_freqs  old value - 3
Cfg.f0_freq_lines=(Cfg.FD+Cfg.FD/4:Cfg.FD/2:Cfg.fs_f0/2-Cfg.FD)';   %Freq grid for analysis

%F0 params (whitening)
Cfg.whiten_frame_smp=Cfg.frame_sub_smp;
Cfg.whiten_step_smp=200;

%Voiced unvoiced decision parameters
Cfg.vu_freq_devision=3000;

%% Instanteniuos correlation estimation params
Cfg.corr_param.FFT_order=4096*4; % Must be even!
Cfg.corr_param.FFT_freq_line_size=floor(Cfg.corr_param.FFT_order/2)+1;

Interp_factor =2;
fs_new = Cfg.fs_f0*Interp_factor;
Cfg.corr_param.Interp_factor = Interp_factor;
Cfg.corr_param.Interp_filter_h_size = Interp_factor*12;
Cfg.corr_param.Interp_filter = fir1(Cfg.corr_param.Interp_filter_h_size*2,1/Interp_factor);

Cfg.corr_param.Left_index_actual = floor(Cfg.fs_f0./(Cfg.f0_limits(2)))+1;
Cfg.corr_param.Right_index_actual= ceil(Cfg.fs_f0./(Cfg.f0_limits(1)))+1;
Sinc_offset=ceil(Cfg.corr_param.Interp_filter_h_size/Interp_factor);

Cfg.corr_param.Left_index=Cfg.corr_param.Left_index_actual - Sinc_offset;
Cfg.corr_param.Right_index=Cfg.corr_param.Right_index_actual + Sinc_offset;
Cfg.corr_param.Actual_indeces = Sinc_offset*Interp_factor + 1: (Cfg.corr_param.Right_index_actual - Cfg.corr_param.Left_index)*Interp_factor+1;
Cfg.corr_param.Actual_freqs = fs_new./(Cfg.corr_param.Actual_indeces+(Cfg.corr_param.Left_index+Sinc_offset-1)*Interp_factor);

Cfg.corr_param.Actual_freqs = Cfg.fs_f0./linspace(Cfg.corr_param.Left_index_actual-1,Cfg.corr_param.Right_index_actual-1,(Cfg.corr_param.Right_index_actual-Cfg.corr_param.Left_index_actual)*Interp_factor+1);
Cfg.corr_param.Actual_freqs_num = length(Cfg.corr_param.Actual_freqs);

n_point=diff(Cfg.f0_limits)+1;

Cfg.corr_param.Window=(0:n_point-1)/(n_point-1)*0.25+0.75;
Cfg.corr_param.Window=MyFit(Cfg.f0_limits(1):Cfg.f0_limits(2),Cfg.corr_param.Window,Cfg.corr_param.Actual_freqs(end:-1:1));
Cfg.corr_param.Window=Cfg.corr_param.Window(end:-1:1);

%% Dynamic programming params

Index_part1=repmat((1:Cfg.f0_max_step*2+1)',1,Cfg.corr_param.Actual_freqs_num);
Index_part2=repmat((0:Cfg.corr_param.Actual_freqs_num-1),Cfg.f0_max_step*2+1,1);
Cfg.dp.Ind_array = Index_part1+Index_part2;
Cfg.dp.Leakage_factor=0.95;

%F0 params (final estimation)
Cfg.f0_final_deviation=0.1; % Max allowed deviation from F0_crude determined by f0_freq_lines (in fractional parts)

%F0 params (refine)
Cfg.max_harmonic_refine=8;
Cfg.freq_tolerance_refine=Cfg.FD;
Cfg.initial_f0_value=0.01;

%Initialization
Signal=resample(Signal,Cfg.fs_f0,Cfg.fs);
Signal_len=length(Signal);
Addon_left_size=floor(Cfg.frame_sub_smp/2);
Addon_right_size=(Cfg.chunk_f0_size-1)*Cfg.step_sub_smp+Addon_left_size-(Signal_len-1-ceil(Signal_len/Cfg.step_sub_smp-1)*Cfg.step_sub_smp);
Signal=[Signal zeros(1,Addon_right_size)];
Signal_len=Signal_len+Addon_right_size;

Samples=1:Cfg.step_sub_smp:Signal_len-Addon_left_size;
Smp_number=length(Samples);
Smp_valued_number=Smp_number-Cfg.chunk_f0_size+1;

Chunk_n=Cfg.chunk_f0_size;

[Amp_bank,Frc_bank,~]=TakeHParams_whole_sig(Signal,Cfg);

%Chunk initialization
Value_matrix=ValueFunc_corr_line_FFT_interp(Amp_bank(1:Chunk_n,:)',Frc_bank(1:Chunk_n,:)',Cfg);

[~,q,D_vec,phi]=dp_pitch_full(Value_matrix,Cfg.f0_max_step);

%Frame by frame processing (initial F0 estimation)
F0=zeros(Smp_valued_number,1);
Voc_value=zeros(Smp_valued_number,1);

F0(1)=Cfg.corr_param.Actual_freqs(q(1));
En=sum(Amp_bank(1:Chunk_n,:).^2,2);
En(En<10^-4)=10^-4;
Voc_value(1:Cfg.chunk_f0_size)=-min(Value_matrix(:,:),[],2)./En;

for N=2:Smp_valued_number %Consider minimization of D_vec

    
    Amp=Amp_bank(Chunk_n+N-1,:)';
    Frc=Frc_bank(Chunk_n+N-1,:)';
    
    Value_vec=ValueFunc_corr_line_FFT_interp(Amp,Frc,Cfg);
    [q,D_vec,phi]=dp_pitch_step(Value_vec,D_vec,phi, Cfg);
    F0(N)=Cfg.corr_param.Actual_freqs(q);
    F0(N)=Get_final_F0(Amp_bank(N,:),Frc_bank(N,:),F0(N),Cfg);
    
    En=sum(Amp.^2);
    En(En<10^-4)=10^-4;
    Voc_value(N+Cfg.chunk_f0_size-1)=-min(Value_vec)./En;
end

end

function [F0_final]=Get_final_F0(Amp_vec,Frc_vec,F0_crude,Cfg)

    Frc_vec_frac=Frc_vec/F0_crude;
    Frc_vec_h_numbers=round(Frc_vec_frac);
    Frc_vec_frac=abs(Frc_vec_frac-Frc_vec_h_numbers);

    Inds_valued=(Frc_vec_frac<=Cfg.f0_final_deviation & Frc_vec_h_numbers>0);
    
    if(isempty(Inds_valued))
        F0_final=F0_crude;
    else
        Amp_vec=Amp_vec(Inds_valued);
        Frc_vec=Frc_vec(Inds_valued);
        Frc_vec_h_numbers=Frc_vec_h_numbers(Inds_valued);
        F0_final=sum(Amp_vec.*(Frc_vec./Frc_vec_h_numbers))/sum(Amp_vec);
    end
end

function [q,D_vec,phi]=dp_pitch_step(Value_vec,D_vec,phi,Cfg)

D_vec=D_vec*Cfg.dp.Leakage_factor;

r =size(phi,1);
c = Cfg.corr_param.Actual_freqs_num;
Max_step = Cfg.f0_max_step;

phi(1:r-1,:) = phi(2:end,:);

D_vec_new=D_vec;
D_vec_new(Max_step+1:end-Max_step)=Value_vec;

[Min_val,Min_ind]=min(D_vec(Cfg.dp.Ind_array),[],1);
D_vec_new(Max_step+1:Max_step+c)=D_vec_new(Max_step+1:Max_step+c)+Min_val;
phi(r,1:c)=Min_ind;

% Traceback from top
[~,j] = min(D_vec_new(Max_step+1:Max_step+c));

for R=r:-1:2
    tb=phi(R,j);
    j=j-(Max_step+1)+tb;
end

q=j;

D_vec = D_vec_new;

end

function [p,q,D_vec,phi]=dp_pitch_full(Value_matrix,Max_step)

[r,c] = size(Value_matrix);

D = zeros(r+1, c+Max_step*2);
D(1,:) = 0;
D(:,[1:Max_step Max_step+c+1:c+Max_step*2]) = NaN;
D(2:(r+1), Max_step+1:Max_step+c) = Value_matrix;

phi = zeros(r,c); %Return indexes

for i = 1:r; 
  for j = Max_step+1:Max_step+c;
    [dmax, tb] = min(D(i, j-Max_step:j+Max_step));
    D(i+1,j) = D(i+1,j)+dmax;
    phi(i,j-Max_step) = tb;
  end
end

% Traceback from top
i = r; 
[V,j] = min(D(r+1,Max_step+1:Max_step+c));
p = i;
q = j;
while i > 1 && j >= 1
  tb = phi(i,j);

  i=i-1;
  j=j-(Max_step+1)+tb;
  
  if(j==0)
      break;
  end
  
  p = [i,p];
  q = [j,q];
end

D_vec = D(r+1,:);
end