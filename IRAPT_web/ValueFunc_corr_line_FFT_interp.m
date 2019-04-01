function [Values]=ValueFunc_corr_line_FFT_interp(Amp,Frc,Cfg) % Amp, Frc, Freq - params x frames
[~,m]=size(Frc);
n_point=length(Cfg.corr_param.Actual_freqs);

Values=zeros(n_point,m);

Amp(Frc<0 | Frc>Cfg.fs_f0/2)=0;
Frc(Frc<0 | Frc>Cfg.fs_f0/2)=0;

FFT_order           =Cfg.corr_param.FFT_order; 
FFT_freq_line_size  =Cfg.corr_param.FFT_freq_line_size;

Interp_factor=Cfg.corr_param.Interp_factor;
Interp_filter_h_size=Cfg.corr_param.Interp_filter_h_size;
Interp_filter=Cfg.corr_param.Interp_filter;


for M=1:m
    Amp_vec=Amp(:,M)';
    Frc_vec=Frc(:,M)';

    if(sum(Amp_vec>0)<=2)
        continue;
    end
    
    FFT_inds=round(Frc_vec/(Cfg.fs_f0/FFT_order)+1);
    FFT_amps=zeros(FFT_freq_line_size,1);

    FFT_amps(FFT_inds)=Amp_vec.^2;
    FFT_amps=[FFT_amps' FFT_amps(end-1:-1:2,:)'];

    Corr_sig=ifft(FFT_amps)*FFT_order/2;

    Corr_sig=Corr_sig(Cfg.corr_param.Left_index:Cfg.corr_param.Right_index);

    Corr_sig_interp=zeros((Cfg.corr_param.Right_index-Cfg.corr_param.Left_index+1)*Interp_factor,1);
    
    Corr_sig_interp(1:Interp_factor:end)=Corr_sig*Interp_factor;
    
    Corr_sig_interp=Filter_no_offset(Interp_filter,Corr_sig_interp,Interp_filter_h_size,0);    
    
    Values(:,M)=Corr_sig_interp(Cfg.corr_param.Actual_indeces);
     
end

Values=-Values'.*repmat(Cfg.corr_param.Window,m,1);

