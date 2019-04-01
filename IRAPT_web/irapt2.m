function [F0_ref]=irapt2(Cfg,Signal,F0)
    
    Signal=resample(Signal,Cfg.fs_f0,Cfg.fs);
    
   
    [Amp,Frc,~]=warpingAnalyser(Signal, F0, Cfg.step_sub_smp, Cfg.fs_f0, ...
                            40, 'Sinc_hash_1000', 2:Cfg.max_harmonic_refine+1,3);
    
    nSamples=size(Amp,2);
                        
    Amp=[repmat(Cfg.initial_f0_value,nSamples,1) Amp'];
    Dev=repmat([1 1:Cfg.max_harmonic_refine],nSamples,1);

    Frc=[F0 Frc']./Dev;

    Valid=(abs(Frc-repmat(F0,1,Cfg.max_harmonic_refine+1))<Cfg.freq_tolerance_refine) & (Frc>0);   
    
    F0_ref=sum(Valid.*Amp.*Frc,2)./sum(Valid.*Amp,2);
end