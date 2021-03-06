% IRAPT demo

% Cite as
% @inproceedings{Azarov-12, 
% 	author={Elias Azarov and Maxim Vashkevich and Alexander A. Petrovsky}, 
% 	booktitle={Proc. of the 20th European Signal Processing Conference (EUSIPCO)}, 
% 	title={Instantaneous pitch estimation based on {RAPT} framework}, 
% 	year={2012}, 
% 	pages={2787-2791}, 	
% 	month={Aug}}


addpath('IRAPT_web');

[Sig,Fs]=audioread('web_src/Demo.wav');

[F0, ~,time_marks] = irapt(Sig, Fs, 'irapt1','speech');

figure; hold on;
plot(time_marks,F0,'LineWidth',1.5);
load('web_src/Demo_true_F0');
plot(time_marks,True_F0,'-.g','LineWidth',1.5);
legend('Estimated F0','True F0');
xlabel('Time, sec','FontSize',12);
ylabel('Frequency, Hz','FontSize',12);
grid on;

%%
[Sig,Fs]=audioread('web_src/001.wav');

[F0_sp, voc_sp,t_mark_sp] = irapt(Sig, Fs, 'irapt2','speech');
[F0_sus, voc_sus,t_mark_sus] = irapt(Sig, Fs, 'irapt2','sustain phonation');

figure;
subplot(311)
plot((0:length(Sig)-1)/Fs,Sig);
xlabel('Time, sec','FontSize',12);

title('Pathological voice (example)');

subplot(312)
plot(t_mark_sp,F0_sp,'LineWidth',1.5);
hold on;
plot(t_mark_sus,F0_sus,'-.r','LineWidth',1.5);
legend('IRAPT (speech)','IRAPT (sustain phonation)');
xlabel('Time, sec','FontSize',12);
ylabel('Frequency, Hz','FontSize',12);
grid on;

subplot(313)
plot(t_mark_sp,voc_sp,'LineWidth',1.5); hold on;
plot(t_mark_sus,voc_sus,'-.r','LineWidth',1.5);


legend('IRAPT (speech)','IRAPT (sustain phonation)');
xlabel('Time, sec','FontSize',12);
ylabel('voiced/unvoiced','FontSize',12);
grid on;
