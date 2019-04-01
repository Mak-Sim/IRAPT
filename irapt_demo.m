% IRAPT demo

addpath('IRAPT_web');

[Sig,Fs]=audioread('web_src/Demo.wav');


[F0, time_marks] = irapt(Sig, Fs, 'irapt1');

figure; hold on;
plot(time_marks,F0,'LineWidth',1.5);
load('web_src/Demo_true_F0');
plot(time_marks,True_F0,'-.g','LineWidth',1.5);
legend('Estimated F0','True F0');
xlabel('Time, sec','FontSize',12);
ylabel('Frequency, Hz','FontSize',12);
grid on;
