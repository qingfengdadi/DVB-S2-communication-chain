clear; close all;
addpath(genpath('Data'));
addpath(genpath('Rapport'));

% load BER_iterHardBPSK.mat; load EB_N0_iterHardBPSK.mat
% BER_hard = BER; EbN0_hard = EbN0;
% load BER_iterSoftBPSK.mat; load EB_N0_iterSoftBPSK.mat
% BER_soft = BER; EbN0_soft = EbN0;

% semilogy(EbN0_hard,BER_hard(:,1),'-',EbN0_hard,BER_hard(:,2),'-',EbN0_soft,BER_soft(:,2),'-'); hold on
% xlabel('E_B/N_0 [dB]');
% ylabel('BER');
% legend('Uncoded','Hard','Soft');
% title('Hard vs Soft decoding for BPSK with 4 iterations')
% grid on;

%% Plusieurs constellation
load hard_4qam.mat;
BER4 = BER; EbN04 = EbN0;
load hard_16qam.mat;
BER16= BER; EbN016 = EbN0;
load hard_64qam.mat;
BER64 = BER; EbN064 = EbN0;

semilogy(EbN04,BER4(:,1),'r--',EbN04,BER4(:,2),'r-',EbN016,BER16(:,1),'g--',EbN016,BER16(:,2),'g-',EbN064,BER64(:,1),'b--',EbN064,BER64(:,2),'b-'); ylim([1e-3,1e0]);hold on
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('4QAM uncoded','4QAM coded','16QAM uncoded','16QAM coded','64QAM uncoded','64QAM coded');
title('Hard decoding')
grid on;
