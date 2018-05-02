clear; close all;
addpath(genpath('Data'));

load BER_iterHardBPSK.mat; load EB_N0_iterHardBPSK.mat
BER_hard = BER; EbN0_hard = EbN0;
load BER_iterSoftBPSK.mat; load EB_N0_iterSoftBPSK.mat
BER_soft = BER; EbN0_soft = EbN0;

semilogy(EbN0_hard,BER_hard(:,1),'-',EbN0_hard,BER_hard(:,2),'-',EbN0_soft,BER_soft(:,2),'-'); hold on
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('Uncoded','Hard','Soft');
title('Hard vs Soft decoding for BPSK with 4 iterations')
grid on;