%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
clear; close all;

%% Parameters
Nbits = 30000; % bit stream length
f_cut = 1e6; % cut off frequency of the nyquist filter [Mhz]
M = 100; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 4; % number of bits per symbol
modulation = 'qam'; % type of modulation 
bits_tx = randi(2,Nbits,1)-1;

tshift_values = [0 2 5 10];

%% Mapping
symbol_tx = mapping(bits_tx,Nbps,modulation);

%% Upsampling
symbol_tx_upsampled = upsample(symbol_tx,M);

%% Implementation of HHRC
RRCtaps = 365;%Nbps*M+101;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrc = HRC(f,Tsymb,beta);
h_t = ifft(Hrc);
h_freq = sqrt(fft(h_t/max(h_t)));
h_time = fftshift(ifft(ifftshift(h_freq)));

%% Convolution
signal_tx = conv(symbol_tx_upsampled, h_time);

%% Noise through the channel
EbN0 = -5:1:50;
BER = zeros(length(EbN0),4);
scatterData = zeros(length(symbol_tx),4);

signal_power = (trapz(abs(signal_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/Nbits; % energy per bit

for m = 1:length(tshift_values)
    for j = 1:length(EbN0)
        N0 = Eb/10.^(EbN0(j)/10);
        NoisePower = 2*N0*fsampling;
        noise = sqrt(NoisePower/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));
           
        signal_rx = signal_tx + noise;
        signal_rx = conv(signal_rx, h_time);
        symbol_rx_upsampled = signal_rx(RRCtaps:end-RRCtaps+1);
        
        %% Time shift
        symbol_rx_upsampled = symbol_rx_upsampled(1+tshift_values(m):end);

        %% Downsampling
        symbol_rx = downsample(symbol_rx_upsampled, M);
        
        %% Demapping
        bits_rx = (demapping(symbol_rx,Nbps,modulation))';
        BER(j,m) = length(find(bits_tx ~= bits_rx'))/length(bits_rx');
        
    end
    scatterData(:,m) = symbol_rx;
end

%% Plot BER results
load time_shift.mat
figure
semilogy(EbN0,BER(:,1),'-',EbN0,BER(:,2),'-o',EbN0,BER(:,3),'-o',EbN0,BER(:,4),'-o');
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('t_0 = 0 Ts','t_0 = 0.02 Ts','t_0 = 0.05 Ts','t_0 = 0.1 Ts');
title('BER degradation for increasing values of the sample time shift')
grid on;

% %% Plot Constellation results for SNR = 20 
% scatterplot(scatterData(:,1),1,0,'r.')          
% title('Time Shift t_0 = 0')
% grid on
%  
% scatterplot(scatterData(:,2),1,0,'r.')     
% title('Time Shift t_0 = 2')
% grid on
%     
% scatterplot(scatterData(:,3),1,0,'r.')         
% title('Time Shift t_0 = 5')
% grid on
%       
% scatterplot(scatterData(:,4),1,0,'r.')         
% title('Time Shift t_0 = 10')
% grid on
