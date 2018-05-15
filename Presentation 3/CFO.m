%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
addpath(genpath('Data'));
clear; close all;

%% Parameters
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 8; % oversampling factor (mettre à 100?)
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 6; % number of bits per symbol
modulation = 'qam'; % type of modulation 

Nbits = 30000; % bit stream length
bits_tx = randi(2,Nbits,1)-1;

fc = 2e+9;
ppm = fc*1e-6;
CFO_values = [0 10 40 70]*ppm;
phi0 = 0;

%% Mapping of encoded signal
symbol_tx = mapping(bits_tx,Nbps,modulation);

%% Upsampling
symbol_tx_upsampled = upsample(symbol_tx,M);

%% Implementation of HHRC
RRCtaps = 365;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrc = HRC(f,Tsymb,beta);
h_t = ifft(Hrc);
h_freq = sqrt(fft(h_t/max(h_t)));
h_time = fftshift(ifft(ifftshift(h_freq)));

%% Convolution
signal_tx = conv(symbol_tx_upsampled, h_time);

%% Noise through the channel coded
EbN0 = 0:1:50;
BER = zeros(length(EbN0),4);
scatterData = zeros(length(symbol_tx),4);

signal_power = (trapz(abs(signal_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/Nbits; % energy per bit

for m = 1:length(CFO_values)
    cfo = CFO_values(m);
    for j = 1:length(EbN0)
        N0 = Eb/10.^(EbN0(j)/10);
        NoisePower = 2*N0*fsampling;
        noise = sqrt(NoisePower/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));
    
        exp_cfo1 = exp(1j*(2*pi*cfo*((0:length(signal_tx)-1)-(RRCtaps-1)/2)*ts+phi0))';
           
        signal_rx = signal_tx + noise;
        signal_rx = signal_rx.*exp_cfo1;
        signal_rx = conv(signal_rx, h_time);
        symbol_rx_upsampled = signal_rx(RRCtaps:end-RRCtaps+1);

        %% Downsampling
        symbol_rx = downsample(symbol_rx_upsampled, M);
        
        exp_cfo2 = exp(-1j*(2*pi*cfo*(0:length(symbol_rx)-1)*ts*M))'; % compensation to only observe ISI
        symbol_rx = symbol_rx.*exp_cfo2;
        
        %% Demapping
        bits_rx = (demapping(symbol_rx,Nbps,modulation))';
        BER(j,m) = length(find(bits_tx ~= bits_rx'))/length(bits_rx');
        
    end
    scatterData(:,m) = symbol_rx;
end

%% Plot BER results
% load phase_drift.mat
% load ISI.mat
figure
semilogy(EbN0,BER(:,1),'-',EbN0,BER(:,2),'-o',EbN0,BER(:,3),'-o',EbN0,BER(:,4),'-o');
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('CFO = 0 ppm','CFO = 10 ppm','CFO = 40 ppm','CFO = 70 ppm');
title('BER degradation for increasing values of the CFO')
grid on;

%% Plot Constellation results for SNR = 20 
scatterplot(symbol_tx,1,0,'r.')          
title('TX Symbols')
grid on

scatterplot(scatterData(:,1),1,0,'r.')          
title('CFO = 0 ppm')
grid on
 
scatterplot(scatterData(:,2),1,0,'r.')     
title('CFO = 10 ppm')
grid on
    
scatterplot(scatterData(:,3),1,0,'r.')         
title('CFO = 40 ppm')
grid on
      
scatterplot(scatterData(:,4),1,0,'r.')         
title('CFO = 70 ppm')
grid on
