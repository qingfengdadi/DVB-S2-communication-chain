%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
clear; close all;

N_values = [20,30,40];
K_values = [1,8,16];
%% Parameters
for k = 1:length(K_values)
for n = 1:length(N_values) 
K = K_values(k);
N = N_values(n);

for iter = 1:2500
    
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 8; % oversampling factor (mettre à 100?)
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 2; % number of bits per symbol
modulation = 'qam'; % type of modulation 

Nbits = 200; % data bit stream length
Npilot = N; % Number of pilot bits
bits_pilot = randi(2,1,Npilot)-1;
bits_data = randi(2,1,Nbits)-1;

fc = 2e+9;
ppm = fc*1e-6;
CFO_values = [2]*ppm;
phi0 = 0;

pilot_pos = 1;

%% Mapping of encoded signal
symbol_pilot = mapping(bits_pilot',Nbps,modulation);
symbol_data = mapping(bits_data',Nbps,modulation);
symbol_tx = [symbol_data(1:pilot_pos-1);symbol_pilot;symbol_data(pilot_pos:end)];

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
EbN0 = 0:16;
BER = zeros(length(EbN0),2);
scatterData = zeros(length(symbol_data),2);

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
        symbol_rx = signal_rx(RRCtaps:end-RRCtaps+1);
        
        %% Downsampling
        symbol_rx = downsample(symbol_rx, M);
%         symbol_rx = [symbol_rx(1:est_n/M);symbol_data(est_n/M:end)];

        %% Frame and frequency acquisition
        [est_n,est_cfo] = cfoEstimate(symbol_rx,symbol_pilot,Tsymb,K);
        
        exp_cfo2 = exp(-1j*(2*pi*cfo*(0:length(symbol_rx)-1)*ts*M))';
        symbol_rx = symbol_rx.*exp_cfo2;
        
        %% Demapping
        bits_rx = (demapping(symbol_rx,Nbps,modulation))';
%         BER(j,m) = length(find(bits_data ~= bits_rx))/length(bits_rx');
        time_error(iter,j,k,n) = est_n - pilot_pos;
        freq_error(iter,j,k,n) = cfo + est_cfo;
    end
%     scatterData(:,m) = symbol_rx;
end
end
end
end

%% Plot results
% load CFO_K_N.mat
time_error_mean = std(time_error);  
freq_error_mean = std(freq_error);

N_values = [20, 30, 40];
K_values = [1, 8, 16];

figure
plot(EbN0,time_error_mean(1,:,2,1),'-r');hold on;
plot(EbN0,time_error_mean(1,:,2,2),'-g');
plot(EbN0,time_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('N = 20, K = 8','N = 30, K = 8','N = 40, K = 8');
title('Time error variances')
grid on;

figure
plot(EbN0,time_error_mean(1,:,1,2),'-r');hold on;
plot(EbN0,time_error_mean(1,:,2,2),'-g');
plot(EbN0,time_error_mean(1,:,3,2),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('N = 30, K = 1','N = 30, K = 8','N = 30, K = 16');
title('Time error variances')
grid on;

figure
plot(EbN0,freq_error_mean(1,:,2,1),'-r');hold on;
plot(EbN0,freq_error_mean(1,:,2,2),'-g');
plot(EbN0,freq_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 20, K = 8','N = 30, K = 8','N = 40, K = 8');
title('Frequency error variances')
grid on;

figure
plot(EbN0,freq_error_mean(1,:,1,2),'-r');hold on;
plot(EbN0,freq_error_mean(1,:,2,2),'-g');
plot(EbN0,freq_error_mean(1,:,3,2),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 30, K = 1','N = 30, K = 8','N = 30, K = 16');
title('Frequency error variances')
grid on;