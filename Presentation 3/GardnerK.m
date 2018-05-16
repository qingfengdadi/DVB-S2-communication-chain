%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
clear; close all;

K_values = [0.025 0.05 0.1 0.5];
time_error = [];
%% Parameters
for m = 1:length(K_values)
K = K_values(m);
for iter = 1:100

Nbits = 5000; % bit stream length
f_cut = 1e6/2; % cut off frequency of the nyquist filter [Mhz]
M = 100; % oversampling factor (mettre à 100?)
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 2; % number of bits per symbol
modulation = 'pam'; % type of modulation 
bits_tx = randi(2,Nbits,1)-1;

tshift = 1;

%% Mapping
symbol_tx = mapping(bits_tx,Nbps,modulation);

%% Upsampling
symbol_tx_upsampled = upsample(symbol_tx,M);

%% Implementation of HHRC
RRCtaps = Nbps*M+1;
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
EbN0 = 10;

signal_power = (trapz(abs(signal_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/Nbits; % energy per bit

        
        N0 = Eb/10.^(EbN0/10);
        NoisePower = 2*N0*fsampling;
        noise = sqrt(NoisePower/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));
           
        signal_rx = signal_tx ;%+ noise;
        signal_rx = conv(signal_rx, h_time);
        symbol_rx_upsampled = signal_rx(RRCtaps:end-RRCtaps+1);
        
        %% Time shift
        symbol_rx_upsampled = symbol_rx_upsampled(1+tshift:end);
        
        %% Gardner
        L=length(symbol_rx_upsampled);
        L=L-mod(L,M);
        
        error=zeros(L/M,1);
        corr=zeros(L/M,1);
        
        prevHoho = symbol_rx_upsampled(1);
        
        for i=1:(L/M)-1
            a=((i-1)*M:M*i);
            b=symbol_rx_upsampled(1+(i-1)*M:i*M+1);
            c=M/2+(i-1)*M-error(i);
            c2=i*M-error(i);
            
            hihi = interp1(a,b,c,'pchip');
            hoho = interp1(a,b,c2,'pchip');
            
            corr(i)=(2*K)*real(hihi*(conj(hoho) - conj(prevHoho)));
            error(i+1) = error(i) + corr(i);
            prevHoho = hoho;
        end
        
        time_error(iter,:,m) = (tshift-error(1:end)).'*Tsymb;
        
end
end

%% Plot BER results
% load gardnerK.mat
time_error_mean = mean(time_error);
time_error_stdv = std(time_error);
mean1 = time_error_mean(1,:,1);
mean2 = time_error_mean(1,:,2);
mean3 = time_error_mean(1,:,3);
mean4 = time_error_mean(1,:,4);
stdv1 = time_error_stdv(1,:,1);
stdv2 = time_error_stdv(1,:,2);
stdv3 = time_error_stdv(1,:,3);
stdv4 = time_error_stdv(1,:,4);

figure
plot(mean1,'-r');hold on;
plot(mean2,'-g')
plot(mean3,'-b')
plot(mean4,'-c')
plot(mean1+stdv1,'--r')
plot(mean1-stdv1,'--r')
plot(mean2+stdv2,'--g')
plot(mean2-stdv2,'--g')
plot(mean3+stdv3,'--b')
plot(mean3-stdv3,'--b')
plot(mean4+stdv4,'--c')
plot(mean4-stdv4,'--c')

xlabel('Symbols');
ylabel('Time error (mean \pm stdv');
legend('K = 0.05','K = 0.2','K = 0.5','K = 1');
title('Convergence of the Gardner algorithm')
grid on;
