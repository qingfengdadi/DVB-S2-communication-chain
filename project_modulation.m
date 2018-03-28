%% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));
clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RANDOM BITS
Nbits = 30000; % bit stream length
bits_tx = randi(2,Nbits,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 8; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols
Nbps = 4; % number of bits per symbol
modulation = 'qam'; % type of modulation 
beta = 0.3; % roll-off factor

%% mapping
signal = mapping(bits_tx,Nbps,modulation);

%% upsampling
signal_tx = upsample(signal,M);

%% implementation of transfer function
RRCtaps = 165;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrrc = HRRC(f,Tsymb,beta);
h_t = ifft(Hrrc);
h_freq = sqrt(fft(h_t/max(h_t)));

% h_freq = sqrt(fft(ifft(Hrrc)/max(ifft(Hrrc))));
% h_freq = (fft(ifft(Hrrc)/max(ifft(Hrrc))));

% figure;
% plot(f,h_freq);

h_time = fftshift(ifft(ifftshift(h_freq)));
deltat = 1/fsampling;
t = (-(RRCtaps-1)/2:(RRCtaps-1)/2)*deltat;

%%plot all shifted signals
% figure;
% plot(t,h_time);
% hold on
% plot(t+Tsymb,h_time);
% hold on
% plot(t+2*Tsymb,h_time);
% hold on
% plot(t+3*Tsymb,h_time);
% hold on
% plot(t+4*Tsymb,h_time);
% grid on;

%% Convolution
signal_hrrc_tx = conv(signal_tx, h_time);
% figure;
% stem(signal_hrrc_tx);

%% Noise through the channel
<<<<<<< HEAD
EbN0 = -5:44;
BER = zeros(length(EbN0),1);
signal_power = (trapz(abs(signal_hrrc_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/Nbits; % energy per bit

for j = 1:length(EbN0)
%     EbN0 = 1000; % SNR (parameter)
    N0 = Eb/10.^(EbN0(j)/10);
    NoisePower = 2*N0*fsampling;
    noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx),1)+1i*randn(length(signal_hrrc_tx),1));

    signal_rx = signal_hrrc_tx + noise;
    signal_hhrc_rx = conv(signal_rx, h_time);

    signal_hhrc_rx_trunc = signal_hhrc_rx(RRCtaps:end-RRCtaps+1);
    %% downsampling
    signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

    %% demapping
    bits_rx = demapping(signal_rx_down,Nbps,modulation);
    
    %% BER
    BER(j) = length(find(bits_tx ~= bits_rx))/length(bits_tx);
end

semilogy(EbN0,BER)
grid on

% if bits_rx == bits_tx
%     disp('ok')
% end
=======
signal_power = (trapz(abs(signal_hrrc_tx)).^2)*(1/fsampling); % total power
Eb = signal_power/N; % energy per bit
% EbN0 = [0.01 0.1 1 10 100 1000 1e+4 1e+5 1e+6 1e+7 1e+8 1e+9 1e+10 1e+11 1e+12 1e+13 1e+14]; % SNR (parameter)
EbN0=1e+15;
BER=[];
for i=1:length(EbN0)
    N0 = Eb/EbN0(i); 
    NoisePower = 2*N0*fsampling;
    noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx),1)+1i*randn(length(signal_hrrc_tx),1));

    signal_rx = signal_hrrc_tx+noise;
    signal_hhrc_rx = conv(signal_rx, h_time);

    signal_hhrc_rx_trunc = signal_hhrc_rx(RRCtaps:end-RRCtaps+1);
    
    %% downsampling
    signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

    %% demapping
    bits_rx = demapping(signal_rx_down,Nbps,modulation);

    if bits_rx ~= bits_tx
        disp('ha...')
    end
    BER(i) = (N-nnz(bits_rx == bits_tx))/N;
end

% plot(20*log10(EbN0),log10(BER))
>>>>>>> 157433bf18f8805880daa43b486cf644f1aab975
