% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RANDOM BITS
N=1000; % bit stream length
bits = randi(2,1,N)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M=8; % oversampling factor
fs=2*f_cut; % symbol frequency
fsampling=M*fs; % sampling frequency
T=1/fs; % time between two symbols
Nbps=10; % number of bits per symbol
beta=0.3; % roll-off factor

%% mapping
signal = mapping(bits.',Nbps,'qam');

%% upsampling
signal_tx = upsample(signal,M);

%% implementation of transfer function
RRCtaps=65;
stepoffset=(1/RRCtaps)*fsampling;
highestfreq=(RRCtaps-1)*stepoffset/2;
f=linspace(-highestfreq,highestfreq,RRCtaps);
Hrrc=HRRC(f,T,beta);
h_freq=sqrt(Hrrc/max(Hrrc));
% figure;
% plot(f,h_freq);

h_time=fftshift(ifft(fftshift(h_freq)));
deltat=1/fsampling;
t=(-(RRCtaps-1)/2:(RRCtaps-1)/2)*deltat;

% figure;
% plot(t,h_time);
%% Convolution
signal_hrrc_tx = conv(signal_tx, h_time);
% figure;
% stem(signal_hrrc_tx);
%% Noise through the channel
signal_power=(trapz(abs(signal_hrrc_tx)).^2)*(1/fsampling); % total power
Eb=signal_power/N; % energy per bit
EbN0=5; % SNR (parameter)
N0=Eb/EbN0; 
NoisePower=2*N0*fsampling;
noise=sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx),1)+1i*randn(length(signal_hrrc_tx),1));

signal_rx = signal_hrrc_tx + noise;
signal_hhrc_rx = conv(signal_rx, h_time);

%% downsampling
signal_rx_down = downsample(signal_hhrc_rx, M);
%% demapping
bits_rx = demapping(signal_rx_down,Nbps,'qam');




