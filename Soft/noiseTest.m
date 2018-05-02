%% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRRC'));
clear; close all;

%% Parameters
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 10; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols
Nbps = 1; % number of bits per symbol
modulation = 'pam'; % type of modulation 
beta = 0.3; % roll-off factor

% Npackets = 60; % choose such that Nbits/Nbps is an integer
% packetLength = 128;
% codedWordLength = 2*packetLength;
% Nbits = Npackets*packetLength; % bit stream length
% NcodedBits = Npackets*codedWordLength; % full coded word length
% bits_tx = randi(2,Nbits,1)-1;
% bits_tx_coded = zeros(NcodedBits,1);

%% LDPC encoder
% H0 = makeLdpc(packetLength, codedWordLength, 0, 1, 3);
% for k=1:Npackets
%     packet_tx = bits_tx(1+(k-1)*packetLength : k*packetLength);
%     [codedbits, H] = makeParityChk(packet_tx , H0, 0);
%     bits_tx_coded(1+(k-1)*codedWordLength : k*codedWordLength) = [codedbits packet_tx];
% end
% tannerGraph = buildTannerGraph(H);

%% Mapping of encoded signal
H = [0 1 0 1 1 0 0 1;
     1 1 1 0 0 1 0 0;
     0 0 1 0 0 1 1 1;
     1 0 0 1 1 0 1 0];
graph = buildTannerGraph(H);
bits_rx = [1 1 0 1 0 1 0 1];
bits_tx_coded = [1 0 0 1 0 1 0 1]';
% signal = mapping(bits_tx,Nbps,modulation);
signal_coded = mapping(bits_tx_coded,Nbps,modulation);

%% Upsampling
signal_tx = upsample(signal_coded,M);

%% Implementation of HHRC
RRCtaps = 365;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrrc = HRRC(f,Tsymb,beta);
h_t = ifft(Hrrc);
h_freq = sqrt(fft(h_t/max(h_t)));
h_time = fftshift(ifft(ifftshift(h_freq)));

%% Plot HHRC
% deltat = 1/fsampling;
% t = (-(RRCtaps-1)/2:(RRCtaps-1)/2)*deltat;
% figure;
% plot(f,h_freq);

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
EbN0 = -1;
signal_power = (trapz(abs(signal_hrrc_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/(length(bits_tx_coded)); % energy per bit

N0 = Eb/10.^(EbN0/10);
NoisePower = 2*N0*fsampling;
noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx),1)+1i*randn(length(signal_hrrc_tx),1));

signal_rx = signal_hrrc_tx + noise;
signal_hhrc_rx = conv(signal_rx, h_time);
signal_hhrc_rx_trunc = signal_hhrc_rx(RRCtaps:end-RRCtaps+1);

%% Downsampling
signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

% %% Demapping
encoded_bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))';
% 
% %% Hard Decoding
% decoded_bits_rx = zeros(Nbits,1);
% for k=1:Npackets
%     packet_rx = encoded_bits_rx(1+(k-1)*codedWordLength : k*codedWordLength);
%     decoded_packet_rx = LdpcHardDecoder(packet_rx, H, tannerGraph, 10);
%     decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
% end
decoded_bits_rx = LdpcSoftDecoder(encoded_bits_rx, signal_rx_down, H, graph, N0, 10)';

%% Compare TX and RX
if bits_tx_coded == decoded_bits_rx
    disp('ok')
else 
    disp('RX not equal to TX')
end
