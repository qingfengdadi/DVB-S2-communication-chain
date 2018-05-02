%% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code decodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
clear; close all;

%% Parameters
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 10; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols
Nbps = 2; % number of bits per symbol
modulation = 'qam'; % type of modulation 
beta = 0.3; % roll-off factor

Npackets = 1; % choose such that Nbits/Nbps is an integer
packetLength = 128;
codedWordLength = 2*packetLength;
Nbits = Npackets*packetLength; % bit stream length
NcodedBits = Npackets*codedWordLength; % full coded word length
bits_tx = randi(2,Nbits,1)-1;
bits_tx_coded = zeros(NcodedBits,1);

%% LDPC encoder
H0 = makeLdpc(packetLength, codedWordLength, 0, 1, 3);
for k=1:Npackets
    packet_tx = bits_tx(1+(k-1)*packetLength : k*packetLength);
    [codedbits, H] = makeParityChk(packet_tx , H0, 0);
    bits_tx_coded(1+(k-1)*codedWordLength : k*codedWordLength) = [codedbits packet_tx];
end
tannerGraph = buildTannerGraph(H);

%% Mapping of encoded signal
% signal = mapping(bits_tx,Nbps,modulation);
signal_coded = mapping(bits_tx_coded,Nbps,modulation);

%% Upsampling
signal_tx = upsample(signal_coded,M);

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
signal_hrrc_tx = conv(signal_tx, h_time);

%% Noise through the channel
EbN0 = 5;
signal_power = (trapz(abs(signal_hrrc_tx).^2))*(1/fsampling); % total power
Eb = signal_power*0.5/(Npackets*codedWordLength); % energy per bit

N0 = Eb/10.^(EbN0/10);
NoisePower = 2*N0*fsampling;
noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx),1)+1i*randn(length(signal_hrrc_tx),1));

signal_rx = signal_hrrc_tx + noise;
signal_hhrc_rx = conv(signal_rx, h_time);
signal_hhrc_rx_trunc = signal_hhrc_rx(RRCtaps:end-RRCtaps+1);

%% Downsampling
signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

%% Demapping
encoded_bits_rx = (demapping(signal_rx_down,Nbps,modulation))';

%% Hard Decoding
decoded_bits_rx = zeros(Nbits,1);
for k=1:Npackets
    packet_rx = encoded_bits_rx(1+(k-1)*codedWordLength : k*codedWordLength);
    decoded_packet_rx = LdpcHardDecoder(packet_rx, H, tannerGraph, 2);
    decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
end

%% Compare TX and RX
if bits_tx == decoded_bits_rx
    disp('ok')
else 
    disp('RX not equal to TX')
end
