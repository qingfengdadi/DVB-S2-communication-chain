%% Projet modulation & coding
addpath(genpath('Code encodeur'));
addpath(genpath('Code mapping-demapping'));
addpath(genpath('Code HRC'));
addpath(genpath('Data'));

%% Parameters
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 10; % oversampling factor
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols
l = 3;
Nbps = 4; % number of bits per symbol
modulation = 'qam'; % type of modulation 
beta = 0.3; % roll-off factor

Npackets = 5000; % choose such that Nbits/Nbps is an integer
packetLength = 128;
codedWordLength = 2*packetLength;
Nbits = Npackets*packetLength; % bit stream length
bits_tx = randi(2,Nbits,1)-1;

%% Mapping of uncoded signal
signal_uncoded = mapping(bits_tx,Nbps,modulation);

%% Upsampling
signal_tx_uncoded = upsample(signal_uncoded,M);

%% Implementation of HRC
RRCtaps = 365;
stepoffset = (1/RRCtaps)*fsampling;
highestfreq = (RRCtaps-1)*stepoffset/2;
f = linspace(-highestfreq,highestfreq,RRCtaps);
Hrc = HRC(f,Tsymb,beta);
h_t = ifft(Hrc);
h_freq = sqrt(fft(h_t/max(h_t)));
h_time = fftshift(ifft(ifftshift(h_freq)));

%% Convolution
signal_hrrc_tx_uncoded = conv(signal_tx_uncoded, h_time);

%% Noise through the channel uncoded
EbN0 = -5:25;
BER = zeros(length(EbN0),5);
signal_power_uncoded = (trapz(abs(signal_hrrc_tx_uncoded).^2))*(1/fsampling); % total power
Eb = signal_power_uncoded*0.5/(Npackets*packetLength); % energy per bit

for j = 1:length(EbN0)
    N0 = Eb/10.^(EbN0(j)/10);
    NoisePower = 2*N0*fsampling;
    noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx_uncoded),1)+1i*randn(length(signal_hrrc_tx_uncoded),1));

    signal_rx = signal_hrrc_tx_uncoded + noise;
    signal_hhrc_rx = conv(signal_rx, h_time);
    signal_hhrc_rx_trunc = signal_hhrc_rx(RRCtaps:end-RRCtaps+1);

    %% Downsampling
    signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

    %% Demapping
    uncoded_bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
    BER(j,l) = length(find(bits_tx ~= uncoded_bits_rx'))/length(uncoded_bits_rx');
end

%% Plot BER results
semilogy(EbN0,BER(:,1),'-',EbN0,BER(:,2),'-',EbN0,BER(:,3),'-',EbN0,BER(:,4),'-');ylim([1e-4;1E0]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('4QAM','16QAM','64QAM','256QAM');
title('BER curves')
grid on;
