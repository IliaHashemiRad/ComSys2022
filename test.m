clc; clear; close all;

% F = 1000;
% t = 0:10e-6:0.003-10e-6;
% x = 2*cos(2*pi*F.*t)+1;

[x, F] = audioread('sound2.wav');

nu = 8;
x = x';
  
fs = 4*F;
[Q_encoded, delta] = ADC(x, nu, fs, F, 1);

beta = 30;
baudRate = 1000;
[x_lineCoded, t2] = lineCoder(Q_encoded', beta, baudRate,10, 1);
    
% Channel
B = 600;
N0 = 2*10e-4;
sigma2 = B*N0;
chnl_out = channel(x_lineCoded, B, sigma2, fs, t2, 1);
    
% Decoder
detected_bits = lineDecoder(chnl_out,t2,baudRate);
data = reshape(detected_bits,[nu length(detected_bits)/nu])';
disp('Matching percentage of extracted data and encoded data :');
disp(sum(data==Q_encoded,'all')/length(Q_encoded)/nu*100);

% DAC
A_data = DAC(data,delta,min(x));
% figure;
% plot(A_data); hold on;
% plot(x(1:fs/F:end));
% title('Reconstructed signal and the original one');
% legend('Reconstructed','Original');
sound(A_data,F^2/fs);