clear; clc; close all;

% ar = audiorecorder;
% duration = 5;
% record(ar,duration);
% pause(duration);
% x = getaudiodata(ar);
% F = ar.SampleRate;

[x, F] = audioread('sound2.wav');

process_data(x,F);
disp('Finish...');

%% Realtime Process
clear; clc; close all;

F = 60000;
fs = 4*F;
micReader = audioDeviceReader(F, F/2,"Driver","DirectSound");
speakerWriter = audioDeviceWriter(F, "Driver","DirectSound");
 
duration = 10;
tic;
while toc < 30
    audio = micReader();
    audio = process_data(audio, F);
    speakerWriter(audio);
end

%%
clear; clc; close all;

[x, F] = audioread('sound2.wav');
B = [100 200 300 400 450];
N0 = 10e-3;
o_SNR = zeros(1,5);
for i=1:5
    out = process_data(x,F,B(i),N0);
    o_SNR(i) = snr(out);
end
disp('    B         SNR(dB)');
disp([B; o_SNR]');

%%
clear; clc; close all;

[x, F] = audioread('sound2.wav');
N0 = [5*10e-5 10e-4 5*10e-4 10e-3 5*10e-3];
B = 100;
o_SNR = zeros(1,5);
for i=1:5
    out = process_data(x,F,B,N0(i));
    o_SNR(i) = snr(out);
end
disp('    N0         SNR(dB)');
disp([N0; o_SNR]');

%%
clear; clc; close all;

[x, F] = audioread('sound2.wav');
sigma2 = 0.5:0.1:4;
o_BER = zeros(1,length(sigma2));
for i=1:length(sigma2)
    o_BER(i) = process_data(x,F,sigma2(i));
end
plot(sigma2,o_BER);
xlabel('sigma2');
ylabel('BER(%)');
title('Bit Error Rate for A=5');

%%
clear; clc; close all;

[x, F] = audioread('sound2.wav');
A = 2:0.5:10;
o_BER = zeros(1,length(A));
for i=1:length(A)
    o_BER(i) = process_data(x,F,A(i));
end
plot(A,o_BER);
xlabel('A');
ylabel('BER(%)');
title('Bit Error Rate for N0=5*10e-4');

%%
clear; clc; close all;

[x, F] = audioread('sound2.wav');
beta = 1:20:500;
o_BER = zeros(1,length(beta));
for i=1:length(beta)
    o_BER(i) = process_data(x,F,beta(i));
end
plot(beta,o_BER);
xlabel('beta');
ylabel('BER(%)');
title('Bit Error Rate for N0=5*10e-4');

%% 

clear; clc; close all;

[x, F] = audioread('sound2.wav');
nu = 2:2:10;
o_SNR = zeros(1,length(nu));
for i=1:length(nu)
    out = process_data(x,F,nu(i));
    o_SNR(i) = snr(out);
end
disp('    nu         SNR(dB)');
disp([nu; o_SNR]');

%% Functions

function A_data = process_data(x,F)
    nu = 8;
    x = x(1:end)';
    
    fs = 4*F;
    [Q_encoded, delta] = ADC(x, nu, fs, F, 0);

    % Line Coding
    beta = 40;
    baudRate = 1000;
    A=10;
    [x_lineCoded, t2] = lineCoder(Q_encoded', beta, baudRate, A, 0);
    
    % Channel
    B = baudRate/2;
%     N0 = 1.380649*10e-23 * (100+273);
    N0 = 10e-4;
    sigma2 = B*N0;
    chnl_out = channel(x_lineCoded, B, sigma2, fs, t2, 0);
    
    % Decoder
    detected_bits = lineDecoder(chnl_out,t2,baudRate);
    data = reshape(detected_bits,[nu length(detected_bits)/nu])';

    % DAC
    A_data = DAC(data,delta,min(x),F,fs);
    % audiowrite('reconstructedSound.wav',A_data,5000);
%     sound(A_data,F^2/fs);
    sound(A_data,F);
%     BER = 100 - (sum(data==Q_encoded,'all'))/(length(Q_encoded)*nu)*100;
end
