clear; clc; close all;

ar = audiorecorder;

duration = 10;
% record(ar,duration);
% pause(duration);
% x = getaudiodata(ar);
% F = ar.SampleRate;
[x, F] = audioread('sound2.wav');
% F = 100;
% t = 0:10e-6:0.03-10e-6;
% x = cos(2*pi*F.*t);
process_data(x,F);
pause(duration);
disp('Finish...');

%% Functions

function process_data(x,F)
    nu = 8;
    x = x(1:end)';
    dt = 1/F;
    t = 0:dt:length(x)/F-dt;

    fs = 4*F;
    % t2 = 0:1/fs:length(x)/fs-1/fs;
    t2 = 0:fs/F*dt:length(x)/F-dt;
    sam_x = x(1:fs/F:end);
    disp(length(sam_x));
    plot(t,x); hold on;
    stairs(t2,sam_x); hold off;

    [Q, delta] = uniform_quan(nu, sam_x);
    figure;
    stem(t2,sam_x); hold on;
    stem(t2,Q);

    Q_encoded = encoder(Q, nu, delta, min(x));

    % Line Coding
    beta = 40;
    baudRate = 1000;
    [x_lineCoded, t2] = lineCoder(Q_encoded', beta, baudRate);
    figure;
    plot(t2,x_lineCoded);

    % Channel
    B = baudRate/2;
    sigma = 0.5;
    chnl_out = channel(x_lineCoded, B, sigma, fs);
    figure;
    plot(t2,chnl_out);

    % Decoder
    detected_bits = lineDecoder(chnl_out,t2,baudRate);
    data = reshape(detected_bits,[nu length(detected_bits)/nu])';

    disp(length(Q_encoded));
    disp(sum(data == Q_encoded, 'all')/length(Q_encoded)/nu);
    
    % DAC
    A_data = bi2de(data,'left-msb')*delta + delta/2 + min(x);
    % audiowrite('reconstructedSound.wav',A_data,5000);
    sound(A_data,F^2/fs);
end

function [Q , delta] = uniform_quan(v, x)
    N = 2^v;
    delta = ceil((max(x) - min(x)))/N;
    Q = zeros(1,length(x));
    for n=0:N-1
       Q(x>=min(x)+n*delta & x<=min(x)+(n+1)*delta) = min(x)+n*delta + delta/2;
    end
%     Q(x>-1/2*delta & x<1/2*delta) = 0;
%     for n=1:N/2-1
%         Q(x>(2*n-1)/2*delta & x<(2*n+1)/2*delta) = n*delta;
%         Q(x<-(2*n-1)/2*delta & x>-(2*n+1)/2*delta) = -n*delta;
%     end
%     Q(x>(N-1)/2*delta) = (N/2-1)*delta;
%     Q(x<-(N-1)/2*delta) = -(N/2-1)*delta;
    
end 

function x_encoded = encoder(x_quan, nu, delta, min_x)
    level = (x_quan-delta/2-min_x) / delta;
    disp(min(level));
%     level = level + (2^nu-2)/2;
    x_encoded = de2bi(level, 'left-msb');
end

function [result, t2] = lineCoder(bit_stream, beta, baudRate)
    A = 10;
    [Rows,Cols] = size(bit_stream);
    res = 10;   %resolution
    t2 = 0:1/baudRate/res:(Rows*Cols-1/res)/baudRate;
    result = zeros(1,length(t2)+400);
    t3 = -200*1/baudRate/res:1/baudRate/res:199*1/baudRate/res;
    p = cos(2*pi*beta.*t3).*sinc(baudRate*t3) ./ (1-(4*beta.*t3).^2) ;
    figure;
    plot(t3,p);
    for i=0:Rows*Cols-1
        if(bit_stream(i+1) == 1)
            result(1+i*res:400+i*res) = result(1+i*res:400+i*res) + A*p;
        end
        if(bit_stream(i+1) == 0)
            result(1+i*res:400+i*res) = result(1+i*res:400+i*res) - A*p;
        end
    end
    result = result(200:end-201);
end

function out = channel(in, B, sigma, Fs)
    filtered_x = lowpass(in,B,Fs);
    noise = sigma * randn(1,length(in));
    figure;
    plot(noise);
    out = filtered_x + noise;
end

function bit = lineDecoder(input, t2, baudRate)
    d = round((1/baudRate / (t2(2) - t2(1))));
    index = 1 : d : length(t2);
    syms = input(index);
    bit = syms>0;
end