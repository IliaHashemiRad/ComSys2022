function [Q_encoded, delta] = ADC(x, nu, fs, F, plot)
    x = x(1:end)';
    sam_x = x(1:fs/F:end);
    [Q, delta] = uniform_quan(nu, sam_x);
    Q_encoded = encoder(Q, delta, min(x));
    if(plot)
        ADC_plot(x,F,fs,sam_x,Q)
    end
end

function [Q , delta] = uniform_quan(v, x)
    N = 2^v;
    delta = ceil((max(x) - min(x)))/N;
    Q = zeros(1,length(x));
    for n=0:N-1
       Q(x>=min(x)+n*delta & x<=min(x)+(n+1)*delta) = min(x)+n*delta + delta/2;
    end
end 

function x_encoded = encoder(x_quan, delta, min_x)
    level = (x_quan-delta/2-min_x) / delta;
    x_encoded = de2bi(level, 'left-msb');
end

function ADC_plot(x,F,fs,sam_x,Q)
    dt = 1/F;
    t = 0:dt:length(x)/F-dt;
    t2 = 0:fs/F*dt:length(x)/F-dt;
    figure;
    plot(t,x); hold on;
    stairs(t2,sam_x); hold off;
    title('Input signal and the sampled Signal');
    xlabel('t');
    ylabel('Amplitude');
    legend('Original input', 'sampled input');
    figure;
    stem(t2,sam_x); hold on;
    stem(t2,Q);
    title('sampled data and its quantized values');
    xlabel('t');
    ylabel('Amplitude');
    legend('Analog', 'Quantized');
end