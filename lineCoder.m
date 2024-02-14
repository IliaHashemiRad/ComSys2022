function [result, t2] = lineCoder(bit_stream, beta, baudRate, A, plot)
    [Rows,Cols] = size(bit_stream);
    res = 10;   %resolution
    t2 = 0:1/baudRate/res:(Rows*Cols-1/res)/baudRate;
    result = zeros(1,length(t2)+400);
    t3 = -200*1/baudRate/res:1/baudRate/res:199*1/baudRate/res;
    p = cos(2*pi*beta.*t3).*sinc(baudRate*t3) ./ (1-(4*beta.*t3).^2) ;
    for i=0:Rows*Cols-1
        if(bit_stream(i+1) == 1)
            result(1+i*res:400+i*res) = result(1+i*res:400+i*res) + A*p;
        end
        if(bit_stream(i+1) == 0)
            result(1+i*res:400+i*res) = result(1+i*res:400+i*res) - A*p;
        end
    end
    result = result(200:end-201);
    
    if(plot)
        lineCoder_plot(t2, result, baudRate, beta);
    end
end

function lineCoder_plot(t2, x_lineCoded, baudRate, beta)
    figure;
    plot(t2,x_lineCoded);
    title('Line coded output');
    xlabel('t');
    ylabel('Amplitude');
    res = 10;
    t3 = -200*1/baudRate/res:1/baudRate/res:199*1/baudRate/res;
    p = cos(2*pi*beta.*t3).*sinc(baudRate*t3) ./ (1-(4*beta.*t3).^2) ;
    figure;
    plot(t3,p)
    title('Nyquist pulse shaping signal');
    xlabel('t');
    ylabel('Amplitude');
end