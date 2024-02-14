function out = channel(in, B, sigma2, Fs, t2, plot)
    filtered_x = lowpass(in,B,Fs);
    noise = sqrt(sigma2) * randn(1,length(in));
    out = filtered_x + noise;
    if(plot)
        channel_plot(t2,out,noise)
    end
end

function channel_plot(t2,chnl_out,noise)
   figure;
   plot(t2,noise);
   title('noise added to channel input');
   figure;
   plot(t2,chnl_out);
   title('Channel noisy output');
   
end