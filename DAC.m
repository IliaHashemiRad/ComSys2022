function filtered_up = DAC(data,delta,min_x,F,fs)
    A_data = bi2de(data,'left-msb')*delta + delta/2 + min_x;
%     A_data = A_data';
%     diff = -[0 A_data] + [A_data 0];
%     diff = diff(2:end-1);
%     diff = upsample([diff 0], 4);
%     A_data = upsample(A_data , 4);
%     A_data1 = [A_data 0 0 0];
%     A_data2 = [0, A_data + diff/4, 0, 0];
%     A_data3 = [0 0, A_data + diff/2, 0];
%     A_data4 = [0 0 0, A_data + 3*diff/4];
%     res = A_data1 + A_data2 + A_data3 + A_data4;
%     Upsampled_data = res(1:end-6);
%     Upsampled_data = Upsampled_data';

    upsampled = upsample(A_data,fs/F);
    filtered_up = lowpass(upsampled,1000,F);
end