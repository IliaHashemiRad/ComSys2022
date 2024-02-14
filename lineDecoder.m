function bit = lineDecoder(input, t2, baudRate)
    d = round((1/baudRate / (t2(2) - t2(1))));
    index = 1 : d : length(t2);
    td = 0;
    ts = 0;
    index = index + td + ts;
    syms = input(index);
    bit = syms>0;
end