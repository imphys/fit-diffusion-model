function g = GetDiffusionGradientDirection(CSAheader)
    pos = findstr(CSAheader, 'DiffusionGradientDirection');
    CSAheader = CSAheader(pos(1):end);
    pos = findstr(CSAheader, 'M');
    CSAheader = CSAheader(pos(1)+8:end);
    pos = findstr(CSAheader, 'M');
    CSAheader = CSAheader(pos(1)+8:end);
    g(1) = str2double(CSAheader(1:1+10));
    g(2) = str2double(CSAheader(1+28:1+28+10));
    g(3) = str2double(CSAheader(1+28+28:1+28+28+10));
end