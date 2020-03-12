% Using IIR filter to realize EQ. (Not real-time realization.)
% Refence: https://github.com/Mosab-Mohamed/Audio-Equalizer-In-MATLAB.git

clear
[data, fs] = audioread('sample.wav');
data = data(:, 1);

filter_n = 9;
% f_center = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] * band + band / 2;
% Q = f_center ./ band;
f = [65, 125, 250, 500, 1000, 2000, 4000, 6000];
gain = 10 .^ ([-10, -20, -20, 0, -6, -6, 20, 20, 20] ./ 20);
order = 2;

Hz = cell(1, filter_n);

% lowpass
Fc2 = f(1);
h = fdesign.lowpass('N,F3dB', order, Fc2, fs);
Hz{1} = design(h, 'butter');

% bandpass
for i = 2:filter_n-1
    Fc1 = Fc2;
    Fc2 = f(i);
    h   = fdesign.bandpass('N,F3dB1,F3dB2', order, Fc1, Fc2, fs);
    Hz{i} = design(h, 'butter');
end

% highpass
Fc = Fc2;
h  = fdesign.highpass('N,F3dB', order, Fc, fs);
Hz{filter_n} = design(h, 'butter');

y = zeros(length(data), filter_n);
for i = 1:filter_n
    y(:,i) = filter(Hz{i}, data) * gain(i);
end

out = sum(y, 2);

audiowrite('eq_out.wav', out, fs)