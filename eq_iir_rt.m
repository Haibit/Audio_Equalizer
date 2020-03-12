% Using IIR filter to realize EQ. (Real-time realization for every sample point.)
% Direct form II

clear
[data, fs] = audioread('sample.wav');
data       = data(:, 1);

filter_n = 9;
f        = [65, 125, 250, 500, 1000, 2000, 4000, 6000];
gain     = 10 .^ ([-10, -20, -20, 0, -6, -6, 20, 20, 20] ./ 20);
order    = 2;
H        = design_iir(filter_n, order, fs, f);

w        = zeros(filter_n, order+1);
a_tmp    = zeros(1, order+1);
len_data = length(data);
out      = zeros(1, len_data);
filter_x = zeros(1, filter_n);
for i = 1:len_data
    w(:, 1:end-1) = w(:, 2:end);
    for j = 1:filter_n
        a_tmp(1) = H{j}(2,1);
        a_tmp(2:end) = -H{j}(2,2:end);
        w(j, end) = [data(i), w(j, end-1:-1:1)] * a_tmp';
        filter_x(j) = gain(j) * (H{j}(1,:) * w(j, end:-1:1)');
    end
    out(i) = sum(filter_x);
end

audiowrite('eq_rt_out.wav', out, fs);
subplot(2,1,1)
plot(data)
subplot(2,1,2)
plot(out)


%% functions to design IIR filter and get 'a' and 'b'
function [H] = design_iir(filter_n, order, fs, f)
Hz = cell(1, filter_n);
H  = cell(1, filter_n);

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

for i = 1:filter_n
    H{i} = zeros(2, order+1);

%     H{i}(1,:) = Hz{i}.sosMatrix(1:order+1) * Hz{i}.ScaleValues(1);       % b
%     H{i}(2,:) = Hz{i}.sosMatrix(order+2:end) * Hz{i}.ScaleValues(2);     % a

%     [H{i}(1,:), H{i}(2,:)] = sos2tf(Hz{i}.sosMatrix);                    % another way
%     H{i}(1,:) = H{i}(1,:) * Hz{i}.ScaleValues(1);
%     H{i}(2,:) = H{i}(2,:) * Hz{i}.ScaleValues(2);

    [H{i}(1,:), H{i}(2,:)] = tf(Hz{i});                                    % third way
end

end