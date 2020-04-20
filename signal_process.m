% load photoacoustic signal
clear; clc;
% data_path = './exp_data/19-06-10 hair tube and vessel';
% load([data_path, '/3MHz/hair/RFData_1.mat']);

data_path = './exp_data/ustc';
data = load([data_path, '/rawdata.mat']);

SpeedOfSound = 1.530;                     % [mm/us]
fs = 125;                                 % [MHz]
fc = 4;                                   % [MHz] % cutting frequency


delta_t_ = SpeedOfSound / (fs * 1000);    % [mm]

% preprocess
% pa_signal = double(RFData(IndexInfo.PAStartSample : IndexInfo.PAEndSample, :));
% t_s = (IndexInfo.PAStartSample : IndexInfo.PAEndSample)' * delta_t_;   % [mm]

pa_signal = data.PARF;
[m, n] = size(pa_signal);
t_s = (1:size(pa_signal,1))' * delta_t_;


% filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=linspace(-fs/2, fs/2-fs/m, m);
Window = zeros(m, 1);
AlgorWindow = 'Hanning';
for i = 1:m

    if abs(f(i)) < (fc)

        switch AlgorWindow
            case 'Hanning'
                Window(i) = 0.5 + 0.5 * cos(pi * (f(i)) / fc);
            case 'Blackman'
                Window(i) = 0.42 + 0.5 * cos(pi * f(i) / fc) + 0.08 * cos(2 * pi * f(i) / fc); % blackman window
        end

    end

end

Wc = Window * ones(1, n);

f_pa = fftshift(fft(pa_signal), 1);

W_f = f_pa .* Wc;

W_f_s = ifftshift(W_f, 1);

f_pa_signal = ifft(W_f_s);

% delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = 1; % time of flight in us,[us]
truncation = abs(round(delay * fs)); % the number of sampling truncated

if delay <= 0
    % zero padding, adding truncation number of zeros before rawdata
    f_pa_signal = [zeros(truncation, n); f_pa_signal];
    f_pa_signal(m + 1:end, :) = [];
else
    % remove redundant signals at the beginning
    f_pa_signal(1:truncation, :) = [];
    f_pa_signal = [f_pa_signal; zeros(truncation, n)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_pa_signal(1:1500,:)=0;
% f_pa_signal(3000:end,:)=0;
pa_d = real(f_pa_signal(1:23:4180, 1:2:end));
bias = mean(mean(pa_d(1:50, :)));
%bias = min(min(pa_d));
pa_b = pa_d-bias;
pa = real(pa_b);
%pa(1:5,:)=0;
save('ustc.mat','pa')

v_p=cumtrapz(pa);
P = -v_p-min(min(v_p));
