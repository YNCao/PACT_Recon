function [pa_sig, neg_vp_sig] = signal_process(data_type, data_path, fs, fc, delay, sensor_radius, downsample_m, downsample_n)
% data_type
% 1: structured data including SpeedOfSound, RFData and Index
% 0: non-structured data
%
% fs = 62.5;                                % [MHz] % sampling rate
% fc = 10;                                  % [MHz] % cutting frequency     
% delay = 1;                                % [us]
% sensor_radius = 25                        % [mm]
% downsample_n = 182


if data_type == 1
%     data_path = './exp_data/19-06-10/3MHz/vessel';
    tag = strrep(data_path(3:end), '/', '_');
    data = load([data_path, '/RFData_1.mat']);
    SpeedOfSound = data.SpeedOfSound/1000;   % [mm/us]
    raw_sig = data.RFData(data.IndexInfo.PAStartSample:data.IndexInfo.PAEndSample,:);
elseif data_type ==0
%     data_path = './exp_data/ustc';
    tag = strrep(data_path(3:end), '/', '_');
    data = load([data_path, '/rawdata.mat']);
    SpeedOfSound = 1.530;                     % [mm/us]
    raw_sig = data.PARF;
end
  
delta_t_ = SpeedOfSound /fs;                % [mm]
[m, n] = size(raw_sig);
% t_s = (0:size(raw_sig,1)-1)' * delta_t_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove the starting pulse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filter
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

f_pa = fftshift(fft(raw_sig), 1);

W_f = f_pa .* Wc;

W_f_s = ifftshift(W_f, 1);

f_raw_sig = ifft(W_f_s);

% delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delay
d_sig = sig_delay(f_raw_sig, fs, delay);

%% downsample
m_index = round((1/2: downsample_m-1/2)/downsample_m * 2*sensor_radius/delta_t_);
n_index = round((0: downsample_n-1)/downsample_n * n);
if m_index(1)<1
    m_index=m_index+1;
end
if n_index(1)<1
    n_index=n_index+1;
end
s = d_sig(m_index, n_index);

bias = mean(mean(s(1:50, :)));
pa_sig = real(s-bias);

save(['./recon_data/',tag], 'pa_sig')

neg_vp_sig=-cumtrapz(pa_sig);
min_vp = min(min(neg_vp_sig));
if min_vp<0
    neg_vp_sig = neg_vp_sig-min_vp;
end
end


% customized function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = sig_delay(sig, fs, delay)
[m, n] = size(sig);
truncation = abs(round(delay * fs)); % the number of sampling truncated
if delay <= 0
    % zero padding, adding truncation number of zeros before rawdata
    sig = [zeros(truncation, n); sig];
    sig(m + 1:end, :) = [];
else
    % remove redundant signals at the beginning
    sig(1:truncation, :) = [];
    sig = [sig; zeros(truncation, n)];
end
end