%PACT_Recon
clear
clc
%%
% set experiment type and signal type
% exp_type: 1 -> num_data;            0 -> exp_data;     
% sig_type: 1 -> velocity potential;  0 -> acoustic pressure

exp_type = 0;
sig_type = 1;

% data path
data_path = 'exp_data/19-06-10/3MHz/hair';
% data_path = 'exp_data/ustc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set sensor parameter
N = 128;                                                 % size of image [pixel]
Radius = floor(N * sqrt(2) / 2 + 2) - 1;                 % Radius [pixel]
sensor_radius = 25;                                      % sensor radius [mm]
sensor_num = 128;   
fs = 62.5;                                               % sampling rate [MHz]
% theta range
theta_start = 0;                                         % [deg]
range = 360;                                             % [deg]
theta_end = range-range/sensor_num;                      % [deg]
theta = linspace(theta_start, theta_end, sensor_num);    % angular distribution of sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal processing
% data_type
% 1: structured data including SpeedOfSound, RFData and Index
% 0: non-structured data
data_type = 1;    
fc = 8;                                                 % cut off frequency when filtering
delay = 0;                                               % time delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate differential matrix
order = 4;
tmp_img = ones(N);
tmp_sig = paradon(tmp_img, theta, Radius, 1);
[m, n] = size(tmp_sig);
clear tmp_img tmp_sig
D = diff_mat(m*n, 4);

%%
% load data
[pa_sig, neg_vp_sig] = signal_process(data_type, data_path, fs, fc, delay, sensor_radius, m, n);

if exp_type == 1
    % simulation
    I0=phantom(N);
    % add noise
    I0 = I0 + 0.03 * randn(N, N);
    % velocity potential
    P = paradon(I0, theta, Radius, 1);
    % velocity potential -> acoustic pressure
    if sig_type == 0
        D_data = diff_mat(m, 4);
        P = D_data * P;
    end
    P = P-min(min(P))+0.01;
elseif exp_type == 0
    % experiment
    if sig_type == 0
        P = pa_sig;
    elseif sig_type == 1
        P = neg_vp_sig;
    end
end

% P = P + 0.03 * randn(m, n);
% P(P<0) = 0;
p = reshape(P, m*n, 1);

%%
% load coefficient matrix
load('coef_mat\CoefMat_128_0_357.1875_128.mat');
if sig_type == 0
    % if using acoustic pressure signal
    A = D * A;
end
%%
%Reconstruction
% choosing reconstruction method
% case 1: ML_EM, case 2: ART, case 3: SART

method = 1;
switch method
    case 1
        % ML_EM
        I = ones(N*N,1);    %initialization
        t = zeros(N*N,100);
        SumA1 = sum(A);
        fmat = moviein(100);
        imshow(reshape(I, N, N))
        fmat(1) = getframe;
        for i = 1:100
            t(:,i) = I;
            p_ = A * I;
            delta = p ./ p_;
            delta(isnan(delta)) = 0;
            delta(isinf(delta)) = 0;
            I=I.*(sum(delta.*A)./SumA1)'; %normalization!!
            I(isinf(I)) = 0;
            % movie
            if mod(i, 1) == 0
                I_ = reshape(I,N,N);
                I_ = full(I_);
                imagesc(I_);
                fmat(i+1) = getframe;
            end
            % end iteration till convergence
            if norm(I-t(:,i))/norm(t(:,i)) < 5e-3
                break
            end
        end
        
    case 2
        % ART
        I = zeros(N*N, 1);    %initialization
        SumA2 = sum(A.^2, 2);
        t = zeros(N*N, 100);
        fmat = moviein(100);
        for ii = 1:30
            t(:, ii) = I;
            for i = 1:m*n
                p_ = A(i, :) * I;
                delta = p(i) - p_;
                C = delta .* A(i,:)' / SumA2(i);
                C(isnan(C)) = 0;
                I = I + C;
                if mod(i, 100) == 0
                    I_ = reshape(I,N,N);
                    imagesc(I_);
                    fmat(i) = getframe;
                end
            end
        end
        
    case 3  
        % SART
        I = ones(N*N, 1);    %initialization
        t = zeros(N*N, 100);
        SumA2 = sum(A, 2);
        fmat = moviein(100);
        for ii = 1:30
            t(:, ii) = I;
            for i = 1:n
                p_ = A((i-1)*m+1:i*m, :) * I;
                delta = p((i-1)*m+1:i*m) - p_;
                Dij = A((i-1)*m+1:i*m, :) .* delta ./ SumA2((i-1)*m+1:i*m);
                Dij(isnan(Dij)) = 0;
                D = sum(Dij);
                C = (D./sum(A((i-1)*m+1:i*m, :)))';
                C(isnan(C)) = 0;
                C(isinf(C)) = 0;
                I = I + C;
                if mod(i,2) == 0
                    I_ = reshape(I,N,N);
                    imagesc(I_);
                    fmat(i) = getframe;
                end
            end
            if norm(I-t(:,ii))/norm(t(:,ii)) < 1e-2
                break
            end
        end
        
    case 4
        I = ones(N*N, 1);
        t = zeros(N*N, 100);
        alpha = 0.01;
        SumA1 = sum(A)*10^17;
        fmat = moviein(100);
        for i = 1:100
            t(:,i) = I;
            res = sum(A*I-p);
            I = I - alpha*res*SumA1';
            if mod(i,2) == 0
                I_ = reshape(I,N,N);
                imagesc(I_);
                fmat(i) = getframe;
            end
        end
end

%%
%========== O U T P U T ========
II=reshape(I,N,N);
II = full(II);
figure
imshow(II)
saveas(gcf, ['method_', num2str(method), '_', num2str(N)], 'tif')
%========== A N A L Y Z E =======
exp = repmat(p,1,size(t,2));
est = A*t;
r = sum((est-exp).^2);
figure
plot(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save resuts
filename = ['results_', num2str(N), '_', num2str(theta_start), '_', num2str(theta_end), '_', num2str(sensor_num), '.mat'];
save(['.\result\', filename], 't');
