edit%PACT_Recon
clear
clc
%%
% load data

data_type = 2;

% set image
N = 128;    %size of image
sensor_radius = floor(N * sqrt(2) / 2 + 2) - 1;
I0=phantom(N);
%I0=makeDisc(N,N,N/2,N/2,10);
% % I0=zeros(N);
% % I0(N/4:3*N/4,N/4:3*N/4)=1;

% add noise
I0 = I0 + 0.03 * randn(N, N);


% set sensor parameter
sensor_num = 128;
theta_start = 0;                                    % [deg]
theta_end = 360-360/sensor_num;                     % [deg]
theta = linspace(theta_start, theta_end, sensor_num);       % angular distribution of sensors

% velocity potential
P = paradon(I0, theta, sensor_radius, 1);
% diffrential
if data_type == 2
    P_ = P;
    P_(end-1:end,:) = [];
    P_ = [zeros(2,sensor_num);P_];
    P = (P-P_)/2;
end

exp = 1;
if exp == 1
    P = load('ustc.mat');
    P = P.pa;
%     P = fliplr(P);
end

[m,n] = size(P);
% P = P + 0.03 * randn(m, n);
% P(P<0) = 0;
p = reshape(P, m*n, 1);

%%
% load coefficient matrix

load('coef_mat\CoefMat_128_0_357.1875_128.mat');

if data_type == 2
    A_ = A;
    A_(end-1:end, :)=[];
    A_ = [zeros(2, size(A,2)); A_];
    A = (A-A_)/2;
end
%%
%Reconstruction

% choosing reconstruction method
% case 1: ML_EM, case 2: ART, case 3: SART

method = 2;
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
                imshow(I_);
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
                    imshow(I_,[]);
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
                imshow(I_);
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
