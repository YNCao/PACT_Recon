% GENERATING SYSTEM MATRIX (COEFFICIENT MATRIX) FOR PACT

% Author: Yunning Cao
% Version: 1.0
% Created date: 04/05/2019
% Last update: 07/28/2019 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Optional parameters
N = 64;                                            % size of image [pixel]
sensor_radius = floor(N * sqrt(2) / 2 + 2) - 1;     % sensor radius [pixel]
sensor_num = 128;
theta_start = 0;                                    % [deg]
range = 360;
theta_end = range-range/sensor_num;                     % [deg]

%%  
theta = linspace(theta_start, theta_end, sensor_num);       % angular distribution of sensors
I0 = zeros(N);
P = paradon(I0, theta, sensor_radius, 1);        % spherical radon transform
[m, n] = size(P);                              % size of the sensor data (P)

% coefficient matrix
A = sparse(m*n, N*N);

% separate system matrix into batches to solve memory overflow
batch = 8;
%%
tic
for i = 1 : batch
    % f = waitbar(0, 'Please wait...');
    A1=zeros(m*n, N*N/batch);    %coefficient matrix
    for ii=1:N*N/batch
        temp=zeros(N);
        temp(ii + (i-1)*N*N/batch)=1;
        A1(:,ii)=reshape(paradon(temp,theta,sensor_radius,1),m*n,1);
        if ~mod(ii, 100)
            % waitbar(ii / N^2, f, sprintf('%.2f%% finished',100 * ii / N^2));
            ii/N/N+(i-1)/batch
        end
    end
    % close(f);
    A(:, (1+(i-1)*N*N/batch : i*N*N/batch)) = A1;
    clear A1
    i
end
toc

%save coef mat
filename = ['CoefMat_', num2str(N), '_', num2str(theta_start), '_', num2str(theta_end), '_', num2str(sensor_num), '.mat'];
save(['.\coef_mat\', filename], 'A');
%clearvars -except A