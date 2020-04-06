function [P] = paradon(I, theta, R, deltaR)
%Spherical Radon Transform for flat phantom in 3D(velocity potential)

%I          ->     Original image;  
%theta      ->     Any array representing detector distribution  0-360 [бу]
%R          ->     Radius of detector circle   [pixel]
%deltaR     ->     Time resolution             [pixel]


%size of original image
[m, n] = size(I);

%coordinate of the image center
Ox = (n + 1) /  2; 
Oy = (m + 1) / 2;  

%coordinate of detectors
Dx = Ox + R * cos(theta * pi / 180); 
Dy = Oy - R * sin(theta * pi /180);  

%size of projection matrix
sensor_num = length(theta);          %sensor number
timebin_num = ceil(2 * R / deltaR);  %time bin number in one sensor
P = zeros(timebin_num, sensor_num);


II = reshape(I,m * n,1);

%grid of original image
xy = (1: m * n)';
x = ceil(xy / m); 
y = mod(xy - 1, m)+ 1;

%pixels outside the detector are set zeros
rr = sqrt((x-Ox) .^ 2 + (y-Oy) .^2);
II(rr > R) = 0;
II(rr == R) = 0;

for j = 1 : sensor_num
    d2 = sqrt((x - Dx(j)) .* (x - Dx(j)) + (y - Dy(j)) .* (y - Dy(j)));  %distance from jth detector to every pixel
    d2_ = d2 / deltaR;
    integer = ceil(d2_);      %integer of distance
    remainder = d2_ - integer +1;
    
    % integer distance
    p = accumarray(integer, II .* (1 - abs(remainder - 0.5))); 
    
    % precise modification (distribute each pixel number into neighbored sample)
    sifter = remainder < 0.5;
    rup = remainder .* ~sifter + sifter * 0.5; 
   % rup(sifter) = 0.5;
    up = accumarray(integer, II .* (rup - 0.5));
    up(end) = up(end) + up(end - 1);
    up(2: end-1) = up(1: end-2); 
    up(1) = 0; 
    
    rdown = remainder .* sifter + ~sifter * 0.5; 
    %rdown(~sifter) = 0.5;
    down = accumarray(integer, II .* (0.5 - rdown));
    down(1) = down(1) + down(2);
    down(2: end - 1) = down(3: end); 
    down(end) = 0; 
    
    %calculate the integral
    P_ = p + up + down;
    lenP_ = length(P_);
    if lenP_ > timebin_num
        P(:, j) = P_(1: timebin_num);       
    else
        P(1: lenP_, j) = P_;      
    end
end

%R_ij = v_s * t_ij
r_ij = (deltaR * (1:timebin_num)' - deltaR / 2);
R_ij = repmat(r_ij, 1, sensor_num);
P = P ./ R_ij;
% end