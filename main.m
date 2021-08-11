function [] = main ()

% Parameters
K = 128;                    % # of subcarriers
Rcom = 6;                   % channel length
Rrad = 64;                  % # of range cells
sigma = 1;                  % std.          (input SNR = 1)

% Communication Channel
hdB = zeros(K, 1);
hdB(1:Rcom, 1) = [-6.0 0.0 -7.0 -22.0 -16.0 -20.0];    % power profile
h = zeros(K, 1);
h(1:Rcom, 1) = 10.^(hdB(1:Rcom, 1)/10);
H = fft(h);                                            % frequency response

hrad = h(1:Rrad, 1);
h_h = hrad * hrad';

% Threshold
SNRmin = 10.^(-30/10);                      % threshold SNR : -30dB
rho = sigma * sqrt(SNRmin) ./ abs(H);       % threshold power of each subcarrier rho = [rho(0) rho(1) ... rho(K-1)]

% Initialization
d = randn(K, 1);
d = d / norm(d);                % OFDM freq domain d = [d(0) d(1) ... d(K-1)]

Rand = orth(randn(K, K));
Q = Rand(:, 1:Rrad);            % semiunitary matrix
s = K * ifft(d);                % OFDM time domain s = [s(0) s(1) ... s(K-1)]
S = zeros(K, Rrad);
for idx = 1 : Rrad
    S(:,idx) = circshift(fliplr(s), K-Rrad+idx);
end

% Objective Function
y = norm(S'*S - Q'*Q, 'fro');

jdx = 0;
z = y;

iterations = zeros(0, 0);
ZZZ = zeros(0, 0);
SNRplot = zeros(0, 0);

while true
    y = z;
    % #1 : Obtain S given Q
    S = alg1 (K, Rrad, rho, Q);
    
    % #2 : Obtain Q given S
    [U, ~, V] = svd(S');
    Vtilde = V(:, 1:Rrad);
    Q = sqrt(K) * Vtilde * U';
    
    z = norm(S'*S - Q'*Q, 'fro');
    
    SS = S' * S;
    
    SNRrad = sum(diag(h_h) ./ diag(inv(SS)));
    SNRraddB = 10 * log(SNRrad) / log(10);
    
    jdx = jdx + 1;
    iterations = horzcat(iterations, jdx);
    ZZZ = horzcat(ZZZ, z);
    SNRplot = horzcat(SNRplot, SNRraddB);
    
    if abs(z-y) < 1e-5      % stopping criterion
        break
    end
    
end

figure
plot(iterations, ZZZ, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Objective Function');
title('Objective Function Convergence');
grid on

figure
plot(iterations, SNRplot, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Radar SNR (dB)');
title('Radar SNR Convergence');
grid on

figure
SS = S' * S;
[xx, yy] = meshgrid(1:Rrad, 1:Rrad);
plot3(xx, yy, abs(SS), 'LineWidth', 1.5);
zlabel('|{S^H}S|')
title('Diagonality of {S^H}S');
grid on

end