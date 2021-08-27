function [initSNR, SNRraddB, initErr, ErrorP] = cyclic2 (SNRin)

% Parameters
K = 128;                        % # of subcarriers
Rcom = 6;                       % channel length
Rrad = 64;                      % # of range cells
sigma = sqrt(1)/sqrt(K);        % std.

% Communication Channel
hdB = zeros(K, 1);
hdB(1:Rcom, 1) = [-6.0 0.0 -7.0 -22.0 -16.0 -20.0];    % power profile
h = zeros(K, 1);
h(1:Rcom, 1) = 10.^(hdB(1:Rcom, 1)/10);
H = fft(h);                                            % frequency response

% Threshold
SNRmin = 10.^(-30/10);                      % threshold SNR : -30dB per each subcarrier
rho = sigma * sqrt(SNRmin) ./ abs(H);       % threshold power of each subcarrier rho = [rho(0) rho(1) ... rho(K-1)]

% Initialization
d = randn(K, 1);
d = sqrt(SNRin) * d / norm(d);                  % OFDM freq domain d = [d(0) d(1) ... d(K-1)]

Rand = orth(randn(K, K));
Q = sqrt(K) * sqrt(SNRin) * Rand(:, 1:Rrad);    % semiunitary matrix
s = K * ifft(d);                                % OFDM time domain s = [s(0) s(1) ... s(K-1)]
S = zeros(K, Rrad);
for idx = 1 : Rrad
    S(:,idx) = circshift(fliplr(s), K-Rrad+idx);
end

% Objective Function
y = norm(S'*S - (Q'* Q), 'fro');

z = y;

% Initial SNR, Initial Error Probability
SS = S' * S;
SNRrad = sum(1 ./ diag(inv(SS))) / Rrad;
initSNR = 10 * log(SNRrad) / log(10);

SNRcom = abs(d).^2 .* abs(H).^2 / sigma^2;
initErr = sum(erfc(SNRcom ./ sqrt(2))) ./ K;

while true
    y = z;
    % #1 : Obtain S given Q
    S = alg1 (SNRin, K, Rrad, rho, Q);
    
    % #2 : Obtain Q given S
    [U, ~, V] = svd(S');
    Vtilde = V(:, 1:Rrad);
    Q = sqrt(K) * sqrt(SNRin) * Vtilde * U';
    
    z = norm(S'*S - (Q' * Q), 'fro');
    
    SS = S' * S;
    
    SNRrad = sum(1 ./ diag(inv(SS))) / Rrad;
    SNRraddB = 10 * log(SNRrad) / log(10);
    
    s = S(:,Rrad);
    d = (1/K) * fft(fliplr(s));
    
    SNRcom = abs(d).^2 .* abs(H).^2 / sigma^2;
    ErrorP = sum(erfc(SNRcom ./ sqrt(2))) ./ K;
    
    if abs(z-y) < 1e-5      % stopping criterion
        break
    end
    
end

end