function [] = convergence ()

iterations = cell(100, 1);
Zplot = cell(100, 1);
SNRplot = cell(100, 1);
ERRORplot = cell(100, 1);
SS = cell(100, 1);
iterations_len = zeros(100, 1);

for idx = 1 : 100
    [iterations{idx}, Zplot{idx}, SNRplot{idx}, ERRORplot{idx}, SS{idx}] = cyclic ();
    iterations_len(idx, 1) = length(iterations{idx});
end

min_len = min(iterations_len);
it = zeros(1, min_len);

for jdx = 1 : min_len
    it(jdx) = jdx;
end

Zpl = zeros(100, min_len);
SNRpl = zeros(100, min_len);
ERRORpl = zeros(100, min_len);

for kdx = 1 : 100
    ZZpl = Zplot{kdx, 1};
    SSNRpl = SNRplot{kdx, 1};
    EERRORpl = ERRORplot{kdx, 1};
    
    Zpl(kdx, :) = ZZpl(1, 1:min_len);
    SNRpl(kdx, :) = SSNRpl(1, 1:min_len);
    ERRORpl(kdx, :) = EERRORpl(1, 1:min_len);

end

Zmean = mean(Zpl);
SNRmean = mean(SNRpl);
ERRORmean = mean(ERRORpl);

figure
plot(it, Zmean, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Objective Function');
title('Objective Function Convergence');
grid on

figure
plot(it, SNRmean, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Radar SNR (dB)');
title('Radar SNR Convergence');
grid on

figure
semilogy(it, ERRORmean, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Error Probability');
title('Error Probability Convergence');
grid on

figure
[xx, yy] = meshgrid(1:64, 1:64);
plot3(xx, yy, abs(SS{1}), 'LineWidth', 1.5);
zlabel('|{S^H}S|')
title('Diagonality of {S^H}S');
grid on

end