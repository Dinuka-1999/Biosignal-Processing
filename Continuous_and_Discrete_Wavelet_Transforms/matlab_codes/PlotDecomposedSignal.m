function PlotDecomposedSignal(signal, wavelet, n)

[C, L] = wavedec(signal, n, wavelet);
A = appcoef(C, L, wavelet);

figure('Name','Approximated Coefficients')
stem(A, 'filled');
title(['Approximated Coefficients when levels = ' num2str(n)]);

figure('Name','Decomposition of the Signal')

for level= 1:n
    d = detcoef(C, L, level);
    subplot(n, 1, level);
    stem(d,'filled');
    title(['Level:' num2str(level) ' Wavelet Function: ' wavelet]);
end
