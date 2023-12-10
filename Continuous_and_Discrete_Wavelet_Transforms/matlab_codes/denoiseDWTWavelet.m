function root_mean_sqr_err = denoiseDWTWavelet(C,L, x, wavelet, threshold, name)

    decomposed_filtering = C;
    
    %Setting the decompositions less than the threshold to 0
    decomposed_filtering(abs(decomposed_filtering) < threshold) = 0;
    
    denoised_signal = waverec(decomposed_filtering, L, wavelet);
    
    % Plotting the results
    figure ('Name',[name ' : ' wavelet]);
    plot(denoised_signal, 'LineWidth', 1);
    title(['Denoised Signal'  ' : ' wavelet ' (To be compared with: ' name ')'])
    xlabel('Sample Number (n)')
    ylabel('Amplitude');
    
    % Checking the disparity
    difference = x - denoised_signal;   
    root_mean_sqr_err = sqrt(sum(abs(difference).^2)/length(difference));

    figure('Name',['Comparison between the Original and Denoised ' name])
    plot(1:1:length(x), x, 'LineWidth', 1)
    hold on;
    plot(1:1:length(denoised_signal), denoised_signal, 'LineWidth', 1)
    hold off;
    xlim([0,length(x)]);
    title(['Comparison between the Original ' name ' and Denoised Signal '  ': ' wavelet ])
    xlabel('Sample Number (n)')
    ylabel('Amplitude');
    legend(['Original ' name], 'Denoised Signal')
end