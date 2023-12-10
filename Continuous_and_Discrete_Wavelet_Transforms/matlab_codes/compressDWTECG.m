function [comp_ratio, paramm_1] = compressDWTECG( decomposed, L, orig_signal, wavelet, name_signal)

    dord_coeff = sort(abs(decomposed(:)),'descend');
    
    % plotting the decomposition
    figure('Name',['Wavelet Coefficients of ' name_signal ' : '  'Wavelet - ' wavelet])
    stem(dord_coeff);
    title(['Wavelet Coefficients of ' name_signal ' : '  'Wavelet ' wavelet]);

    energy = 0;
    paramm_1 = 0;
   
    total_energy = sum(abs(dord_coeff).^2);
    for i = 1:length(dord_coeff)
        energy = energy + dord_coeff(i)^2;
        % Finding the coeffient whose cumulative sum is 99%
        if (round(energy/total_energy, 2) == 0.99)
            paramm_1 = i;
            break;
        end
    end
    thresh_arb = dord_coeff(paramm_1);

    decomposed_filtering = decomposed;
    decomposed_filtering(abs(decomposed_filtering) < thresh_arb) = 0;
    
    compressed_signal = waverec(decomposed_filtering, L, wavelet);
    
    % Plotting the results
    figure ('Name',[name_signal ' : ' wavelet]);
    plot(compressed_signal,'LineWidth',1);
    title(['Compressed Signal'  ' : ' wavelet ' (To be compared with: ' name_signal ')'])
    xlabel('Samples (n)')
    ylabel('Amplitude');
    
    % Checking the disparity
    comp_ratio=paramm_1/length(orig_signal);
    
    figure('Name',['Comparison between the Original and Compressed ' name_signal])
    plot(1:length(orig_signal), orig_signal,'LineWidth',1)
    hold on;
    plot(1:length(compressed_signal), compressed_signal,'LineWidth',1)
    hold off;
    xlim([0,length(orig_signal)]);
    title(['Comparison between the Original ' name_signal ' and Compressed Signal '  ': ' wavelet ])
    xlabel('Samples (n)')
    ylabel('Amplitude');
    legend(['Original ' name_signal], 'Compressed Signal')
end