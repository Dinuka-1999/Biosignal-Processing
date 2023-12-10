function estimated_signal = wienerFilteredSignal(signal, weights)
    m = length(weights);
    estimated_signal = zeros(size(signal));
    
    % implementing convolution
    for i = 1: length(signal)- m
        estimated_signal(i) = signal(i:i+m-1)* weights;
    end
    
end