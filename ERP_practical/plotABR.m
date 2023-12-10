function plotABR(Data1,data2,Limit)

    peaks=findPeaks(Data1,Limit);
    j = 0;
    epochs=[];
    
    for i=1:length(peaks)-1
        j = j + 1;
        epochs(:,j) = data2(peaks(i)-50:peaks(i)+149);
    end

    ensmbl_avg = mean(epochs(:,(1:length(peaks)-1)),2);
    %disp(length(ensmbl_avg));
    figure;
    plot((-50:149)/10,ensmbl_avg*(10^6),'LineWidth',1)
    title(['Ensemble averaged ABR from ',num2str(length(epochs)),' epochs'])
    ylim([-15 20]);
    ylabel('Amplitude(uV)');
    xlabel('Time(ms)')
end