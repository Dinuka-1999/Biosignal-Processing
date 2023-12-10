function MSE=CalcMSE(N,nECG,ECG_template)
    b=ones(1,N+1)*(1/(N+1));
    filtered=filter(b,(1),nECG);
    g_delay=floor(N/2);
    filtered=[filtered(g_delay+1:end), ones(1,g_delay).*filtered(end)];
    MSE=sum(((filtered-ECG_template).^2))/length(ECG_template);
end