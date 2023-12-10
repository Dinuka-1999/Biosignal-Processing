function SGMSE=CalcSGMSE(N,L,nECG,ECG_template)
    filtered1=sgolayfilt(nECG,N,2*L+1);
    SGMSE=sum(((filtered1-ECG_template).^2))/length(ECG_template);
end