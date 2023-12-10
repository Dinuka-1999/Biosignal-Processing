function MSEk = CalcMSEk(N,k,template,epochs)
    yk=mean(epochs(:,1:k),2);
    y=sum((template-yk).^2);
    MSEk=(y/N)^0.5;
end