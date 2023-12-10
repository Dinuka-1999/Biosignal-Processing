load ECG_template.mat;
t =linspace(0,length(ECG_template)/500,length(ECG_template));
plot(t,ECG_template);
title('ECG template');
xlabel('time(s)')
ylabel('Amplitude(mV)')

nECG=awgn(ECG_template,5,'measured');
[pxx,f]=periodogram(nECG,[],[],500);
plot(f,10*log10(pxx))

new_signal=[0,0 nECG];
ma3ECG_1=ones(1,length(nECG));
for r=3:length(new_signal)
    value=0;
    for i=0:2
        value=value+new_signal(r-i)/3;
    end
    ma3ECG_1(r-2)=value;
end
Group delay = 2/2 = 1

compensated_signal=[ma3ECG_1(2:end) 0 ];
plot(t,ma3ECG_1,t, nECG)

[pxx1,f1]=periodogram(ma3ECG_1,[],[],500);
plot(f1,10*log10(pxx1),f,10*log10(pxx))

b=ones(1,3)*(1/3);
ma3ECG_2=filter(b,[1],nECG);
plot(t,ma3ECG_2,t,nECG,t,ECG_template)

fvtool(b,[1])

b=ones(1,10)*0.1;
fvtool(b,[1])

