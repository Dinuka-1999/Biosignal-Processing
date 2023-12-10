function [y,t]=getLinearModel(start,stop,Fs)
t=(1:stop-start)/Fs;
y=zeros(1,stop-start);
for r=1:stop-start
    if (r<=6)
        y(r)=-0.105;
    elseif (r<=9)
        y(r)=(r-9)*(0.035+0.065)/2+0.035;
    elseif (r<=13)
        y(r)=(r-13)*(0.015+0.095)/(-3)-0.095;
    elseif (r<=26)
        y(r)=(-0.095-0.115)/2;
    elseif (r<=28)
        y(r)=(r-28)*(-0.345+0.375)/(-1)-0.375;
    elseif (r<=31)
        y(r)=(r-31)*(1.805+0.025)/2+1.805;
    elseif (r<=34)
        y(r)=(r-34)*(1.045+0.585)/(-2)-0.585;
    elseif (r<=37)
        y(r)=(r-37)*(-0.105+0.405)/2-0.105;
    elseif (r<=54)
        y(r)=(-0.105-0.045)/2;
    elseif (r<=63)
        y(r)=(r-63)*(0.185+0.005)/8+0.185;
    elseif (r<=69)
        y(r)=(0.205+0.235)/2;
    elseif (r<=77)
        y(r)=(r-77)*(0.185+0.045)/(-7)-0.045;
    else
        y(r)= -0.045;
    end
end
end