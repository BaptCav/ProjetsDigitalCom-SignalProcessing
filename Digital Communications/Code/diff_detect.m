function bhat = diff_detect(r_diff)
% Differential QPSK demodulation
n= length(r_diff);

r=zeros(1,n);
r(1)=r_diff(1);
% de-differentiate the symbols
for k = 2:n
    r(k)=r_diff(k)/r_diff(k-1);
end
    
bhat=zeros(1,2*n);
%Simple case mapping based on the decisions regions mentionned above
for i= 1:n
    if(real(r(i))>0 && imag(r(i))>0) %00
        bhat(2*i-1)=0;
        bhat(2*i)=0;
    elseif(real(r(i))<0 && imag(r(i))>0) %10
        bhat(2*i-1)=1;
        bhat(2*i)=0;
    elseif(real(r(i))<0 && imag(r(i))<0) %11
        bhat(2*i-1)=1;
        bhat(2*i)=1;
    else %01
        bhat(2*i-1)=0;
        bhat(2*i)=1;
    end
end

