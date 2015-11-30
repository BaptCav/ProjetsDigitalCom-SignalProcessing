function d = qpsk(b)
% d = qpsk(b)
%
% Map the bits to be transmitted into QPSK symbols using Gray coding. The
% resulting QPSK symbol is complex-valued, where one of the two bits in each
% QPSK symbol affects the real part (I channel) of the symbol and the other
% bit the imaginary part (Q channel). Each part is subsequently PAM
% modulated to form the complex-valued QPSK symbol. The energy per QPSK
% symbol is normalized to unity.
%
% The mapping resulting from the two PAM branches are:
%
% complex part (Q channel)
%         ^
%         |
%  10 x   |   x 00   (odd bit, even bit)
%         |
%  -------+------->  real part (I channel)
%         |
%  11 x   |   x 01
%         |
%
%
%
% Input:
%   b = bits {0, 1} to be mapped into QPSK symbols
%
% Output:
%   d = complex-valued QPSK symbols

%Add a 0 the end if the sequence is not even
if(mod(length(b),2)==1)
    b=[b,0];
end

N= length(b);
d=zeros(1,N/2);
for i =1:(N/2)
    if(b(2*i-1)==1 && b(2*i)==1) %11
        d(i)= -1-1i;
    elseif(b(2*i-1)==0 && b(2*i)==1) %01
        d(i)=1-1i;
    elseif (b(2*i-1)==0 && b(2*i)==0) %00
        d(i)=1+1i;
    else %10
        d(i)=1i-1;
    end
end
d=1/sqrt(2)*d;

