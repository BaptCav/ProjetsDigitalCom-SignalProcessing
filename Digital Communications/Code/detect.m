function bhat = detect(r)
% bhat = detect(r)
%
% Computes the received bits given a received sequence of (phase-corrected)
% QPSK symbols. Gray coding of the individual bits is assumed. Hence, the
% two bits for each symbol can be detected from the real and imaginary
% parts, respectively. The first of the two bits below is output first in
% the bhat-sequence.
%
% Assumed mapping:
%
%  10 x   |   x 00
%         |
%  -------+-------
%         |
%  11 x   |   x 01
%
% Input:
%   r  = sequence of complex-valued QPSK symbols
%
% Output:
%   bhat  = bits {0,1} corresponding to the QPSK symbols

n= length(r);
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

