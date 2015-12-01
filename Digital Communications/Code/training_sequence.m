function b = training_sequence(nr_training_bits)
% b = training_sequence(nr_training_bits)
%
% Generate a training sequence consisting of n bits. Currently, a random
% sequence is used, but a deterministic sequence with better
% autocorrelation properties could be used instead.
%
% Input
%   nr_training_bits = length of training sequence (should be an even number)
%
% Output
%   b = training sequence {0/1}
%
% History:
%   2000-09-28  written /Stefan Parkvall
%   2001-10-21  modified / George J�ngren

%b = (rand(1, nr_training_bits) > .5);


%Finding a good training sequence
number_test = 100;
min_corr=10000; %initial value
L=nr_training_bits/2;
for k = 1 :number_test
    s=(rand(1, 2*L) > .5); %create rand
    s_qpsk=qpsk(s); 
    
    %Calculate autocorrelation difference from a Dirac
    corr=0;
    for i=2:(L-1)
        corr=corr+abs(sum(conj(s_qpsk(i:L)).*s_qpsk(1:(L+1-i)),2)); 
    end
    
    %Find minimum
    if(corr<min_corr)
        min_corr=corr;
        b=s;
    end 
    
end