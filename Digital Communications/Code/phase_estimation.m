function phihat = phase_estimation(r, b_train)
% phihat = phase_estimation(r, b_train)
%
% Phase estimator using the training sequence. The phase estimate is
% obtained by minimizing the norm of the difference between the known
% transmitted QPSK-modulated training sequence and the received training
% part. NB! There are other ways of estimating the phase, this is just
% one example.
%
% Input:
%   r       = received baseband signal
%   b_train = the training sequence bits
%
% Output:
%   phihat     = estimated phase
bqpsk=qpsk(b_train);
rt=r(1:length(bqpsk));%% We consider only the bits of the training sequence
%%Initialisation of the things to consider in the for loop
e=dot(rt-bqpsk,rt-bqpsk);
phihat=0;
%%for loop in order to estimate the phase that minimizes the distance
%%beetween the QPSK Mapped training sequence and the first bits received in
%%the synchronized signal
for p=0:0.01:2*pi
    rprime=rt*exp(-j*p);
    eprime=dot(rprime-bqpsk,rprime-bqpsk);
    if(e>eprime)%%If our new error if lower we take the considered phase
        e=eprime;
        phihat=p;
    end
end
     



