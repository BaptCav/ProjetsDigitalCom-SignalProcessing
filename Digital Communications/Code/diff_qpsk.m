function d = diff_qpsk(b)
% Differential QPSK modulation

d=qpsk(b);

for symbol = 2:length(d)  
    d(symbol)=d(symbol)*d(symbol-1);
end