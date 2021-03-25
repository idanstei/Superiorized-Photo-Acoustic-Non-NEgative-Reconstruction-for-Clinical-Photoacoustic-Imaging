function out = Quantization(in,Bits,Utilization)

% Input range
i1 = min(in(:));  
i2 = max(in(:));

% Output range
o1 = ceil(-2^(Bits-1)*Utilization); 
o2 = floor(2^(Bits-1)*Utilization);
 
% Transform input to output range, round and return to the original range
A = (o2-o1)/(i2-i1); B = o1-i1*A;
out = (round(A*in+B)-B)./A;