function C = Create_C(H_max,P,smooth_all)

L = zeros(H_max-2,H_max);
for l=1:H_max-2
    L(l,l:l+2) = [1 -2 1];
end 

C2 = L'*L;
% C = L'*L;

XX_size = H_max*(P);
C = zeros(XX_size,XX_size);
%knoneker product for C
I_p = eye(P);
if smooth_all == 1
    C = kron(I_p,C2);
elseif smooth_all == 0    
    C(1:H_max,1:H_max) = C2;
else 
    error('Unknown smoothing')
end    

