function IRF = irf_function(A,B_0_inv,str)

[n,m] = size(A);
if m ~= n
    A = [A; [eye(m-n) zeros(m-n,n)]];;
end
K = length(B_0_inv);

%Calculating Phi
IRF = nan(K,K,str.H_max-str.H_min+1);
for i =str.H_min:str.H_max
    First_stage  = A^i;
    IRF(:,:,i+1-str.H_min) = First_stage(1:K,1:K)*B_0_inv;
end

end