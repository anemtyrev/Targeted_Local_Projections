function [D, L, K] = select_matrix(n)

m   = n * (n + 1) / 2;
nsq = n^2;
r   = 1;
a   = 1;
v   = zeros(1, nsq);
cn  = cumsum(n:-1:2);   % [EDITED, 2021-08-04], 10% faster
for i = 1:n
   % v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
   v(r:r + i - 2) = i - n + cn(1:i - 1);   % [EDITED, 2021-08-04]
   r = r + i - 1;
   
   v(r:r + n - i) = a:a + n - i;
   r = r + n - i + 1;
   a = a + n - i + 1;
end
D2 = sparse(1:nsq, v, 1, nsq, m);
D = full(D2);

%% L matrix
T = tril(ones(n)); % Lower triangle of 1's
f = find(T(:)); % Get linear indexes of 1's
k = n*(n+1)/2; % Row size of L
m2 = n*n; % Colunm size of L
L = zeros(m2,k); % Start with L'
x = f + m2*(0:k-1)'; % Linear indexes of the 1's within L'
L(x) = 1; % Put the 1's in place
L = L'; % Now transpose to actual L

%% K matrix
% determine permutation applied by K
A = reshape(1:n*n, n, n);
v = reshape(A', 1, []);

% apply this permutation to the rows (i.e. to each column) of identity matrix
P = eye(n*n);
K = P(v,:);