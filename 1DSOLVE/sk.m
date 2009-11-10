m = 10;
n = m;

A = normrnd(0,1,m,n);

b = normrnd(0,1,m,1);

xstar = A\b;

x = normrnd(0,1,n,1);
xnew = x;
maxiter = 100;
err = zeros(maxiter, 1);
for i = 1 : m
  A1(i,:) = A(i,:) / norm(A(i,:),2);
  b1(i) = b(i) / norm(A(i,:),2);
end

for k = 1 : maxiter
  i = mod(k-1, m)+1;
  % i = ceil(rand(1)*m);
  xnew = x + (b1(i) - A1(i,:)*x) * A1(i,:).';
  x = xnew;
  err(k) = norm(A*x-b,2);
  if( mod(k,1000) == 0 )
    k
    % plot(x-xstar)
    % pause
  end
end

semilogy(1:maxiter, err);
