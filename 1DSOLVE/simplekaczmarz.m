Nh = 63;
h = 1/(Nh+1);
eh = ones(Nh,1);
% k0 = 5*pi;
k0 = 2*pi;
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) + k0^2 * spdiags(eh, 0, Nh, Nh);

hmesh = (h:h:1-h)';
% f = exp( -((hmesh-1/2)/(1/16)).^2 );
f = 5*pi^2 * sin(pi*hmesh);
xstar = A \ f;
plot(hmesh,xstar-sin(pi*hmesh));
pause

x = xstar + rand(Nh,1)*1e-4;
maxiter = 10000;
err = zeros(maxiter,1);

B = zeros(size(A));
g = zeros(size(f));
for i = 1 : Nh
  B(i,:) = A(i,:) ./ norm(A(i,:));
  g(i) = f(i) ./ norm(A(i,:)); 
end

for k = 1 : maxiter
  i = mod(k-1, Nh)+1;
  for i = 1 : Nh
    x = x + (g(i) - B(i,:)*x) * B(i,:)';
  end
  err(k) = norm(x-xstar,2);
  if( mod(k,100) == 0 )
    k
    clf
    hold on
    % % plot(hmesh, x)
    % % plot(hmesh, xstar, 'r--');
    plot( x - xstar )
    hold off
    pause
  end
end

semilogy(1:maxiter, err);
