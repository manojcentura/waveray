Nh = 49;
h = 1/(Nh+1);
eh = ones(Nh,1);
hmesh = (h:h:1-h)';
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh); 

f = pi^2 * sin(pi*hmesh);
xstar = A \ f;
plot(hmesh,xstar-sin(pi*hmesh));
pause

x0 = xstar + rand(Nh,1)*1e-4;
x = x0;
maxiter = 100;
err = zeros(maxiter,1);

for k = 1 : maxiter
  x = gssmooth(Nh,h,x0,f,10);
  x0 = x;
  err(k) = norm(x-xstar,2);
  % err(k) = norm(A*x-f,2);
  if( mod(k,1) == 0 )
    k
    clf
    hold on
    % plot(hmesh, x)
    % plot(hmesh, xstar, 'r--');
    plot( x - xstar )
    hold off
    pause
  end
end

semilogy(1:maxiter, err);

