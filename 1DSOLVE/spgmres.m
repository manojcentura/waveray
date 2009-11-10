Nh = 49;
h = 1/(Nh+1);
eh = ones(Nh,1);
k0 = 5.01*pi;
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - k0^2 * spdiags(eh, 0, Nh, Nh);

hmesh = (h:h:1-h)';
% f = exp( -((hmesh-1/2)/(1/16)).^2 );
f = -24*pi^2 * sin(pi*hmesh);
xstar = A \ f;
plot(hmesh,xstar-sin(pi*hmesh));
pause

x0 = xstar + rand(Nh,1)*1e-4;
x = x0;
maxiter = 1000;
err = zeros(maxiter,1);
params=[1d-13 5];
atv = @(x_)A*(x_);

for k = 1 : maxiter
  x = gmresb(x0,f,atv,params);
  x0 = x;
  % err(k) = norm(x-xstar,2);
  err(k) = norm(A*x-f,2);
  if( mod(k,10) == 0 )
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

