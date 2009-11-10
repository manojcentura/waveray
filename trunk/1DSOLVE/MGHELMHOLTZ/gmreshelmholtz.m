Nh = 31;
h = 1/(Nh+1);
eh = ones(Nh,1);
hmesh = (h:h:1-h)';
k0 = 5.1*pi;
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - ...
  k0^2 * spdiags(eh, 0, Nh, Nh);

f = (1-k0^2)*pi^2 * sin(pi*hmesh);
xstar = A \ f;
plot(hmesh,xstar-sin(pi*hmesh));
pause

x0 = xstar + rand(Nh,1)*1e-2;
x = x0;
maxiter = 100;
err = zeros(maxiter,1);
atv = @(x_)A*(x_);
params=[1d-13 10];


for k = 1 : maxiter
  x = gmresb(x0,f,atv,params);
  x0 = x;
  err(k) = norm(x-xstar,2);
  % err(k) = norm(A*x-f,2);
  if( mod(k,1) == 0 )
    k
    clf
    hold on
    plot( x - xstar )
    hold off
    pause
  end
end

semilogy(1:maxiter, err);
