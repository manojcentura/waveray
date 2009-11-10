global xstar;
Nh = 2^10-1;
h = 1/(Nh+1);
eh = ones(Nh,1);
hmesh = (h:h:1-h)';
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh); 

f = pi^2 * sin(pi*hmesh);
xstar = A \ f;
plot(hmesh,xstar-sin(pi*hmesh));
pause

x0 = xstar + rand(Nh,1)*1e-1;
x = x0;
maxiter = 100;
err = zeros(maxiter,1);

Niter = 5;

for i = 1 : 10
  x = mgpoisson(Nh, h, h*16, x0, f, Niter);
  x0 = x;
  plot(x-xstar)
  pause
end
