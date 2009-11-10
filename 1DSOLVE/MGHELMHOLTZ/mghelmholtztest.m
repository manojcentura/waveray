global xstar;
Nh = 2^9-1;
h = 1/(Nh+1);
eh = ones(Nh,1);
hmesh = (h:h:1-h)';
k0 = 5.1*pi;
k0round = round(k0/pi)*pi;
k0vec = sin(k0round*hmesh);
k0vec = k0vec / norm(k0vec);
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - ...
  k0^2 * spdiags(eh, 0, Nh, Nh);


f = ( pi^2 - k0^2 ) * sin(pi*hmesh);
f = f - (f' * k0vec) * k0vec;
xstar = A \ f;
xstar = xstar - (xstar' * k0vec) * k0vec
plot(hmesh,xstar-sin(pi*hmesh));
pause

x0 = xstar + rand(Nh,1)*1e-2;
x0 = x0 - (x0' * k0vec) * k0vec;
x = x0;
maxiter = 1000;
err = zeros(maxiter,1);
atv = @(x_)A*(x_);
params=[1d-13 50];

for i = 1 : maxiter
  x = mghelmholtz(Nh, h, k0, h*16, x0, f, params);
  x0 = x;
  plot(x-xstar)
  pause
end
