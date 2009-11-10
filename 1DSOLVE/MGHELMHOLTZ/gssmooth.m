function u = gssmooth( Nh, h, k0, u0, f, Niter )
u = u0;
for i = 1 : Niter
  u(1) = 1/(2-(k0*h)^2)*(u(2)+h*h*f(1));
  for l = 2 : Nh-1
    u(l) = 1/(2-(k0*h)^2)*(u(l-1)+u(l+1)+h*h*f(l));
  end
  u(Nh) = 1/(2-(k0*h)^2)*(u(Nh-1)+h*h*f(Nh));
end
