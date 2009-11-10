function u = gssmooth( Nh, h, u0, f, Niter )
u = u0;
for i = 1 : Niter
  u(1) = 0.5*(u(2)+h*h*f(1));
  for l = 2 : Nh-1
    u(l) = 0.5*(u(l-1)+u(l+1)+h*h*f(l));
  end
  u(Nh) = 0.5*(u(Nh-1)+h*h*f(Nh));
end
