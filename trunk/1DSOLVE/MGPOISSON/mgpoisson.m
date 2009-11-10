function uh = mgpoisson(Nh, h, HMAX, u0, f, Niter)
global xstar;
uh = u0;
uh = gssmooth(Nh, h, uh, f, Niter);
% rh = f - Ah * u
rh = zeros(Nh,1);
rh(1) = f(1) - 1/h^2*(2*uh(1)-uh(2));
rh(2:end-1) = f(2:end-1) - ...
  1/h^2*(2*uh(2:end-1)-uh(1:end-2)-uh(3:end));
rh(end) - f(end) - 1/h^2*(2*uh(end)-uh(end-1));
% rH = I_h^H rh
H = h * 2;
NH = (Nh-1) / 2;
rH = 1/4*(rh(1:2:end-2)+rh(3:2:end)+2*rh(2:2:end-1));

if( H >= HMAX )
  eH = ones(NH,1);
  AH = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH); 
  deltaH = AH \ rH;
else
  deltaH = mgpoisson(NH, H, HMAX, zeros(NH,1), rH, Niter);
end

% deltah = I_H^h deltaH
deltah = zeros(Nh,1);
deltah(1) = 0.5*deltaH(1);
deltah(3:2:end-2) = 0.5*(deltaH(1:end-1)+deltaH(2:end));
deltah(end) = 0.5*deltaH(end);
deltah(2:2:end-1) = deltaH;

uh = uh + deltah;

uh = gssmooth(Nh, h, uh, f, Niter);
% disp(HMAX/h);
% plot(uh,'r')
% pause
