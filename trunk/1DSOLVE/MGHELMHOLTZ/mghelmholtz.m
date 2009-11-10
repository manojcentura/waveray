function uh = mghelmholtz(Nh, h, k0, HMAX, u0, f, params)
global xstar;
k0round = round(k0/pi)*pi;
uh = u0;
eh = ones(Nh,1);
hmesh = (h:h:1-h)';
Ah = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - ...
  k0^2 * spdiags(eh, 0, Nh, Nh);
atv = @(x_)Ah*(x_);

uh = gmresb(u0, f, atv, params);

rh = f - Ah * uh;

% rH = I_h^H rh
H = h * 2;
NH = (Nh-1) / 2;
rH = 1/4*(rh(1:2:end-2)+rh(3:2:end)+2*rh(2:2:end-1));


if( H >= HMAX )
  HMAX
  2/k0round
  Hmesh = (H:H:1-H)';
  eH = ones(NH, 1);
  AH = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) - ...
    k0^2 * spdiags(eH, 0, NH, NH); 
  w0 = 1/(2*(1-cos(k0round*H)));
  w1 = -cos(k0*H)/(1-cos(k0round*H));

  Ahatp = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) - ...
    2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH) + ...
    (k0round^2-k0^2) * spdiags(eH, 0, NH, NH);

  Ahatm = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) + ...
    2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH) + ...
    (k0round^2-k0^2) * spdiags(eH, 0, NH, NH);

  rpH = rH .* exp(-1i*k0round*Hmesh);
  rpH = w0 * [rpH(2:end);0] + ...
    w0 * [0; rpH(1:end-1)] + ...
    w1 * rpH;
  vpH = Ahatp \ rpH;

  rmH = rH .* exp(+1i*k0round*Hmesh);
  rmH = w0 * [rmH(2:end);0] + ...
    w0 * [0; rmH(1:end-1)] + ...
    w1 * rmH;
  vmH = Ahatm \ rmH;

  % plot(Hmesh, real(rH))
  % pause
  % plot(Hmesh, real(vmH))
  % pause
  deltaH = real(vpH .* exp(1i*k0round*Hmesh) + ...
    vmH .* exp(-1i*k0round*Hmesh));
  deltaH = zeros(NH,1);
else
  deltaH = mghelmholtz(NH, H, k0, HMAX, zeros(NH,1), rH, params);
end

% deltah = I_H^h deltaH
deltah = zeros(Nh,1);
deltah(1) = 0.5*deltaH(1);
deltah(3:2:end-2) = 0.5*(deltaH(1:end-1)+deltaH(2:end));
deltah(end) = 0.5*deltaH(end);
deltah(2:2:end-1) = deltaH;

% clf
% subplot(1,2,1);
% plot(hmesh, rh)
% subplot(1,2,2);
% plot(hmesh, deltah)
% pause

uh = uh + deltah;

uh = gmresb(u0, f, atv, params);
% disp(HMAX/h);
% plot(uh,'r')
% pause
