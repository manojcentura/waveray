Nh = 499;
h = 1/(Nh+1);
eh = ones(Nh,1);
k0 = 10*pi;
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - k0^2 * spdiags(eh, 0, Nh, Nh);
% f = rand(Nh,1);
hmesh = (h:h:1-h)';
f = exp( -((hmesh-1/2)/(1/16)).^2 );
g = A' * f;
B = A' * A;
L = tril(B, 0);
U = triu(B, 1);
xstar = B \ g;
Linv = inv(L);

x = rand(Nh,1);
maxiter = 1000000;
err = zeros(maxiter,1);
H = pi/3.5/k0;
NH = round(1 / H);
H = 1/(NH+1);
Hmesh = (H:H:1-H)';

w0 = 1/(2*(1-cos(k0*H)));
w1 = -cos(k0*H)/(1-cos(k0*H));

eH = ones(NH, 1);
Ahatp = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) - ...
  2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH);
Bhatp = Ahatp' * Ahatp;
Lhatp = tril(Bhatp);
Uhatp = triu(Bhatp, 1);
Lhatpinv = inv(Lhatp);

Ahatm = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) + ...
  2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH);
Bhatm = Ahatm' * Ahatm;
Lhatm = tril(Bhatm);
Uhatm = triu(Bhatm, 1);
Lhatminv = inv(Lhatm);

for i = 1 : maxiter
  x = Linv * (g - U*x);
  % ferr = abs(fft(x-xstar));
  
  if( mod(i,10000) == 0 )
    % rp = (f - A*x) .* exp(-1i*k0*hmesh);
    rph = (f-A*x) .*exp(-1i*k0*hmesh);
    rpH = interp1( [0;hmesh;1], [0;rph;0], Hmesh, 'linear');
    rpH = w0 * [rpH(2:end);0] + ...
      w0 * [0; rpH(1:end-1)] + ...
      w1 * rpH;
    vph = x .* exp(-1i*k0*hmesh);
    vpH = interp1( [0;hmesh;1], [0;vph;0], Hmesh, 'linear');
    gpH = Ahatp' * rpH;
    for j = 1 : 50
      vpH = Lhatpinv * (gpH - Uhatp * vpH);
    end
    vph = interp1( [0;Hmesh;1], [0;vpH;0], hmesh, 'linear' );
    vph = zeros(size(vph));

    rmh = (f - A*x) .* exp(1i*k0*hmesh);
    rmH = interp1( [0;hmesh;1], [0;rmh;0], Hmesh, 'linear');
    rmH = w0 * [rmH(2:end);0] + ...
      w0 * [0; rmH(1:end-1)] + ...
      w1 * rmH;
    vmh = x .* exp(1i*k0*hmesh);
    vmH = interp1( [0;hmesh;1], [0;vmh;0], Hmesh, 'linear');
    gmH = Ahatm' * rmH;
    for j = 1 : 50
      vmH = Lhatminv * (gmH - Uhatm * vmH);
    end
    vmh = interp1( [0;Hmesh;1], [0;vmH;0], hmesh, 'linear' );
    vmh = zeros(size(vmh));


    x = x + vph .* exp(1i*k0*hmesh) + vmh .* exp(-1i*k0*hmesh);
    % 
    % x = x + vph .* exp(1i*k0*hmesh);

    % plot( Hmesh, real(vpH) )
    % pause
    % plot( Hmesh, real(vmH) )
    % pause
    plot(abs(x-xstar))
    pause  
  end
  % if( mod(i,100) == 0 )
    % plot( abs(Weightx) )
    % pause
  % end
  % err(i) = norm(x - xstar);
end

