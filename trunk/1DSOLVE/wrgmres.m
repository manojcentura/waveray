Nh = 499;
h = 1/(Nh+1);
eh = ones(Nh,1);
k0 = 5.01*pi;
A = 1/h^2 * spdiags([-eh 2*eh -eh], -1:1, Nh, Nh) - k0^2 * spdiags(eh, 0, Nh, Nh);
% f = rand(Nh,1);
hmesh = (h:h:1-h)';
f = -24*pi^2 * sin(pi*hmesh);;
xstar = A \ f;

x0 = xstar + 1d-3*rand(Nh,1);
maxiter = 1000;
err = zeros(maxiter,1);
H = 1/k0;
NH = round(1 / H);
H = 1/(NH+1);
Hmesh = (H:H:1-H)';

w0 = 1/(2*(1-cos(k0*H)));
w1 = -cos(k0*H)/(1-cos(k0*H));

eH = ones(NH, 1);
Ahatp = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) - ...
  2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH);

Ahatm = 1/H^2 * spdiags([-eH 2*eH -eH], -1:1, NH, NH) + ...
  2i*k0 * 1 / (2*H) * spdiags([-eH eH], [-1 1], NH, NH);

params = [1d-13 10];

atv = @(x_)A*(x_);

for i = 1 : maxiter
  x = gmresb(x0,f,atv,params);
  x0 = x;
  
  if( mod(i,10) == 0 )
    plot(x-xstar)
    title('before')
    pause  
    % rp = (f - A*x) .* exp(-1i*k0*hmesh);
    rph = (f-A*x) .*exp(-1i*k0*hmesh) ;
    % rpH = interp1( [0;hmesh;1], [0;rph;0], Hmesh, 'linear');
    rpH = spline( [0;hmesh;1], [0;rph;0], Hmesh);
    rpH = w0 * [rpH(2:end);0] + ...
      w0 * [0; rpH(1:end-1)] + ...
      w1 * rpH;
    plot( hmesh, real(rph) )
    pause
    plot( Hmesh, real(rpH) )
    pause

    % vph = x .* exp(-1i*k0*hmesh);
    % vpH = interp1( [0;hmesh;1], [0;vph;0], Hmesh, 'linear');
    vpH = Ahatp \ rpH;
    vph = interp1( [0;Hmesh;1], [0;vpH;0], hmesh, 'linear' );
    % vph = zeros(size(vph));

    rmh = (f - A*x) .* exp(1i*k0*hmesh) ;
    rmH = interp1( [0;hmesh;1], [0;rmh;0], Hmesh, 'linear');
    rmH = w0 * [rmH(2:end);0] + ...
      w0 * [0; rmH(1:end-1)] + ...
      w1 * rmH;
    % vmh = x .* exp(1i*k0*hmesh);
    % vmH = interp1( [0;hmesh;1], [0;vmh;0], Hmesh, 'linear');
    vmH = Ahatm \ rmH;
    vmh = interp1( [0;Hmesh;1], [0;vmH;0], hmesh, 'linear' );
    % vmh = zeros(size(vmh));


    x = x + vph .* exp(1i*k0*hmesh) + vmh .* exp(-1i*k0*hmesh);

    % plot( Hmesh, real(vpH) )
    % pause
    % plot( Hmesh, real(vmH) )
    % pause
    plot(real(x-xstar))
    title('after')
    pause  
  end
  % err(i) = norm(x - xstar);
end
semilogy(err)
