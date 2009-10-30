N = 501;
h = 1/(N-1);
e = ones(N,1);
k0 = 100*pi;
A = 1/h^2 * spdiags([-e 2*e -e], -1:1, N, N) - k0^2 * spdiags(e, 0, N, N);
% f = rand(N,1);
mesh = (0:h:1)';
% f = cos(4*pi*mesh);
f = exp( -((mesh-1/2)/(1/16)).^2 );
g = A' * f;
B = A' * A;
L = tril(B, 0);
U = triu(B, 1);
xstar = B \ g;
Linv = inv(L);

x = rand(N,1);
maxiter = 1000;
err = zeros(maxiter,1);
Nh = floor(pi/2.5/k0 / h);
% % hcoarse = Nh * h;
% hcoarse = pi/1.5/k0;
meshcoarse = (0:hcoarse:1)';

w0 = 1/(2*(1-cos(k0*hcoarse)));
w1 = -cos(k0*hcoarse)/(1-cos(k0*hcoarse));

for i = 1 : maxiter
  x = Linv * (g - U*x);
  ferr = abs(fft(x-xstar));
  Weightx = Wx( (x-xstar) .* cos(k0*mesh), Nh, w0, w1 )
  % Weightx = w0 * [interperr(2:end);0] + ...
    % w0 * [0; interperr(1:end-1)] + ...
    % w1 * interperr;
  % plot( ferr(1:100)  )
  % plot( (x-xstar) .* cos(k0*mesh) )
  % plot( (x-xstar) .* sin(k0*mesh) )
  % plot( (x-xstar) )
  % plot( abs(fft(Weightx)) )
  plot( Weightx )
  pause
  err(i) = norm(x - xstar);
end

