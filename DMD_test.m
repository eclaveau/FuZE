clear all
close all
clc

b_th180 = [
    0.2047    0.2039    0.2040    0.2046    0.2058    0.2060    0.2056    0.2053    0.2056;
    0.2373    0.2362    0.2355    0.2347    0.2340    0.2329    0.2316    0.2300    0.2286;
    0.2698    0.2686    0.2670    0.2648    0.2623    0.2598    0.2575    0.2548    0.2517;
    0.2406    0.2496    0.2565    0.2593    0.2598    0.2567    0.2513    0.2472    0.2445;
    0.1409    0.1500    0.1637    0.1802    0.1989    0.2227    0.2420    0.2586    0.2706;
    0.1242    0.1200    0.1164    0.1127    0.1071    0.1046    0.1083    0.1163    0.1308;
    0.1311    0.1292    0.1279    0.1252    0.1225    0.1198    0.1170    0.1156    0.1156;
    0.0909    0.0909    0.0909    0.0909    0.0909    0.0881    0.0862    0.0862    0.0844;
    0.1412    0.1399    0.1372    0.1335    0.1301    0.1272    0.1247    0.1216    0.1185];

 % BY: ROY TAYLOR, UNIVERSITY OF WASHINGTON // rktaylor@uw.edu
  % Performs the dynamic mode decomposition.
  %
  % PARAMETERS
  % ----------
  % data : matrix
  % r : int                 (truncate after r modes)
  %                         (FOR NO SVD TRUNCATION, call as [] `empty`)
  % t : array               (timebase)
  % fs : int                (how many extra steps for future-state)
  %
  % RETURNS
  % ----------
  % mud : matrix            (u_modes array)
  % kip : matrix            (phi array)
  % omg : array             (omega; frequency/growth/decay constants)
  % sig : array             (sigma; energy weights per mode)
  %
  % NOTES
  % ----------
  % To reconstruct fully, call
  %      >> reconstruction = mud * kip
  % To truncate for modes a:b, call
  %      >> truncd = mud(:, a:b) * kip(a:b, :)
  % Reconstructions often include spurious imaginary components.
  % Only th`e real component matters.

  b_th180_T = b_th180;
  rank = 8;
  t=0:9;
  
  
  
X1 = b_th180_T(:,1:end-1);
X2 = b_th180_T(:,2:end);
u = b_th180_T(:,1);
dt = 1;

[U1, Sigma1, V1] = svd(X1, 'econ');
U = U1(:,1:rank);
V = V1(:,1:rank);
Sigma = Sigma1(1:rank,1:rank);
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/dt;
Phi = U*eV;
y0 = Phi\u;

for iter=1:(size(b_th180_T,2))
    u_modes(:,iter) = (y0.*exp(omega*t(iter)));

end
  reconstruction = Phi*u_modes;
  figure(1)
  surf(real(reconstruction))