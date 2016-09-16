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


rank = 8;
b_untrans = b_th180;
b_th180_T = b_th180';

 raw = b_th180_T;
  r = rank;
  t = linspace(0,9,9);
  fs = 0;
  
  x = raw(:, 1:end-1);
  y = raw(:, 2:end);
  first_frame = raw(:, 1);

  % finalize our timebase
  dt = mean(diff(t));
  tfinal = t(end) + dt*fs;
  time = [t(1):dt:tfinal];

  % truncated SVD
  [U1, S1, V1] = svd(x, 'econ');

  if isempty(r)==1
    U = U1;
    V = V1;
    S = S1;
  else
    U = U1(:, 1:r);
    V = V1(:, 1:r);
    S = S1(1:r, 1:r);
  end

  % linear map s.t. y = Ax
  Atilde = U'*y*V*diag(1./diag(S));
  [eV, D] = eig(Atilde); mu = diag(D); omg = log(mu)/dt;
  mud = U*eV; y0 = mud\first_frame;

  % compute modes
  for iterate = 2:length(time)-1
    kip(:, iterate) = (y0.*exp(omg*time(iterate)));
  end

  sig = diag(S);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% u_dmd = Phi*u_modes;
u_dmd = mud*kip;
% err_dmd = mean(mean(abs(real(u_dmd)-b_th180_T)));
% X_DMD = Phi*X_inter;
% rel_err_dmd = err_dmd/err_svd;

figure(11)
surf(real(u_dmd))
xlabel('t')
ylabel('x')
title('DMD generated data')
figure(12)
surf(real(u_dmd)-b_th180_T(:,1:end-1))
xlabel('t')
ylabel('x')
title('Error')
