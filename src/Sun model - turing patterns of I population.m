%% Figure 2: Turing patterns dependence on D_S for D_I=1
% left panel D_S=3.5, central panel D_S=5, right panel D_S=10. 

close all; clc; 
clearvars -except pphome;
lx = 50; ly = 50; 
nx = 101; % nx s.t. dx=lx/(nx-1)=0.5 as the paper
Ds_values = [3.5, 5.0, 10.0]; 

% Evaluation Endemic Equilibrium
beta = 35; mu = 1; nu = 1.8; % same as the paper, notice that beta is high, so infection rate is really high
coeff = [beta*(mu+nu), -beta*mu, mu*(mu+nu)];
I_ee = max(roots(coeff));
S_ee = (mu+nu)/(beta*I_ee);

figure(20); clf;
set(gcf, 'Position', [100, 200, 1400, 400]);

for i = 1:3
    fprintf('=== Evaluation PDE Ds = %.1f ===\n', Ds_values(i));
    
    par_init = [beta, mu, nu, 2, 0, 0, Ds_values(i), 1];
    p = [];
    p = SIRinit(p, [lx, ly], nx, par_init);
    close(1); close(2); close(6); figure(20);
    
    u0 = zeros(p.nu, 1);
    
    % fix the seed s.t. the values of rho_S and rho_I are the same for all
    % the three panels
    rng(42); 
    
    % Uniform random noise in the interval (-0.01, 0.01)
    rho_S = -0.01 + 0.02 * rand(p.np, 1);
    rho_I = -0.01 + 0.02 * rand(p.np, 1);
    
    % Initial conditions
    u0(1:p.np) = S_ee * (1 + rho_S); 
    u0(p.np+1:2*p.np) = I_ee * (1 + rho_I)
    
    % ODE
    M_single = p.mat.M(1:p.np, 1:p.np);
    Jpat = kron(ones(2), M_single ~= 0); 
    
    opts = odeset('Mass', p.mat.M, 'JPattern', Jpat, 'RelTol', 1e-3, 'AbsTol', 1e-4);
    tspan = [0, 45]; % already at 45 it has nearly reached the equilibrium
    
    odefun = @(t, u) -p.fuha.sG(p, [u; par_init']);
    
    % Resolution of the system
    [~, U] = ode15s(odefun, tspan, u0, opts);
    
    % Extraction final results for I
    u_final = U(end, :)';
    I_final = u_final(p.np+1:2*p.np);
    
    % Plot
    subplot(1, 3, i);
    trisurf(p.mesh.bt(1:3,:)', p.mesh.bp(1,:), p.mesh.bp(2,:), I_final, 'EdgeColor', 'none');
    view(2); shading interp; colormap parula; colorbar; axis square;
    
    caxis([min(I_final), max(I_final)]); 
    title(sprintf('D_S = %.1f', Ds_values(i)));
    xlabel('$x \rightarrow$', 'Interpreter', 'latex'); 
    ylabel('$y \rightarrow$', 'Interpreter', 'latex');
    
    xticks([-50 -25 0 25 50]); xticklabels({'0','25','50','75','100'});
    yticks([-50 -25 0 25 50]); yticklabels({'0','25','50','75','100'});
    
    drawnow;
end
fprintf('\n=== Done!! Turing patterns available :) ===\n');