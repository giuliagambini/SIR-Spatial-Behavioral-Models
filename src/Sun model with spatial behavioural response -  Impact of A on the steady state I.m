%% Figure 4 - Impact of the spatial distancing parameter A
close all; clc; clear;

% Parameters
beta = 35; mu = 1; nu = 1.8;
coeff = [beta*(mu+nu), -beta*mu, mu*(mu+nu)];
I_ee = max(roots(coeff));
S_ee = (mu+nu)/(beta*I_ee);

% Diffusion
D_I = 1; 
D_S = 10; 

% Evaluation of the critical value of A
J11 = -mu - beta*I_ee^2;
J12 = -2*beta*S_ee*I_ee;
J21 = beta*I_ee^2;
J22 = 2*beta*S_ee*I_ee - (mu+nu);

b0_0 = J11*J22 - J12*J21;
A_cr = (D_I*J11 + D_S*J22 - 2*sqrt(D_S*D_I*b0_0)) / J21;
fprintf('A_cr evaluated: %.4f (Paper indicates ~7.54)\n', A_cr);

% Paper panels: 0.97 * A_cr (1° and 2°), 1.03 * A_cr (3°)
A_1 = 0.97 * A_cr;
A_2 = 1.03 * A_cr;

% We need to consider div(A_bar * S * grad(I)), where A = A_bar * S_ee
A_bar_1 = A_1 / S_ee;
A_bar_2 = A_2 / S_ee;
A_bar_array = [A_bar_1, A_bar_2];

% Domain and resolution
Lx = 200; Ly = 200;  
Nx = 201; Ny = 201;  
dx = Lx / (Nx-1); 
dy = Ly / (Ny-1);
x_vec = linspace(0, Lx, Nx);
y_vec = linspace(0, Ly, Ny);

rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initializing the matrices
S = zeros(Ny, Nx, 2);
I = zeros(Ny, Nx, 2);

% Apply to both S and I
for k = 1:2
    S(:,:,k) = S_ee * (1 + rho_S);
    I(:,:,k) = I_ee * (1 + rho_I);
end

% Temporal growth 
dt = 0.02;    
t_end = 107;  
n_steps = round(t_end / dt);

iN = [1, 1:Ny-1];  iS = [2:Ny, Ny];
jW = [1, 1:Nx-1];  jE = [2:Nx, Nx];

fprintf('Calcolo in corso per t = %d... \n', t_end);

for step = 1:n_steps
    for k = 1:2
        Sk = S(:,:,k);
        Ik = I(:,:,k);
        
        % Laplacian
        lapS = (Sk(iN,:) + Sk(iS,:) + Sk(:,jW) + Sk(:,jE) - 4*Sk) / dx^2;
        lapI = (Ik(iN,:) + Ik(iS,:) + Ik(:,jW) + Ik(:,jE) - 4*Ik) / dx^2;
        
        % Flussi per la Cross-Diffusione
        S_E = (Sk + Sk(:, jE))/2;  S_W = (Sk + Sk(:, jW))/2;
        S_S = (Sk + Sk(iS, :))/2;  S_N = (Sk + Sk(iN, :))/2;
        
        dI_E = (Ik(:, jE) - Ik) / dx;  dI_W = (Ik - Ik(:, jW)) / dx;
        dI_S = (Ik(iS, :) - Ik) / dy;  dI_N = (Ik - Ik(iN, :)) / dy;
        
        cross_diff = A_bar_array(k) * ( (S_E .* dI_E - S_W .* dI_W)/dx + (S_S .* dI_S - S_N .* dI_N)/dy );
        
        % Reaction part
        ReacS = mu*(1 - Sk) - beta * (Ik.^2) .* Sk;
        ReacI = beta * (Ik.^2) .* Sk - (mu + nu) * Ik;
        
        % Sum
        S(:,:,k) = Sk + dt * (D_S * lapS + cross_diff + ReacS);
        I(:,:,k) = Ik + dt * (D_I * lapI + ReacI);
    end
end

%% PLOT 
figure(40); clf;
set(gcf, 'Position', [100, 200, 1400, 400]);

% Left panel: 2D map for A = 0.97 Acr 
subplot(1, 3, 1);
imagesc(x_vec, y_vec, I(:,:,1));
set(gca, 'YDir', 'normal'); 
shading interp; colormap parula; colorbar; axis square;
title(sprintf('I(x,y) for A \\approx 0.97 A_{cr}'));
xlabel('x $\rightarrow$', 'Interpreter', 'latex');
ylabel('y $\rightarrow$', 'Interpreter', 'latex');
caxis([0.0, 0.35]);

% Central panel: Mesh 3D per A = 0.97 Acr
subplot(1, 3, 2);
mesh(x_vec, y_vec, I(:,:,1));
colormap parula; axis square;
title(sprintf('Mesh plot for A \\approx 0.97 A_{cr}'));
xlabel('x $\rightarrow$', 'Interpreter', 'latex');
ylabel('y $\rightarrow$', 'Interpreter', 'latex');
zlim([0 1.0]); 
view(-37.5, 30);

% Right panel: Mesh 3D for A = 1.03 Acr 
subplot(1, 3, 3);
mesh(x_vec, y_vec, I(:,:,2));
colormap parula; axis square;
title(sprintf('Mesh plot for A \\approx 1.03 A_{cr} (Suppressed)'));
xlabel('x $\rightarrow$', 'Interpreter', 'latex');
ylabel('y $\rightarrow$', 'Interpreter', 'latex');
zlim([0 1.0]); 
view(-37.5, 30);

fprintf('Done! :) \n');
sound(sin(1:3000));