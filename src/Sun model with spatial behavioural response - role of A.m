%% Figure 5 - Morphologic transition of patterns
close all; clc; clear;

% Biological parameters
beta = 35; mu = 1; nu = 1.8;
coeff = [beta*(mu+nu), -beta*mu, mu*(mu+nu)];
I_ee = max(roots(coeff));
S_ee = (mu+nu)/(beta*I_ee);

% Diffusion
D_I = 1; 
D_S = 10; 

% A_cr evaluation
J11 = -mu - beta*I_ee^2;
J12 = -2*beta*S_ee*I_ee;
J21 = beta*I_ee^2;
J22 = 2*beta*S_ee*I_ee - (mu+nu);

b0_0 = J11*J22 - J12*J21;
A_cr = (D_I*J11 + D_S*J22 - 2*sqrt(D_S*D_I*b0_0)) / J21;

% 3 scenarios for the figure
A_1 = 0.30 * A_cr;
A_2 = 0.50 * A_cr;
A_3 = 0.78 * A_cr;

% Array of the parameter A_bar = A / S_ee
A_bar_array = [A_1/S_ee, A_2/S_ee, A_3/S_ee];

% Domain and resolution
Lx = 200; Ly = 200;  
Nx = 200; Ny = 200;  % Delta x = 0.5
dx = Lx / (Nx-1); 
dy = Ly / (Ny-1);
x_vec = linspace(0, Lx, Nx);
y_vec = linspace(0, Ly, Ny);

rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initializing the matrices
S = zeros(Ny, Nx, 3);
I = zeros(Ny, Nx, 3);

% Apply to both S and I
for k = 1:3
    S(:,:,k) = S_ee * (1 + rho_S);
    I(:,:,k) = I_ee * (1 + rho_I);
end

% Temporal growth
dt = 0.02;   % Small temporal step to have no convergence problem
t_end = 100;  % Big enough time
n_steps = round(t_end / dt);

iN = [1, 1:Ny-1];  iS = [2:Ny, Ny];
jW = [1, 1:Nx-1];  jE = [2:Nx, Nx];

fprintf('=== Avvio calcolo Figura 5 (3 Scenari HD in parallelo) ===\n');
fprintf('Elaborazione in corso per t = %d... (Richiederà diversi minuti. Pazienza!)\n', t_end);

for step = 1:n_steps
    for k = 1:3
        Sk = S(:,:,k);
        Ik = I(:,:,k);
        
        % Laplacian
        lapS = (Sk(iN,:) + Sk(iS,:) + Sk(:,jW) + Sk(:,jE) - 4*Sk) / dx^2;
        lapI = (Ik(iN,:) + Ik(iS,:) + Ik(:,jW) + Ik(:,jE) - 4*Ik) / dx^2;
        
        % Cross-Diffusion
        S_E = (Sk + Sk(:, jE))/2;  S_W = (Sk + Sk(:, jW))/2;
        S_S = (Sk + Sk(iS, :))/2;  S_N = (Sk + Sk(iN, :))/2;
        
        dI_E = (Ik(:, jE) - Ik) / dx;  dI_W = (Ik - Ik(:, jW)) / dx;
        dI_S = (Ik(iS, :) - Ik) / dy;  dI_N = (Ik - Ik(iN, :)) / dy;
        
        cross_diff = A_bar_array(k) * ( (S_E .* dI_E - S_W .* dI_W)/dx + (S_S .* dI_S - S_N .* dI_N)/dy );
        
        % Reaction terms
        ReacS = mu*(1 - Sk) - beta * (Ik.^2) .* Sk;
        ReacI = beta * (Ik.^2) .* Sk - (mu + nu) * Ik;
        
        % Increase
        S(:,:,k) = Sk + dt * (D_S * lapS + cross_diff + ReacS);
        I(:,:,k) = Ik + dt * (D_I * lapI + ReacI);
    end
end

% Plot
figure(50); clf;
set(gcf, 'Position', [100, 200, 1400, 400]);

labels = {'0.30', '0.50', '0.78'};

for k = 1:3
    subplot(1, 3, k);
    imagesc(x_vec, y_vec, I(:,:,k));
    set(gca, 'YDir', 'normal'); 
    shading interp; colormap parula; colorbar; axis square;
    title(sprintf('I(x,y) for A \\approx %s A_{cr}', labels{k}));
    xlabel('$x \rightarrow$', 'Interpreter', 'latex');
    ylabel('$y \rightarrow$', 'Interpreter', 'latex');
    
    % Adapt
    I_min = min(min(I(:,:,k)));
    I_max = max(max(I(:,:,k)));
    caxis([I_min, I_max]); 
end

fprintf('=== Figura 5 Completata! Osserva la transizione morfologica! ===\n');
sound(sin(1:3000)); 