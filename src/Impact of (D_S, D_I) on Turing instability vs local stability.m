%% Figure 1 - Spatial Behavioural Responses (d'Onofrio et al.) 
% Impact of (D_S,D_I) on Turing Instability vs Local Stability
clear; clc; close all;

% Parameters
beta = 35;
mu = 1;
nu = 1.8;

% Evaluation of the Endemic Equilibrium
% From the steady state conditions, we solve the quadratic for I:
a_coef = (beta * (mu + nu)) / mu;
b_coef = -beta;
c_coef = mu + nu;

% Find roots and select the appropriate Endemic Equilibrium
I_roots = roots([a_coef, b_coef, c_coef]);
I_e = max(I_roots); % We take the stable EE
S_e = (mu + nu) / (beta * I_e);

% Evaluate the Jacobian at the EE
J11 = -mu - beta * I_e^2;
J12 = -2 * beta * S_e * I_e;
J21 = beta * I_e^2;
J22 = 2 * beta * S_e * I_e - (mu + nu);

detJ = J11 * J22 - J12 * J21; 
trJ = J11 + J22; 

% Set up common spatial frequency axis (k)
k_vec = linspace(0, 2, 500); 

% Plot
figure('Name', 'Figure 1: Turing Instability Regions', 'Position', [100, 100, 1200, 400]);

% Colours
dark_blue = [0.16, 0.32, 0.53];
beige     = [0.96, 0.91, 0.76];
colormap([dark_blue; beige]); 

%% Panel 1: Ratio r = Ds/Di vs k (A = 0)
subplot(1, 3, 1);
r_vec = linspace(0, 10, 500);
Di_1 = 1;
[K1, R] = meshgrid(k_vec, r_vec);
Ds_1 = R .* Di_1; 

a0_1 = Ds_1 .* Di_1 .* K1.^4 - (Di_1 .* J11 + Ds_1 .* J22) .* K1.^2 + detJ;
a1_1 = (Ds_1 + Di_1) .* K1.^2 - trJ;
turing_1 = (a0_1 < 0) & (a1_1 > 0);

% (Uso contourf solo per richiedere che ci siano 2 livelli di riempimento)
contourf(K1, R, double(turing_1), [0 0.5 1], 'LineStyle', 'none');

title('r vs k (A = 0)');
xlabel('k');
ylabel('r = D_S / D_I');

% Linea
hold on;
% Cerco le coordinate del punto più basso della regione di turing
[r_idx, k_idx] = find(turing_1);
[~, pos] = min(r_idx); % Indice della riga più bassa
r_min_val = r_vec(r_idx(pos));
k_min_val = k_vec(k_idx(pos));

% linea da x=0 fino alla k del minimo
plot([0, k_min_val], [r_min_val, r_min_val], '--r', 'LineWidth', 1.5);

% aggiungo il valore all'asse y
ax = gca;
yt = yticks;
% rimozione 'tick' troppo vicini per evitare sovrapposizioni
yt(abs(yt - 3.11) < 0.5) = []; 
ax.YTick = sort([yt, 3.11]);
hold off;

%% PANEL 2: Ds vs k (Limit Case D_I = 0, A = 0)
subplot(1, 3, 2);
Ds_vec = linspace(0, 10, 500);
Di_2 = 0;
[K2, DS2] = meshgrid(k_vec, Ds_vec);

a0_2 = DS2 .* Di_2 .* K2.^4 - (Di_2 .* J11 + DS2 .* J22) .* K2.^2 + detJ;
a1_2 = (DS2 + Di_2) .* K2.^2 - trJ;

turing_2 = (a0_2 < 0) & (a1_2 > 0);

contourf(K2, DS2, double(turing_2), [0 0.5 1], 'LineColor', 'none');
title('D_S vs k (D_I = 0, A = 0)');
xlabel('k');
ylabel('D_S');

%% PANEL 3: Spatial Distancing A vs k (D_S = 10, D_I = 1)
subplot(1, 3, 3);
A_vec = linspace(0, 8, 500);
Ds_3 = 10;
Di_3 = 1;
[K3, A_mat] = meshgrid(k_vec, A_vec);

% Base polynomial (without A)
a0_base = Ds_3 * Di_3 * K3.^4 - (Di_3 * J11 + Ds_3 * J22) * K3.^2 + detJ;
a1_base = (Ds_3 + Di_3) * K3.^2 - trJ;

% Now add the spatial distancing term (A) to the determinant
b0 = a0_base + J21 .* A_mat .* K3.^2;

turing_3 = (b0 < 0) & (a1_base > 0);

contourf(K3, A_mat, double(turing_3), [0 0.5 1], 'LineColor', 'none');
title('A vs k (D_S = 10, D_I = 1)');
xlabel('k');
ylabel('A');

% --- AGGIUNTA LINEA A = 7.54 ---
% --- AGGIUNTA LINEA AL MASSIMO E TICK SULL'ASSE Y ---
hold on;
% Cerchiamo le coordinate del punto più alto della regione di Turing
[A_idx, k_idx] = find(turing_3);
[~, pos] = max(A_idx); % Indice della riga più alta
A_max_val = A_vec(A_idx(pos));
k_max_val = k_vec(k_idx(pos));
% Disegniamo la linea da x=0 fino alla k del massimo
plot([0, k_max_val], [A_max_val, A_max_val], '--r', 'LineWidth', 1.5);

% Aggiungiamo il valore esatto come "tacca" sull'asse 
ax = gca;
yt = yticks;
yt(abs(yt - 7.54) < 0.5) = []; 
ax.YTick = sort([yt, 7.54]);
% -------------------------------

% Legend
hold on;
% 2 squares to put the legend with right colors
h1 = plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', dark_blue, 'MarkerEdgeColor', 'none');
h2 = plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', beige, 'MarkerEdgeColor', 'none');

lgd = legend([h1, h2], {'Stable Region', 'TI Region'}, 'Location', 'northeast');
lgd.FontSize = 10;
hold off;