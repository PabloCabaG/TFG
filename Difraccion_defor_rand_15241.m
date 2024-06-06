clear;

%% Difracción sistema 1 %%

% Cargamos datos necesarios.
load ("1.5241Info_atom.mat");
% load ("1.5241Info_cambio.mat");

% Datos iniciales.
Parmtr_ret = 3.08e-10;
lambda = 1.54e-10; % Longitud de onda.
k = 2 * pi / lambda; % Módulo del vector de onda.
alfa = 1.147; % Ángulo incidencia 65.74 grados.
theta_0 = 20.7 * pi/180; % Mitad del ángulo entre rayo indicente y difractado situación inicial 20.7 grados.
k_x = k * cos(alfa); % Componente x de k.
k_y = - k * sin(alfa); % Componente y de k.
k_z = 0; % Plano (1 1 0)
N = size(Info_atomo,1); % Número de atomos.


% Creamos bucle para la difracción de cada átomo no relajado.
Desv = 1 * pi/180;
Paso = 0.01 * pi/180;
N_paso_T = 2 * Desv / Paso + 1;
Matriz_0 = zeros(N_paso_T, 4);
theta_valores = (theta_0 - Desv) : Paso : (theta_0 + Desv); % Valores de theta para cada paso
Rep = 30;

for m = 1 : N_paso_T

    theta = theta_valores(m);
    N_paso = m; % Utiliza el índice del bucle como N_paso
    F_0 = 0;

    for x = 0 : 4 : 4*Rep
        for y = 0 : 4 : 4*Rep
            for i = 1 : N

                k_x_prima_0 = k * cos(alfa - 2 * theta);
                k_y_prima_0 = - k * sin(alfa - 2 * theta);
                dk_x_0 = k_x_prima_0 - k_x;
                dk_y_0 = k_y_prima_0 - k_y;
                x_0 = (Info_atomo (i,5) + x) * Parmtr_ret;
                y_0 = (Info_atomo (i,6) + y) * Parmtr_ret;
                z_0 = Info_atomo (i,7);
                expnt_0 = dk_x_0 * x_0 + dk_y_0 * y_0;

                if Info_atomo (i,2) == 1 % W, Z = 74
                    f_i = 74;
                elseif Info_atomo (i,2) == 2 % Nb, Z = 41
                    f_i = 41;
                elseif Info_atomo (i,2) == 3 % Mo, Z = 42
                    f_i = 42;
                elseif Info_atomo (i,2) == 4 % Cr, Z = 24
                    f_i = 24;
                elseif Info_atomo (i,2) == 5 % V, Z = 23
                    f_i = 23;
                end

                F_i = f_i * exp(1i * expnt_0);
                F_0 = F_0 + F_i;
                Matriz_0 (N_paso, :) = [N_paso, theta, F_0, abs(F_0)^2];

            end
        end
    end
end


% Creamos bucle para la difracción de cada átomo relajado.

Matriz_relj = zeros(N_paso_T, 4);

e_paso = 0.07;
e_max = 0.56;
dH_frente_e = zeros (e_max/e_paso, 5);
H_0 = max(Matriz_0(:,4));
Random_mod = zeros (N, 4);

for e = 0 : e_paso : e_max

    for j = 1 : N
        Random_mod (j,1) = Info_atomo (j,5);
        Random_mod (j,2) = Info_atomo (j,5) + e * (2 * rand() - 1);
        Random_mod (j,3) = Info_atomo (j,6);
        Random_mod (j,4) = Info_atomo (j,6) + e * (2 * rand() - 1);
    end

    for m = 1 : N_paso_T

        theta = theta_valores(m);
        N_paso = m; % Utiliza el índice del bucle como N_paso
        F_rljd = 0;

        for x = 0 : 4 : 4*Rep
            for y = 0 : 4 : 4*Rep
                for j = 1 : N

                    k_x_prima = k * cos(alfa - 2 * theta);
                    k_y_prima = - k * sin(alfa - 2 * theta);
                    dk_x = k_x_prima - k_x;
                    dk_y = k_y_prima - k_y;
                    x_rljd = (Random_mod (j,1) + x) * Parmtr_ret;
                    y_rljd = (Random_mod (j,2) + y) * Parmtr_ret;
                    expnt = dk_x * x_rljd + dk_y * y_rljd;

                    if Info_atomo (j,2) == 1 % W, Z = 74
                        f_j = 74;
                    elseif Info_atomo (j,2) == 2 % Nb, Z = 41
                        f_j = 41;
                    elseif Info_atomo (j,2) == 3 % Mo, Z = 42
                        f_j = 42;
                    elseif Info_atomo (j,2) == 4 % Cr, Z = 24
                        f_j = 24;
                    elseif Info_atomo (j,2) == 5 % V, Z = 23
                        f_j = 23;
                    end

                    F_j = f_j * exp(1i * expnt);
                    F_rljd = F_rljd + F_j;
                    Matriz_relj (N_paso, :) = [N_paso, theta, F_rljd, abs(F_rljd)^2];

                end
            end
        end
    end

    dx_2 = sum((Random_mod (:,2) - Random_mod (:,1)).^2);
    dy_2 = sum((Random_mod (:,4) - Random_mod (:,3)).^2);
    dr_2 = dx_2 + dy_2;
    H_relj = max(Matriz_relj(:,4));
    dH_frente_e(e/e_paso + 1, :) = [e, dr_2 / N, H_0, H_relj, (H_0 - H_relj)/H_0];

end

% plot(Matriz_0(:,2) * 180 / pi, Matriz_0(:,4), '-b', "LineWidth", 1.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
% set(gca, 'FontSize', 15); % Tamaño índices de los ejes
% xlabel ('\Theta (\circ)')
% ylabel ('|F|^2')
% title ('Sin relajar, U_0 = -95.14eV')
% grid on
% 
% figure;
% 
% plot(Matriz_relj(:,2) * 180 / pi, Matriz_relj(:,4), '-b', "LineWidth", 1.5)%,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
% set(gca, 'FontSize', 15); % Tamaño índices de los ejes
% xlabel ('\Theta (\circ)')
% ylabel ('|F|^2')
% title ('Relajada, U_0 = -95.14eV')
% grid on
% 
% figure;
% 
% plot(Matriz_0(:,2) * 180 / pi, Matriz_0(:,4), '-b', "LineWidth", 1.5)
% hold on
% plot(Matriz_relj(:,2) * 180 / pi, Matriz_relj(:,4), '-r', "LineWidth", 1.5)
% set(gca, 'FontSize', 15); % Tamaño índices de los ejes
% xlabel ('\Theta (\circ)')
% ylabel ('|F|^2')
% title ('U_0 = -95.14eV')
% legend (' Sin relajar',' Relajada')
% lgd = legend;
% lgd.FontSize = 20;
% grid on

plot(dH_frente_e(:,1),dH_frente_e(:,5), '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('Deformación')
ylabel ('\DeltaH/H_0')

figure;

plot(dH_frente_e(:,2),dH_frente_e(:,5), '-b', "LineWidth", 2.5,"Marker",".","MarkerSize",30,"MarkerEdgeColor",'r','MarkerFaceColor','r')
set(gca, 'FontSize', 15); % Tamaño índices de los ejes
xlabel ('\deltar^2')
ylabel ('\DeltaH/H_0')

