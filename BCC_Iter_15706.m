clear;

% Datos anteriores necesarios.

N_celda_max = 3; % N_celda_max(3) = N_celda_lado(4) - 1.
L = N_celda_max + 1;
Num_atomo_T0 = (N_celda_max + 1)^3 * 2;
Num_par_vec = 1+2*(L-1)*3+4*3*(L-1)^2+8*(L-1)^3;

Prop_1 = 34/128;
Prop_2 = 30/128;
Prop_3 = 24/128;
Prop_4 = 22/128;
Prop_5 = 18/128;

Radio_1 = sqrt(3)*3.16e-10/4; % a = 3.16 (W 1e-10 m)
Radio_2 = sqrt(3)*3.30e-10/4; % a = 3.30 (Nb 1e-10 m)
Radio_3 = sqrt(3)*3.15e-10/4; % a = 3.15 (Mo 1e-10 m)
Radio_4 = sqrt(3)*2.88e-10/4; % a = 2.88 (Cr 1e-10 m)
Radio_5 = sqrt(3)*3.03e-10/4; % a = 3.03 (V 1e-10 m)

A_1 = 12.56e-27; % A = 12.56 1e-27 Jm
B_1 = 16.19e-37; % B = 16.19 1e-37 Jm^2
A_2 = 7.87e-27; % A = 7.87 1e-27 Jm
B_2 = 11.24e-37; % B = 11.24 1e-37 Jm^2
A_3 = 8.86e-27; % A = 8.86 1e-27 Jm
B_3 = 12.09e-37; % B = 12.09 1e-37 Jm^2
A_4 = 4.29e-27; % A = 4.29 1e-27 Jm
B_4 = 5.35e-37; % B = 5.35 1e-37 Jm^2
A_5 = 5.32e-27; % A = 5.32 1e-27 Jm
B_5 = 6.98e-37; % B = 6.98 1e-37 Jm^2

A_promedio = Prop_1 * A_1 + Prop_2 * A_2 + Prop_3 * A_3 + Prop_4 * A_4 + Prop_5 * A_5;
B_promedio = Prop_1 * B_1 + Prop_2 * B_2 + Prop_3 * B_3 + Prop_4 * B_4 + Prop_5 * B_5;
r0_min = 2*B_promedio/A_promedio;
U0_min = - A_promedio/r0_min + B_promedio/(r0_min^2);
Parmtr_ret = 2 * r0_min / sqrt(3);

% Cargamos datos de BCC_sin_Iter.
load("1.5706Info_atom.mat");
load("1.5706Info_pareja_ampl.mat");
U_T = sum(Info_pareja_ampl(:,11))/(Num_par_vec + 1);

% Buscar las posiciones de menor energía.
Orden_atomo_aleat = randperm(Num_atomo_T0);
Info_cambio = zeros (Num_atomo_T0, 9); % 8 columnas: Número del átomo, Tipo del átomo, phi_min, theta_min, rho_min, [x, y, z] finales, U_i.
Iter = 50;
U_iter = zeros (Iter, 4);
dr_2 = zeros (Iter, 2); % Matriz que recoge la desviación de las posición de los átomos de cada iteración respecto a la inicial.
rho_sat = zeros (Iter, 3); % Matriz para comprobar si rho satura.

for Iteracion = 1:Iter

    rho_promd = 0;
    U = 0;

    for i = Orden_atomo_aleat

        Min_U_i = Inf;
        A_i = Info_atomo (i, 3);
        B_i = Info_atomo (i, 4);

        if Iteracion == 1

            x_0 = Info_atomo (i, 5) * Parmtr_ret;
            y_0 = Info_atomo (i, 6) * Parmtr_ret;
            z_0 = Info_atomo (i, 7) * Parmtr_ret;

        elseif Iteracion > 1

            x_0 = Info_cambio (i, 6) * Parmtr_ret;
            y_0 = Info_cambio (i, 7) * Parmtr_ret;
            z_0 = Info_cambio (i, 8) * Parmtr_ret;

        end
    
        % Seleccionamos los vecinos j del átomo i.
        % Creamos una matriz con todas las parejas i-j.
        Colum_interes = 1;
        Valor_interes = i;
        Filas_valor_interes = Info_pareja_ampl(:, Colum_interes) == Valor_interes;
        Filas_selec = Info_pareja_ampl(Filas_valor_interes, :);
        Num_vec_i = size(Filas_selec, 1); % Número de vecinos del átomo i.
    
        for rho = 0 : (0.05e-11 / Iteracion) : (2e-11 / Iteracion)
            for theta = 0 : (pi/36) : pi
                for phi = 0 : (pi/36) : (2*pi)
    
                    U_i = 0;
                    x_esf = rho * sin(theta) * cos(phi);
                    y_esf = rho * sin(theta) * sin(phi);
                    z_esf = rho * cos(theta);
                    x_i = x_0 + x_esf;
                    y_i = y_0 + y_esf;
                    z_i = z_0 + z_esf;
                
                    for j = 1 : Num_vec_i
    
                        Vec_j = Filas_selec (j, 2);
                        A_j = Info_atomo (Vec_j, 3);
                        B_j = Info_atomo (Vec_j, 4);
                        A_ij = (A_i + A_j)/2;
                        B_ij = (B_i + B_j)/2;
    
                        if (Info_cambio (Vec_j, 6) == 0) && (Info_cambio (Vec_j, 7) == 0) && (Info_cambio (Vec_j, 8) == 0)
                            x_j = Info_atomo(Vec_j,5) * Parmtr_ret;
                            y_j = Info_atomo(Vec_j,6) * Parmtr_ret;
                            z_j = Info_atomo(Vec_j,7) * Parmtr_ret;
                        else
                            x_j = Info_cambio(Vec_j,6) * Parmtr_ret;
                            y_j = Info_cambio(Vec_j,7) * Parmtr_ret;
                            z_j = Info_cambio(Vec_j,8) * Parmtr_ret;
                        end
        
                        r_ij = sqrt((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2);
                        U_ij = (- A_ij/r_ij + B_ij/(r_ij^2)) / Num_vec_i; % Para que cada pareja contribuya según su número de vecinos.
                        U_i = U_i + U_ij;
    
                    end
                
                    if  U_i < Min_U_i

                        Min_U_i = U_i;
                        phi_min = phi;
                        theta_min = theta;
                        rho_min = rho;

                    end
                
                end 
    
            end

        end

        x_min = x_0 + (rho_min * sin(theta_min) * cos(phi_min));
        y_min = y_0 + (rho_min * sin(theta_min) * sin(phi_min));
        z_min = z_0 + (rho_min * cos(theta_min));
        Info_cambio (i, :) = [i, Info_atomo(i, 2), phi_min * 180/pi , theta_min * 180/pi, rho_min, x_min/Parmtr_ret, y_min/Parmtr_ret, z_min/Parmtr_ret, U_i];
        Nombre_colum_cambio = {'Número átomo', 'Tipo átomo', 'Phi mínimo', 'Theta mínimo', 'Rho mínimo', 'X nuevo', 'Y nuevo', 'Z nuevo', 'Energía del átomo'};
        Matriz_info_cambio = [Nombre_colum_cambio; num2cell(Info_cambio)];
        U = U + U_i;
        dU = U/(Num_atomo_T0) - U_T;
        rho_promd = rho_promd + rho_min;

    end
    
    dr_2_Iter = (sum((Info_cambio(:,6) - Info_atomo(:,5)).^2 + (Info_cambio(:,7) - Info_atomo(:,6)).^2 + (Info_cambio(:,8) - Info_atomo(:,7)).^2)) * (2/sqrt(3)) / (Num_atomo_T0);
    U_iter (Iteracion,:) = [Iteracion, U/(Num_atomo_T0), U/(Num_atomo_T0)/1.602e-19, dU/1.602e-19];
    dr_2 (Iteracion, :) = [Iteracion, dr_2_Iter];
    rho_sat (Iteracion, :) = [Iteracion, 2e-11 / Iteracion, rho_promd / Num_atomo_T0];

end