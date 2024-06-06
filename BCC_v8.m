clear;

%% ESTRUCTURA BCC %%

% Número de celda unidad por lado.
N_celda_max = 3; % N_celda_max(3) = N_celda_lado(4) - 1.
L = N_celda_max + 1;

% Proporciones de cada tipo de átomo.
Prop_1 = 52/128;
Prop_2 = 38/128;
Prop_3 = 38/128;

% Recuento número de átomos.
Num_atomo_T0 = 0;
Num_atomo_T = (N_celda_max + 1)^3 * 2; % Recuento inverso.
Num_atomo_1 = Num_atomo_T * Prop_1;
Num_atomo_2 = Num_atomo_T * Prop_2;
Num_atomo_3 = Num_atomo_T * Prop_3;

% Radios.
Radio_1 = sqrt(3)*3.16e-10/4; % a = 3.16 (W 1e-10 m)
Radio_2 = sqrt(3)*3.30e-10/4; % a = 3.30 (Nb 1e-10 m)
Radio_3 = sqrt(3)*3.15e-10/4; % a = 3.15 (Mo 1e-10 m)
       % Tántalo(Ta), Cromo(Cr), Vanadio(V) %

% Constantes A y B para cada tipo de átomo.
A_1 = 12.56e-27; % A = 12.56 1e-27 Jm
B_1 = 16.19e-37; % B = 16.19 1e-37 Jm^2
A_2 = 7.87e-27; % A = 7.87 1e-27 Jm
B_2 = 11.24e-37; % B = 11.24 1e-37 Jm^2
A_3 = 8.86e-27; % A = 8.86 1e-27 Jm
B_3 = 12.09e-37; % B = 12.09 1e-37 Jm^2

% Cálculo de energía antes de las microdeformaciones.
A_promedio = Prop_1 * A_1 + Prop_2 * A_2 + Prop_3 * A_3;
B_promedio = Prop_1 * B_1 + Prop_2 * B_2 + Prop_3 * B_3;
r0 = 0:0.01e-10:3.50e-10;
U0 = (- A_promedio./r0 + B_promedio./(r0.^2)); %*6.022e3; % Esto estaría en unidades de kJ/mol.
hold on
plot (r0, U0, 'LineWidth', 1.5);
r0_min = 2*B_promedio/A_promedio;
U0_min = - A_promedio/r0_min + B_promedio/(r0_min^2);
plot (r0_min,U0_min,'.', 'MarkerSize', 25)
xlabel ('r_0')
ylabel ('U_0')
hold off

% Distancia entre dos átomos consecutivos en las esquinas.
Parmtr_ret = 2 * r0_min / sqrt(3);

% Distancia a los primeros vecinos.
Dist_vecino = r0_min;

% Recuento energía.
Energia_T = 0;

% Matriz donde se almacena la información de los átomos.
Info_atomo = zeros((N_celda_max + 1)^3 * 2, 7); % 7 columnas: Número, Tipo, A, B, X, Y, Z

% Representar la estructura.
figure;
hold on;
axis equal;

% ------------------------------------------------------------------------
% Colocar átomos en las esquinas de cada celda unidad.

for x = 0:Parmtr_ret:N_celda_max*Parmtr_ret
    for y = 0:Parmtr_ret:N_celda_max*Parmtr_ret
        for z = 0:Parmtr_ret:N_celda_max*Parmtr_ret

            % Posición de los átomos.
            Pos_esquina = [x, y, z];
            
            % Genera número aleatorio.
            Num_aleatorio = rand();

            % Determinar el tipo de átomo según las proporciones dadas.
            if Num_aleatorio < Prop_1
                % Tipo 1.
                Tipo_atom = 1;
                Color_esq = 'r';
                Radio_esq = Radio_1;
                A = A_1;
                B = B_1;
                Num_atomo_1 = Num_atomo_1 - 1;
            elseif Num_aleatorio < (Prop_1+Prop_2)
                % Tipo 2.
                Tipo_atom = 2;
                Color_esq = 'g';
                Radio_esq = Radio_2;
                A = A_2;
                B = B_2;
                Num_atomo_2 = Num_atomo_2 - 1;
            else
                % Tipo 3.
                Tipo_atom = 3;
                Color_esq = 'b';
                Radio_esq = Radio_3;
                A = A_3;
                B = B_3;
                Num_atomo_3 = Num_atomo_3 - 1;
            end

            % Cálculo de las proporciones reales.
            Prop_1 = Num_atomo_1/Num_atomo_T;
            Prop_2 = Num_atomo_2/Num_atomo_T;
            Prop_3 = Num_atomo_3/Num_atomo_T;

            % Almacenar información del átomo en la matriz.
            Nombre_colum_atom = {'Número del átomo', 'Tipo del átomo', 'A', 'B', 'X', 'Y', 'Z'};
            Info_atomo(Num_atomo_T0 + 1, :) = [Num_atomo_T0 + 1, Tipo_atom, A, B, Pos_esquina/Parmtr_ret];
            Matriz_info_atom = [Nombre_colum_atom; num2cell(Info_atomo)];       

            % Acumular número de átomos.
            Num_atomo_T0 = Num_atomo_T0 + 1;

            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_esq + Pos_esquina(1);
            y_sphere = y_sphere * Radio_esq + Pos_esquina(2);
            z_sphere = z_sphere * Radio_esq + Pos_esquina(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_esq, 'EdgeColor', 'k');
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T - 1;

        end
    end
end

% ------------------------------------------------------------------------
% Colocar átomos en el centro de cada celda unidad.
for x = 0:Parmtr_ret:N_celda_max*Parmtr_ret
    for y = 0:Parmtr_ret:N_celda_max*Parmtr_ret
        for z = 0:Parmtr_ret:N_celda_max*Parmtr_ret

            % Posición del átomo en el centro de la celda.
            Pos_centro = [x, y, z] + Parmtr_ret * (1/2) * [1, 1, 1];
            
            % Genera número aleatorio.
            Num_aleatorio = rand();

            % Determinar el tipo de átomo según las proporciones dadas.
            if Num_aleatorio < Prop_1
                % Tipo 1.
                Tipo_atom = 1;
                Color_centro = 'r';
                Radio_centro = Radio_1;
                A = A_1;
                B = B_1;
                Num_atomo_1 = Num_atomo_1 - 1;
            elseif Num_aleatorio < (Prop_1+Prop_2)
                % Tipo 2.
                Tipo_atom = 2;
                Color_centro = 'g';
                Radio_centro = Radio_2;
                A = A_2;
                B = B_2;
                Num_atomo_2 = Num_atomo_2 - 1;
            else
                % Tipo 3.
                Tipo_atom = 3;
                Color_centro = 'b';
                Radio_centro = Radio_3;
                A = A_3;
                B = B_3;
                Num_atomo_3 = Num_atomo_3 - 1;
            end
            
            % Cálculo de las proporciones reales.
            Prop_1 = Num_atomo_1/Num_atomo_T;
            Prop_2 = Num_atomo_2/Num_atomo_T;
            Prop_3 = Num_atomo_3/Num_atomo_T;

            % Almacenar la información de cada átomo en la matriz.
            Info_atomo(Num_atomo_T0 + 1, :) = [Num_atomo_T0 + 1, Tipo_atom, A, B, Pos_centro/Parmtr_ret];
            Matriz_info_atom = [Nombre_colum_atom; num2cell(Info_atomo)];

            % Acumular número de átomos.
            Num_atomo_T0 = Num_atomo_T0 + 1;


            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_centro + Pos_centro(1);
            y_sphere = y_sphere * Radio_centro + Pos_centro(2);
            z_sphere = z_sphere * Radio_centro + Pos_centro(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_centro, 'EdgeColor', 'none');
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T - 1;

        end
    end
end

% ------------------------------------------------------------------------
% Sumar energía de cada pareja de vecinos.
U_T = 0;
Num_parejas_1 = 0;
Num_par_vec = 1+2*(L-1)*3+4*3*(L-1)^2+8*(L-1)^3;

% Matriz donde se almacena la información de las parejas de vecinos.
Info_pareja_1 = zeros(Num_par_vec, 11); % 11 columnas: Número del átomo, Número vecino, Tipo átomo, Tipo vecino, Coordenadas átomo, Coordenadas vecino, Energía pareja.

for x = 0:N_celda_max
    for y = 0:N_celda_max
        for z = 0:N_celda_max

            % Localizamos el átomo.
            N = (N_celda_max + 1)^3 - (N_celda_max - x) * (N_celda_max + 1)^2 - (N_celda_max - y) * (N_celda_max + 1) - (N_celda_max - z);

            % Definimos los parámetros del átomo.
            A_N = Info_atomo (N, 3);
            B_N = Info_atomo (N, 4);

            % Localizamos a los vecinos.
            for h= -1/2:1/2
                for k = -1/2:1/2
                    for m = -1/2:1/2

                        % Coordenas del vecino.
                        [x_vec, y_vec, z_vec] = deal(x + h, y + k, z + m);

                        % Número del vecino.
                        N_vec = ((N_celda_max + 1)^3) * 2 - (N_celda_max + 1/2 - x_vec) * (N_celda_max + 1)^2 - (N_celda_max + 1/2 - y_vec) * (N_celda_max + 1) - (N_celda_max + 1/2 - z_vec);

                        % Comprobamos que el vecino exista.
                        if (x_vec >= 1/2 && x_vec <= (N_celda_max + 1/2)) && (y_vec >= 1/2 && y_vec <= (N_celda_max + 1/2)) && (z_vec >= 1/2 && z_vec <= (N_celda_max + 1/2))
                            
                            % Definimos constantes del vecino.
                            A_N_vec = Info_atomo (N_vec, 3);
                            B_N_vec = Info_atomo (N_vec, 4);
                            Num_parejas_1 = Num_parejas_1 + 1; 
                            
                            % Calculamos las constantes promedio de cada pareja de vecinos.
                            A_mn = (A_N + A_N_vec)/2;
                            B_mn = (B_N + B_N_vec)/2;
                            
                            % Calculamos la energía de cada pareja de átomos.
                            U_mn =  (- A_mn/r0_min + B_mn/(r0_min^2))/ Num_par_vec; % Comprobar si se soluciona dividiéndolo entre 2.
                            
                            % Sumamos la energía de cada pareja de átomos.
                            U_T = U_T + U_mn; 
                            
                            % Almacenamos la información de cada pareja de vecinos en la matriz.
                            Nombre_colum_pareja = {'Número del átomo', 'Número vecino', 'Tipo átomo', 'Tipo vecino', 'X', 'Y', 'Z', 'X_VEC', 'Y_VEC', 'Z_VEC', 'Energía Pareja'};
                            Info_pareja_1(Num_parejas_1, :) = [Info_atomo(N,1), Info_atomo(N_vec,1), Info_atomo(N,2), Info_atomo(N_vec,2), x, y, z, x_vec, y_vec, z_vec, U_mn];
                            Matriz_info_pareja = [Nombre_colum_pareja; num2cell(Info_pareja_1)];

                        else
                            continue;
                        end   
                    end
                end
            end
        end
    end
end

% Nuevo bucle para tener ordenados los vecinos de los átomos centrales.
Num_parejas_2 = 0;

% Matriz donde se almacena la información de las parejas de vecinos.
Info_pareja_2 = zeros(Num_par_vec, 11); % 11 columnas: Número del átomo, Número vecino, Tipo átomo, Tipo vecino, Coordenadas átomo, Coordenadas vecino, Energía pareja.

for x = 1/2:(N_celda_max + 1/2)
    for y = 1/2:(N_celda_max + 1/2)
        for z = 1/2:(N_celda_max + 1/2)

            % Localizamos el átomo.
            N = ((N_celda_max + 1)^3) * 2 - (N_celda_max + 1/2 - x) * (N_celda_max + 1)^2 - (N_celda_max + 1/2 - y) * (N_celda_max + 1) - (N_celda_max + 1/2 - z);

            % Definimos los parámetros del átomo.
            A_N = Info_atomo (N, 3);
            B_N = Info_atomo (N, 4);

            % Localizamos a los vecinos.
            for h= -1/2:1/2
                for k = -1/2:1/2
                    for m = -1/2:1/2

                        % Coordenas del vecino.
                        [x_vec, y_vec, z_vec] = deal(x + h, y + k, z + m);

                        % Número del vecino.
                        N_vec = (N_celda_max + 1)^3 - (N_celda_max - x_vec) * (N_celda_max + 1)^2 - (N_celda_max - y_vec) * (N_celda_max + 1) - (N_celda_max - z_vec);

                        % Comprobamos que el vecino exista.
                        if (x_vec >= 0 && x_vec <= N_celda_max) && (y_vec >= 0 && y_vec <= N_celda_max) && (z_vec >= 0 && z_vec <= N_celda_max)
                            
                            % Definimos constantes del vecino.
                            A_N_vec = Info_atomo (N_vec, 3);
                            B_N_vec = Info_atomo (N_vec, 4);
                            Num_parejas_2 = Num_parejas_2 + 1; 

                            
                            % Calculamos las constantes promedio de cada pareja de vecinos.
                            A_mn = (A_N + A_N_vec)/2;
                            B_mn = (B_N + B_N_vec)/2;
                            
                            % Calculamos la energía de cada pareja de átomos.
                            U_mn =  - A_mn/r0_min + B_mn/(r0_min^2);
                            
                            % Almacenamos la información de cada pareja de vecinos en la matriz.
                            Info_pareja_2(Num_parejas_2, :) = [Info_atomo(N,1), Info_atomo(N_vec,1), Info_atomo(N,2), Info_atomo(N_vec,2), x, y, z, x_vec, y_vec, z_vec, U_mn];
                            Info_pareja_ampl = [Info_pareja_1; Info_pareja_2];
                            Matriz_info_pareja_ampl = [Matriz_info_pareja; num2cell(Info_pareja_2)];

                        else
                            continue;
                        end   
                    end
                end
            end
        end
    end
end


% ------------------------------------------------------------------------
% Buscar las posiciones de menor energía.
Orden_atomo_aleat = randperm(Num_atomo_T0);
Info_cambio = zeros (Num_atomo_T0, 9); % 8 columnas: Número del átomo, Tipo del átomo, phi_min, theta_min, rho_min, [x, y, z] finales, U_i.
Iter = 5;
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


% Calcular energía media por átomo.
Energia_media_atomo = Energia_T / Num_atomo_T0;

% Dar la energía media por átomo.
disp(['Energía total: ' num2str(Energia_T)]);
disp(['Número de átomos: ' num2str(Num_atomo_T0)]);
disp(['Energía media por átomo: ' num2str(Energia_media_atomo)]);

% Especificaciones de la representación.
title('BCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
