clear;

%% ESTRUCTURA BCC %%

% Número de celda unidad por lado.
N_celda_extra_lado = 3; % N_celda_extra_lado = N_celda_lado - 1.

% Proporciones de cada tipo de átomo.
Prop_1 = 52/128;
Prop_2 = 38/128;
Prop_3 = 38/128;

% Recuento número de átomos.
Num_atomo_T0 = 0;
Num_atomo_T = (N_celda_extra_lado + 1)^3 * 2; % Recuento inverso.
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
Parmtr_reticular = 2 * r0_min / sqrt(3);

% Distancia entre el centro y las esquinas de la celda.
Dist_centro_esquinas = r0_min;

% Recuento energía.
Energia_T = 0;

% Matriz donde se almacena la información de los átomos.
Info_atomo = zeros((N_celda_extra_lado + 1)^3 * 2, 5); % 5 columnas: Número, Tipo, X, Y, Z

% Representar la estructura.
figure;
hold on;
axis equal;

% Colocar átomos en las esquinas de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
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
            Nombre_columnas = {'Número del átomo', 'Tipo del átomo', 'X', 'Y', 'Z'};
            Info_atomo(Num_atomo_T0 + 1, :) = [Num_atomo_T0 + 1, Tipo_atom, Pos_esquina/Parmtr_reticular];
            Matriz_info = [Nombre_columnas; num2cell(Info_atomo)];       

            % Acumular número de átomos.
            Num_atomo_T0 = Num_atomo_T0 + 1;

            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_esq + Pos_esquina(1);
            y_sphere = y_sphere * Radio_esq + Pos_esquina(2);
            z_sphere = z_sphere * Radio_esq + Pos_esquina(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_esq, 'EdgeColor', 'k');
            
            % Distancia a los primeros vecinos.
            Dist_vecino = r0_min;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T - 1;

        end
    end
end

% Colocar átomos en el centro de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
            % Posición del átomo en el centro de la celda.
            Pos_centro = [x, y, z] + Parmtr_reticular * (1/2) * [1, 1, 1];
            
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
            Info_atomo(Num_atomo_T0 + 1, :) = [Num_atomo_T0 + 1, Tipo_atom, Pos_centro/Parmtr_reticular];
            Matriz_info = [Nombre_columnas; num2cell(Info_atomo)];

            % Acumular número de átomos.
            Num_atomo_T0 = Num_atomo_T0 + 1;


            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_centro + Pos_centro(1);
            y_sphere = y_sphere * Radio_centro + Pos_centro(2);
            z_sphere = z_sphere * Radio_centro + Pos_centro(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_centro, 'EdgeColor', 'none');
            
            % Distancia a los primeros vecinos.
            Dist_vecino = r0_min;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T - 1;

        end
    end
end

% Calcular energía media por átomo.
Energia_media_atomo = Energia_T / Num_atomo_T0;

% Dar la energía media por átomo.
disp(['Energía total: ' num2str(Energia_T)]);
disp(['Número de átomos: ' num2str(Num_atomo_T0)]);
disp(['Energía media por átomo: ' num2str(Energia_media_atomo)]);

% Especificaciones de la representación.
set(gca, 'FontSize', 15);
title('BCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
