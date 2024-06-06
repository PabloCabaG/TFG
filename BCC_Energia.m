%% ESTRUCTURA BCC %% 

% Proporciones de cada tipo de átomo.
Prop_1 = 0.4;
Prop_2 = 0.3;
Prop_3 = 0.3;

% Radios.
Radio_1 = 0.05; % a = 2.87 (Fe 1e-10 m)
Radio_2 = 0.08; % a = 3.30 (Nb 1e-10 m)
Radio_3 = 0.1;  % a = 5.585 (Rb 1e-10 m)

% Constantes A y B para cada tipo de átomo.
A_1 = 0.1; % A = 4.49 1e-27 Jm
B_1 = 0.2; % B = 5.58 1e-37 Jm^2
A_2 = 0.1; % A = 7.87 1e-27 Jm
B_2 = 0.2; % B = 11.24 1e-37 Jm^2
A_3 = 0.1; % A = 0.72 1e-27 Jm
B_3 = 0.2; % B = 1.74 1e-37 Jm^2

% Número de celda unidad por lado.
N_celda_extra_lado = 3; % N_celda_extra_lado = N_celda_lado - 1.

% Cálculo de energía antes de las microdeformaciones.
A_promedio = Prop_1 * A_1 + Prop_2 * A_2 + Prop_3 * A_3;
B_promedio = Prop_1 * B_1 + Prop_2 * B_2 + Prop_3 * B_3;
r0 = 2:0.01:7;
U0 = - A_promedio./r0 + B_promedio./(r0.^2);
hold on
plot (r0, U0, 'LineWidth', 1.5);
r0_min = 2*B_promedio/A_promedio;
U0_min = - A_promedio/r0_min + B_promedio/(r0_min^2);
plot (r0_min,U0_min,'.', 'MarkerSize', 25)
xlabel ('r_0')
ylabel ('U_0')
hold off

% Obtener el radio máximo.
Radio_max = max([Radio_1, Radio_2, Radio_3]);
Radio_min = min([Radio_1, Radio_2, Radio_3]);

% Distancia entre dos átomos consecutivos en las esquinas.
Parmtr_reticular = 4 / sqrt(3) * Radio_max;

% Distancia entre el centro y las esquinas de la celda.
Dist_centro_esquinas = 2 * Radio_max;

% Recuento número de átomos.
Num_atomo_T = 0;

% Recuento energía.
Energia_T = 0;

% Representar la estructura.
figure;
hold on;
axis equal;

% Colocar átomos en las esquinas de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
            % Posición de los átomos
            Pos_esquina = [x, y, z];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_esq = rand();
            
            % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_esq < Prop_1
                % Tipo 1.
                Color_esq = 'r';
                Radio_esq = Radio_1;
                A = A_1;
                B = B_1;
            elseif N_aleatorio_esq < (Prop_1+Prop_2)
                % Tipo 2.
                Color_esq = 'g';
                Radio_esq = Radio_2;
                A = A_2;
                B = B_2;
            else
                % Tipo 3.
                Color_esq = 'b';
                Radio_esq = Radio_3;
                A = A_3;
                B = B_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_esq + Pos_esquina(1);
            y_sphere = y_sphere * Radio_esq + Pos_esquina(2);
            z_sphere = z_sphere * Radio_esq + Pos_esquina(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_esq, 'EdgeColor', 'k');
            
            % Distancia a los primeros vecinos.
            Dist_vecino = 2 * Radio_max;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T + 1;

        end
    end
end

% Colocar átomos en el centro de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda_extra_lado*Parmtr_reticular
            % Posición del átomo en el centro de la celda.
            Pos_centro = [x, y, z] + Parmtr_reticular * (1/2) * [1, 1, 1];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_centro = rand();
            
            % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_centro < Prop_1
                % Tipo 1.
                Color_centro = 'r';
                Radio_centro = Radio_1;
                A = A_1;
                B = B_1;
            elseif N_aleatorio_centro < (Prop_1+Prop_2)
                % Tipo 2.
                Color_centro = 'g';
                Radio_centro = Radio_2;
                A = A_2;
                B = B_2;
            else
                % Tipo 3.
                Color_centro = 'b';
                Radio_centro = Radio_3;
                A = A_3;
                B = B_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_centro + Pos_centro(1);
            y_sphere = y_sphere * Radio_centro + Pos_centro(2);
            z_sphere = z_sphere * Radio_centro + Pos_centro(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_centro, 'EdgeColor', 'none');
            
            % Distancia a los primeros vecinos.
            Dist_vecino = 2 * Radio_max;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;

            % Número de átomos total.
            Num_atomo_T = Num_atomo_T + 1;

        end
    end
end

% Calcular energía media por átomo.
% Num_atomo_T = N_celda^3 * 2; % Número de átomos por celda unidad.
Energia_media_atomo = Energia_T / Num_atomo_T;

% Dar la energía media por átomo.
disp(['Energía total: ' num2str(Energia_T)]);
disp(['Número de átomos: ' num2str(Num_atomo_T)]);
disp(['Energía media por átomo: ' num2str(Energia_media_atomo)]);

% Especificaciones de la representación.
title('BCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
