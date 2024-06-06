%% ESTRUCTURA FCC %%

% Proporciones de cada tipo de átomo.
Prop_1 = 0.4;
Prop_2 = 0.3;
Prop_3 = 0.3;

% Radios.
Radio_1 = 0.05;
Radio_2 = 0.08;
Radio_3 = 0.1;

% Constantes A y B para cada tipo de átomo.
A_1 = 0.1; 
B_1 = 0.2; 
A_2 = 0.1; 
B_2 = 0.2; 
A_3 = 0.1; 
B_3 = 0.2; 

% Número de celda unidad por lado.
N_celda = 1;

% Obtener el radio máximo.
Radio_max = max([Radio_1, Radio_2, Radio_3]);

% Distancia entre dos átomos consecutivos en las esquinas.
Parmtr_reticular = 2 * sqrt(2) * Radio_max;

% Distancia entre el centro y las esquinas de la celda.
Dist_centro_esquinas = sqrt(3) * sqrt(2) * Radio_max;

% Recuento energía.
Energia_T = 0;

% Representar la estructura.
figure;
hold on;
axis equal;

% Colocar átomos en las esquinas de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda*Parmtr_reticular
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
        end
    end
end

% Colocar el átomo C1 (1,1,0).
for x = 0:Parmtr_reticular:N_celda*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda*Parmtr_reticular
            % Posición del átomo C1.
            Pos_C1 = [x, y, z] + Parmtr_reticular * (1/2) * [1, 1, 0];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_C1 = rand();
            
            % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_C1 < Prop_1
                % Tipo 1.
                Color_C1 = 'r';
                Radio_C1 = Radio_1;
                A = A_1;
                B = B_1;
            elseif N_aleatorio_C1 < (Prop_1+Prop_2)
                % Tipo 2.
                Color_C1 = 'g';
                Radio_C1 = Radio_2;
                A = A_2;
                B = B_2;
            else
                % Tipo 3.
                Color_C1 = 'b';
                Radio_C1 = Radio_3;
                A = A_3;
                B = B_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_C1 + Pos_C1(1);
            y_sphere = y_sphere * Radio_C1 + Pos_C1(2);
            z_sphere = z_sphere * Radio_C1 + Pos_C1(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_C1, 'EdgeColor', 'none');
            
            % Distancia a los primeros vecinos.
            Dist_vecino = 2 * Radio_max;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;
        end
    end
end

% Colocar el átomo C2 (1,0,1).
for x = 0:Parmtr_reticular:N_celda*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda*Parmtr_reticular
            % Posición del átomo C2.
            Pos_C2 = [x, y, z] + Parmtr_reticular * (1/2) * [1, 0, 1];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_C2 = rand();
            
            % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_C2 < Prop_1
                % Tipo 1.
                Color_C2 = 'r';
                Radio_C2 = Radio_1;
                A = A_1;
                B = B_1;
            elseif N_aleatorio_C2 < (Prop_1+Prop_2)
                % Tipo 2.
                Color_C2 = 'g';
                Radio_C2 = Radio_2;
                A = A_2;
                B = B_2;
            else
                % Tipo 3.
                Color_C2 = 'b';
                Radio_C2 = Radio_3;
                A = A_3;
                B = B_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_C2 + Pos_C2(1);
            y_sphere = y_sphere * Radio_C2 + Pos_C2(2);
            z_sphere = z_sphere * Radio_C2 + Pos_C2(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_C2, 'EdgeColor', 'none');

            % Distancia a los primeros vecinos.
            Dist_vecino = 2 * Radio_max;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;
        end
    end
end

% Colocar el átomo C3 (0,1,1).
for x = 0:Parmtr_reticular:N_celda*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda*Parmtr_reticular
            % Posición del átomo C3.
            Pos_C3 = [x, y, z] + Parmtr_reticular * (1/2) * [0, 1, 1];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_C3 = rand();
            
            % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_C3 < Prop_1
                % Tipo 1.
                Color_C3 = 'r';
                Radio_C3 = Radio_1;
                A = A_1;
                B = B_1;
            elseif N_aleatorio_C3 < (Prop_1+Prop_2)
                % Tipo 2.
                Color_C3 = 'g';
                Radio_C3 = Radio_2;
                A = A_2;
                B = B_2;
            else
                % Tipo 3.
                Color_C3 = 'b';
                Radio_C3 = Radio_3;
                A = A_3;
                B = B_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_C3 + Pos_C3(1);
            y_sphere = y_sphere * Radio_C3 + Pos_C3(2);
            z_sphere = z_sphere * Radio_C3 + Pos_C3(3);
             
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_C3, 'EdgeColor', 'none');

            % Distancia a los primeros vecinos.
            Dist_vecino = 2 * Radio_max;
            
            % Calcular la energía de cada átomo.
            Energia_atomo = -(A / Dist_vecino) + (B / Dist_vecino^2);
            
            % Sumar la energía de todos los átomos.
            Energia_T = Energia_T + Energia_atomo;
        end
    end
end

% Calcular energía media por átomo.
Num_atomo_T = N_celda^3 * 4; % Número de átomos por celda unidad.
Energia_media_atomo = Energia_T / Num_atomo_T;

% Dar la energía.
disp(['Energía total: ' num2str(Energia_T)]);
disp(['Energía media por átomo: ' num2str(Energia_media_atomo)]);

% Especificaciones de la representación.
title('FCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
