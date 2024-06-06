%% ESTRUCTURA BCC %%

% Proporciones de cada tipo de átomo.
Prop_1 = 0.4;
Prop_2 = 0.3;
Prop_3 = 0.3;

% Radios.
Radio_1 = 0.05;
Radio_2 = 0.08;
Radio_3 = 0.1;

% Número de celda unidad por lado.
N_celda = 3;

% Obtener el radio máximo.
Radio_max = max([Radio_1, Radio_2, Radio_3]);

% Distancia entre dos átomos consecutivos en las esquinas.
Parmtr_reticular = 4 / sqrt(3) * Radio_max;

% Distancia entre el centro y las esquinas de la celda.
Dist_centro_esquinas = 2 * Radio_max;

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
            elseif N_aleatorio_esq < (Prop_1+Prop_2)
                % Tipo 2.
                Color_esq = 'g';
                Radio_esq = Radio_2;
            else
                % Tipo 3.
                Color_esq = 'b';
                Radio_esq = Radio_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_esq + Pos_esquina(1);
            y_sphere = y_sphere * Radio_esq + Pos_esquina(2);
            z_sphere = z_sphere * Radio_esq + Pos_esquina(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_esq, 'EdgeColor', 'none');
        end
    end
end

% Colocar átomos en el centro de cada celda unidad.
for x = 0:Parmtr_reticular:N_celda*Parmtr_reticular
    for y = 0:Parmtr_reticular:N_celda*Parmtr_reticular
        for z = 0:Parmtr_reticular:N_celda*Parmtr_reticular
            % Posición del átomo en el centro de la celda.
            Pos_centro = [x, y, z] + Parmtr_reticular * (1/2) * [1, 1, 1];
            
            % Generar un número aleatorio para determinar el tipo de átomo.
            N_aleatorio_centro = rand();
            
           % Determinar el tipo de átomo según las proporciones dadas.
            if N_aleatorio_centro < Prop_1
                % Tipo 1.
                Color_centro = 'r';
                Radio_centro = Radio_1;
            elseif N_aleatorio_centro < (Prop_1+Prop_2)
                % Tipo 2.
                Color_centro = 'g';
                Radio_centro = Radio_2;
            else
                % Tipo 3.
                Color_centro = 'b';
                Radio_centro = Radio_3;
            end
            
            % Representar el átomo con el color y radio correspondiente.
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * Radio_centro + Pos_centro(1);
            y_sphere = y_sphere * Radio_centro + Pos_centro(2);
            z_sphere = z_sphere * Radio_centro + Pos_centro(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', Color_centro, 'EdgeColor', 'none');
        end
    end
end

% Especificaciones de la representación.
title('BCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
