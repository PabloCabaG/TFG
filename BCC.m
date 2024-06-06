% Parámetros
radio_tipo1 = 0.05;
radio_tipo2 = 0.08;
radio_tipo3 = 0.1;

% Obtener el radio máximo
radio_maximo = max([radio_tipo1, radio_tipo2, radio_tipo3]);

% Distancia entre dos átomos consecutivos en las esquinas
distancia_entre_esquinas = 4 / sqrt(3) * radio_maximo;

% Distancia entre el centro y las esquinas de la celda
distancia_centro_esquinas = 2 * radio_maximo;

% Crear una figura
figure;
hold on;
axis equal;

% Bucle para colocar átomos en cada posición entera de la red
for x = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
    for y = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
        for z = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
            % Posición de los átomos
            pos = [x, y, z];
            
            % Generar un número aleatorio para determinar el tipo de átomo
            tipo_aleatorio = rand();
            
            % Determinar el tipo de átomo según las proporciones especificadas
            if tipo_aleatorio < 0.4
                % Tipo 1
                color = 'r';
                radio = radio_tipo1;
            elseif tipo_aleatorio < 0.7
                % Tipo 2
                color = 'g';
                radio = radio_tipo2;
            else
                % Tipo 3
                color = 'b';
                radio = radio_tipo3;
            end
            
            % Representar átomo como una esfera con el color y radio correspondiente
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * radio + pos(1);
            y_sphere = y_sphere * radio + pos(2);
            z_sphere = z_sphere * radio + pos(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', color, 'EdgeColor', 'none');
        end
    end
end

% Bucle para colocar átomos en el centro de cada cubo
for x = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
    for y = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
        for z = 0:distancia_entre_esquinas:3*distancia_entre_esquinas
            % Posición del átomo en el centro del cubo
            pos_centro_cubo = [x, y, z] + distancia_entre_esquinas * (1/2) * [1, 1, 1];
            
            % Generar un número aleatorio para determinar el tipo de átomo en el centro
            tipo_aleatorio_centro = rand();
            
            % Determinar el tipo de átomo en el centro según las proporciones especificadas
            if tipo_aleatorio_centro < 0.4
                % Tipo 1
                color_centro_cubo = 'r';
                radio_centro_cubo = radio_tipo1;
            elseif tipo_aleatorio_centro < 0.7
                % Tipo 2
                color_centro_cubo = 'g';
                radio_centro_cubo = radio_tipo2;
            else
                % Tipo 3
                color_centro_cubo = 'b';
                radio_centro_cubo = radio_tipo3;
            end
            
            % Representar átomo en el centro como una esfera con el color y radio correspondiente
            [x_sphere, y_sphere, z_sphere] = sphere;
            x_sphere = x_sphere * radio_centro_cubo + pos_centro_cubo(1);
            y_sphere = y_sphere * radio_centro_cubo + pos_centro_cubo(2);
            z_sphere = z_sphere * radio_centro_cubo + pos_centro_cubo(3);
            
            surf(x_sphere, y_sphere, z_sphere, 'FaceColor', color_centro_cubo, 'EdgeColor', 'none');
        end
    end
end

% Configuración de la figura
title('BCC');
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
