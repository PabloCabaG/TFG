% Parámetros iniciales
num_atomos = 144; % Número total de átomos
radio_atomos = [0.1, 0.15, 0.2]; % Radios de los tres tipos de átomos
proporciones = [0.4, 0.3, 0.3]; % Proporción de átomos de cada tipo (ejemplo: 40% de Tipo 1, 30% de Tipo 2, 30% de Tipo 3)
num_tipos_atomos = length(radio_atomos); % Número de tipos de átomos
lado_cuadrado = sqrt(num_atomos); % Lado de la red cuadrada

% Crea una matriz para almacenar las coordenadas de los átomos
coordenadas = zeros(num_atomos, 2);
tipos_atomos = [];

% Calcula las coordenadas y tipos de átomos
for tipo = 1:num_tipos_atomos
    num_atomos_tipo = round(num_atomos * proporciones(tipo));
    tipos_atomos = [tipos_atomos; repmat(tipo, num_atomos_tipo, 1)];
end

% Mezcla los tipos de átomos para asignación aleatoria
tipos_atomos = tipos_atomos(randperm(num_atomos));

% Calcula las coordenadas de los átomos sin nodos libres
for i = 1:num_atomos
    tipo = tipos_atomos(i);
    [row, col] = ind2sub([lado_cuadrado, lado_cuadrado], i);
    x = (col - 1) * (2 * max(radio_atomos)); % Multiplica por 2 para asegurar la misma arista
    y = (row - 1) * (2 * max(radio_atomos)); % Multiplica por 2 para asegurar la misma arista
    coordenadas(i, :) = [x, y];
end

% Dibuja la red cuadrada simple con átomos de diferentes colores
figure;

% Define una matriz de colores para cada tipo de átomo
colores = hsv(num_tipos_atomos);

for tipo = 1:num_tipos_atomos
    tipo_indices = find(tipos_atomos == tipo);
    scatter(coordenadas(tipo_indices, 1), coordenadas(tipo_indices, 2), radio_atomos(tipo)^2 * 1000, colores(tipo, :), 'filled');
    hold on;
end

axis equal;
title('Red Cuadrada Simple con Átomos de Tres Tipos, Proporciones Definidas y Arista Igual para Todos los Átomos');
xlabel('Coordenada X');
ylabel('Coordenada Y');

% Crear una leyenda para los tipos de átomos
leyenda = cell(num_tipos_atomos, 1);
for tipo = 1:num_tipos_atomos
    leyenda{tipo} = sprintf('Tipo %d', tipo);
end
legend(leyenda);

