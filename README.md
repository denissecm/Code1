% Simulación para calcular y visualizar el campo magnético con la Ley Biot Savart 

%Limpia el área de trabajo 
clear 
clc  
clf 

% Parámetros del selenoide y espiras 

nl = 5; % Número de espiras 
N = 20; % Número de puntos por espira 
R = 1.5; % Radio por espira 
sz = 1; % Separación entre espiras 
I = 300; % Corriente 
mo = 4*pi*1e-7; % Permeabilidad magnética 
km = mo * I / (4*pi); % Constante de Biot-Savart
rw = 0.2; % Grosor del alambre 

% Ángulos para la espira
dtheta = 2*pi / N; % Incremento angular 
ang = 0:dtheta:(2*pi - dtheta); % Vector de ángulos 

% Posiciones y componentes de los vectores 
s = 1; % Índice de inicio 
for i = 1:nl % Bucle en la espira
% Coordenadas 
    Px(s:s+N-1) = R * cos(ang);
    Py(s:s+N-1) = R * sin(ang);
    Pz(s:s+N-1) = -nl/2*sz + (i-1)*sz;

    % Componentes X e Y del vector dl
    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;
    dy(s:s+N-1) =  Px(s:s+N-1) * dtheta;

    s = s + N; % Actualiza el índice de inicio 
end
dz = zeros(1, N*nl); % Componente Z de dl es cero, porque se calculan individual las espiras 

% Visualización de la corriente en las espiras
figure(1) 
quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 
% '-r' define las flechas en color rojo, 'LineWidth', 2 define el grosor de la línea.
view(-34, 33) % ángulo de visión de la figura 3D
title('Corriente en espiras')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal % proporcion de ejes 

%Calculo del campo magnetico 
% Malla del campo 
ds = 0.1; % Tamaño de malla 
x = -5:ds:5; % Rango X
y = x; % Rango Y 
z = x; % Rango Z 
Lx = length(x); % Número de puntos en el eje X
Ly = length(y); % Número de puntos en el eje Y
Lz = length(z); % Número de puntos en el eje Z

% Inicialización de los componentes del campo magnético 
dBx = zeros(Lx, Ly, Lz);
dBy = zeros(Lx, Ly, Lz);
dBz = zeros(Lx, Ly, Lz);

% Cálculo del campo en cada punto del espacio 
tic % cronometro para el tiempo de calculo 
for i = 1:Lx % Bucle sobre las coordenadas X
    for j = 1:Ly % Bucle sobre las coordenadas Y 
        for k = 1:Lz % Bucle sobre las coordenadas Z 

            for l = 1:nl*N % Bucle sobre cada elemento de corriente de las espiras 
                % rx, ry, rz: Componentes del vector 'r' 
                rx = x(i) - Px(l);
                ry = y(j) - Py(l);
                rz = z(k) - Pz(l);

                % r: Magnitud del vector 'r'. 
                r = sqrt(rx^2 + ry^2 + rz^2 + rw^2);
                r3 = r^3; 
                
                dBx(i,j,k) = dBx(i,j,k) + km * (dy(l) * rz - dz(l) * ry) / r3;
                dBy(i,j,k) = dBy(i,j,k) + km * (dz(l) * rx - dx(l) * rz) / r3;
                dBz(i,j,k) = dBz(i,j,k) + km * (dx(l) * ry - dy(l) * rx) / r3;
            end
        end
    end
end
toc % Detiene el cronómetro y muestra el tiempo 

% Magnitud del campo magnético total 
Bmag = sqrt(dBx.^2 + dBy.^2 + dBz.^2);

% Corte en el plano XZ para visualización
centery = round(Ly / 2); %Encuentra el índice del plano Y más cercano a Y=0.
Bx_xz = squeeze(dBx(:, centery, :)); % Extrae el plano XZ de la componente Bx
Bz_xz = squeeze(dBz(:, centery, :)); % Extrae el plano XZ de la componente Bz
Bxz = squeeze(Bmag(:, centery, :)); % Extrae el plano XZ de la magnitud del campo

% Visualización del campo en el plano XZ
figure(2) 
hold on
pcolor(x, z, (Bxz').^(1/3)); shading interp; colormap jet; colorbar
h1 = streamslice(x, z, Bx_xz', Bz_xz', 3); % Dibuja líneas de campo magnético en el plano XZ.
set(h1, 'Color', [0.8 1 0.9]);
xlabel('x'); ylabel('z'); 
title('Campo magnético generado por un solenoide') 
% (Bxz') se usa para transponer la matriz
% shading interp: Suaviza los colores entre los puntos de la malla.
% colormap jet: Establece el mapa de colores.
% colorbar: Añade una barra de color para interpretar los valores.
% streamslice: Dibuja líneas de campo
% h1: Handle a los objetos de línea de campo
