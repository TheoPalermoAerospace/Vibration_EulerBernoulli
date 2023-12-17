clear
clc

% Implementação do MEF para obter as frequências naturais

% Definindo os parâmetros
I = 0.0008; % Momento de inércia
rho = 7800; % Densidade
L = 3; % Comprimento da viga
E = 200000000000; % Módulo de elasticidade
A = 0.02; % Área de seção transversal
n = 20; % Número de elementos
p_0 = 250000; % Amplitude da força
t = 0; % Tempo atual

% Matriz de rigidez do elemento
k = [12 6*L -12 6*L; 6*L 4*L^2 -6*L 2*L^2; -12 -6*L 12 -6*L; 6*L 2*L^2 -6*L 4*L^2]*E*I/L^3;
% Matriz de massa do elemento
m = [156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2; 54 13*L 156 -22*L; -13*L -3*L^2 -22*L 4*L^2]*rho*A*L/420;
% Inicialização das matrizes globais
M = zeros(2*n+2);
K = zeros(2*n+2);

% Montagem das matrizes globais
for i = 1:n
  M(2*i-1:2*i+2, 2*i-1:2*i+2) = M(2*i-1:2*i+2, 2*i-1:2*i+2) + m;
  K(2*i-1:2*i+2, 2*i-1:2*i+2) = K(2*i-1:2*i+2, 2*i-1:2*i+2) + k;
endfor

K_modificado = K(3:end, 3:end);
M_modificado = M(3:end, 3:end);

[V,W] = eig(K_modificado, M_modificado);
omega = sqrt(W);

% Resposta do Sistema
u0 = zeros(40);
u0_dot = zeros(40);
u0_ddot = inv(M_modificado)*(K_modificado*(-u0));
h = 0.4;
a1 = ((4/(h^2))*M_modificado);
a2 = ((4/h)*M_modificado);
a3 = M_modificado;
k_til = a1+K_modificado;

vetorDeslocamentos = zeros(101, n); % Alterado para armazenar apenas os deslocamentos
uAntigo = u0;
udotAntigo = u0_dot;
uddotAntigo = u0_ddot;

i = 1;
for t = 0:1:190
  % Vetor de força variante no tempo
  F = zeros(2*n, 1);
  if t <= 6.35234
    F(end) = -p_0 * sin(pi*t/6.35234); % Força aplicada no último nó
  endif
  p_til = F + a1*uAntigo + a2*udotAntigo + a3*uddotAntigo;
  uNovo = inv(k_til)*(p_til);
  udotNovo = (-2/h)*uAntigo - udotAntigo + (2/h)*uNovo;
  uddotNovo = (-4/(h^2))*uAntigo - (4/h)*udotAntigo - uddotAntigo + (4/(h^2))*uNovo;
  uAntigo = uNovo;
  udotAntigo = udotNovo;
  uddotAntigo = uddotNovo;
  vetorDeslocamentos(i, :) = transpose(uAntigo(1:2:end))(1:size(vetorDeslocamentos, 2)); % Alterado para considerar apenas os deslocamentos
  i = i + 1;
endfor

% Selecionando os três primeiros modos de vibração
modos = V(:,1:3);
% Normalizando os modos de vibração
modos = modos ./ max(abs(modos));
% Selecionando apenas os elementos ímpares (deslocamentos)
modos = modos(1:2:end,:);
% Gerando um conjunto de pontos mais densamente espaçados
x_fino = linspace(1, size(modos, 1), size(modos, 1));
% Plotando os modos de vibração
figure;
for i = 1:3
    subplot(3,1,i);
    scatter(x_fino, modos(:,i), 'filled');
    title(['Modo de vibração ', num2str(i)]);
end

tempo = 0:1:190;
figure;
hold on;
% Plota cada coluna de vetorDeslocamentos em relação ao tempo
for j = 1:size(vetorDeslocamentos, 2)
  plot(tempo, vetorDeslocamentos(:, j));
endfor

xlabel('Tempo');
ylabel('Deslocamentos', 'FontSize', 14);
title('Deslocamentos vs Tempo', 'FontSize', 14);
hold off;

x = linspace(0, L, n);
for i = 1:size(vetorDeslocamentos, 1)
    plot(x, vetorDeslocamentos(i, :));
    hold on;
endfor

x_nos = linspace(0, L, n);
for i = 1:size(vetorDeslocamentos, 1)
    plot(x_nos, vetorDeslocamentos(i, 1:n), 'o');
    hold on;
endfor

title('Deflexão da viga em função do comprimento');
xlabel('Comprimento da viga (m)', 'FontSize', 16);
ylabel('Deflexão (m)', 'FontSize', 16);
hold off;

x = linspace(0, L, n);
figure;
for i = 1:size(vetorDeslocamentos, 1)
    plot(x, vetorDeslocamentos(i, :), 'b-');
    hold on;
    plot(x, vetorDeslocamentos(i, :), 'bo');
    grid on;
    line([min(x) max(x)], [0 0], 'Color', 'red');
    hold off;
    title('Deflexão da viga em função do comprimento');
    xlabel('Comprimento da viga (m)', 'FontSize', 14);
    ylabel('Deflexão (m)', 'FontSize', 14);
    ylim([-max(abs(vetorDeslocamentos(:))) max(abs(vetorDeslocamentos(:)))]);
    pause(0.1); % Pausa por 0.1 segundos antes de plotar o próximo gráfico
endfor
