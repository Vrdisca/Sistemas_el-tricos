%Fluxo linearizado+ perdas
clc
clear all
%Dados da Geração
%[Barra Geração(pu) Carga (pu)]
BusData = [1 1.9 0.0;
 2 0.0 0.6;
 3 0.0 0.7;
 4 0.0 0.8;
 5 0.85 0.65];
LinData = [1 2 0.0108 0.0649;
 1 4 0.0235 0.0941;
 2 5 0.0118 0.0471;
 3 5 0.0147 0.0588;
 4 5 0.0118 0.0529];

%Número de barras
[NBus, ncol] = size(BusData);
%Número de linhas
[NLin, ncol] = size(LinData)

% Construção da matriz de susceptância
BBus(1:NBus,1:NBus) = 0;
g(1:NLin) = 0;
for il = 1:NLin
 k = LinData(il,1);
 m = LinData(il,2);
 r = LinData(il,3);
 x = LinData(il,4);
 g(il) = r/(r^2 + x^2);
 b = -1/x;
 BBus(k,k) = BBus(k,k) + b;
 BBus(k,m) = -b;
 BBus(m,k) = -b;
 BBus(m,m) = BBus(m,m) + b;
 end
% Eliminando a barra de referência
B = BBus(1:NBus-1,1:NBus-1);
% Vetor de potências injetadas
P(1:NBus-1,1) = 0;
for ib = 1:NBus-1
 P(ib) = BusData(ib,2) - BusData(ib,3);
end
% Calculando os ângulos
Theta(1:NBus,1) = 0;
Theta(1:NBus-1,1) = -inv(B)*P;
PLOld(1:NBus) = 0;
erro = 1e10;
tol = 1e-3;
ThetaOld= Theta
while erro >= tol
 %Calcular as perdas
 PL(1:NBus) = 0;
 for il = 1:NLin
 k = LinData(il,1);
 m = LinData(il,2);
 Loss = g(il)*(Theta(k)-Theta(m))^2;
 PL(k) = PL(k) + Loss/2;
PL(m) = PL(m) + Loss/2;
 end
 % Calculando o vetor de injeções considerando as perdas
 for ib = 1:NBus-1
 P(ib) = BusData(ib,2)- (BusData(ib,3) + PL(ib));
 end
 % Calculando os desvios das perdas
 DeltaLoss = PL(2:NBus) - PLOld(2:NBus);
 erro = max(abs(DeltaLoss));
 % Calculando os novos ângulos
 Theta(1:NBus-1,1) = -inv(B)*P;
 PLOld = PL;
end
% Fluxos nas Linhas
Flow(1:NLin,1) = 0;
for il = 1:NLin
 k = LinData(il,1);
 m = LinData(il,2);
b = -1/LinData(il,4);
 Flow(il) = -b*(Theta(k)-Theta(m));
end
sum(PL)% Soma das perdas