function Phase_Space_Main(Axes, PlotParameters)

SystemParameters = PlotParameters.SystemParameters;
FigureParameters = PlotParameters.FigureParameters;
NumericalParameters = PlotParameters.NumericalParameters;
ListOfComponents = PlotParameters.ListOfComponents;

% --- Задание параметров графика --- %

InitializePlot(Axes, SystemParameters, FigureParameters)

% --- Нанесение областей D --- %

AreaOfD(Axes, SystemParameters);

% --- Нанесение основных траекторий --- %
if ListOfComponents.MainTrajectories == 1
    MainTrajectories(Axes, SystemParameters, FigureParameters, NumericalParameters);
end

% --- Нанесение изоклин --- %
if ListOfComponents.Isocline == 1
    Isocline(Axes, SystemParameters)
end

%Нанесение векторного поля
if ListOfComponents.VectorField == 1
    VectorField(Axes, SystemParameters);
end

% --- Нанесение сепаратрис --- %

if ListOfComponents.Separatrix == 1
   Separatrix(Axes, SystemParameters, FigureParameters, NumericalParameters);
end

% --- Нанесение вспомогательной траектории --- %
if ListOfComponents.SupportTrajectory == 1
   SupportTrajectory(Axes, SystemParameters, FigureParameters, NumericalParameters); 
end
% --- Нанесение состояний равновесия --- %
if ListOfComponents.Equilibria == 1
   Equilibria(Axes, SystemParameters); 
end

end

function InitializePlot(ax, SP, FP)

%figure;
cla(ax);
hold(ax, 'on');
grid(ax, 'on');
%LineWidth = 3;

ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';



ax.XTick = [0, pi / 2, pi, 3*pi/2, 2 * pi];
ax.YTick = [0, pi / 2, pi, 3*pi/2, 2 * pi];
ax.XTickLabel = {'0', '$\pi / 2$', '$\pi$', '$3\pi / 2$', '$2\pi$'};
ax.YTickLabel = {'0', '$\pi / 2$', '$\pi$', '$3\pi / 2$', '$2\pi$'};
%ax.XTickLabel = {'0', '$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'};
%ax.YTickLabel = {'0', '$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'};

FontSize = 24;
ax.FontSize = FontSize;

ax.XLim = [0, 2*pi];
ax.YLim = [0, 2*pi];

AxisLabelsFontSize = 24;

ax.Box = 'on';

ax.XLabel.String = '$\phi_1$';
ax.XLabel.FontSize = AxisLabelsFontSize;

ax.YLabel.String = '$\phi_2$';
ax.YLabel.FontSize = AxisLabelsFontSize;
%ax.YLabel.Rotation = 0;

ax.XTickLabelRotation = 0;

end

function AreaOfD(ax, SP)

Sigma = SP.Sigma;
% 
% Center = SP.Gamma(1) - SP.d / 2;
% if Center > 1
%     Center = pi / 2;
% elseif Center < -1
%         Center = 3 * pi / 2;
% else
%     Center = asin(Center);
% end

Center = pi/2;

D = [Center - Sigma, Center + Sigma];

fill(ax, [D(1), D(2), D(2), D(1)] + 2*pi, [0, 0, 2 * pi, 2 * pi], [159/255, 1, 160/255], 'EdgeColor', 'none');
fill(ax, [0, 0, 2 * pi, 2 * pi], [D(1), D(2), D(2), D(1)] + 2*pi, [159/255, 1, 160/255], 'EdgeColor', 'none');
fill(ax, [D(1), D(2), D(2), D(1)] + 2*pi, [D(1), D(1), D(2), D(2)] + 2*pi, [97/255, 1, 98/255], 'EdgeColor', 'none');

fill(ax, [D(1), D(2), D(2), D(1)], [0, 0, 2 * pi, 2 * pi], [159/255, 1, 160/255], 'EdgeColor', 'none');
fill(ax, [0, 0, 2 * pi, 2 * pi], [D(1), D(2), D(2), D(1)], [159/255, 1, 160/255], 'EdgeColor', 'none');
fill(ax, [D(1), D(2), D(2), D(1)], [D(1), D(1), D(2), D(2)], [97/255, 1, 98/255], 'EdgeColor', 'none');

plot(ax, [Center, Center] + 2*pi, [0, 2 * pi], 'Color', [0.1 1 0.1]);
plot(ax, [0, 2 * pi], [Center, Center] + 2*pi, 'Color', [0.1 1 0.1]);

plot(ax, [Center, Center], [0, 2 * pi], 'Color', [0.1 1 0.1]);
plot(ax, [0, 2 * pi], [Center, Center], 'Color', [0.1 1 0.1]);
end

function Isocline(ax, SP)

Gamma = SP.Gamma;
Sigma = SP.Sigma;
d = SP.d;

D = [pi/2 - Sigma, pi/2 + Sigma];

NOfPhi = 10000;
Phi = [0 : 10/NOfPhi : 10];
Phi1 = zeros(length(Phi), 1);
Phi2 = zeros(length(Phi), 1);

k = -500;
Alpha = D(1);
Delta = 2 * Sigma;


for i = 1 : NOfPhi + 1
    Arg = Gamma(1) - d * F2func(Phi(i), k, Delta, Alpha);
    if abs(Arg) <= 1
    Phi1(i) =  asin(Arg);
    Phi2(i) = pi - asin(Arg);
    else
        Phi1(i) = NaN;
        Phi2(i) = NaN;
    end  
end

       Phi1 = mod(Phi1, 2 * pi);
       Phi2 = mod(Phi2, 2 * pi);
        for i = 1:length(Phi1) - 1
           if (abs(Phi1(i) - Phi1(i + 1)) > pi)
               Phi1(i) = NaN;
           end
           if (abs(Phi2(i) - Phi2(i + 1)) > pi)
               Phi2(i) = NaN;
           end
        end
        
plot(ax, Phi1, Phi, 'm', 'LineWidth', 2);
plot(ax, Phi2, Phi, 'm', 'LineWidth', 2);


for i = 1 : NOfPhi + 1
    Arg = Gamma(2) - d * F2func(Phi(i), k, Delta, Alpha);
    if abs(Arg) <= 1
    Phi1(i) =  asin(Arg);
    Phi2(i) = pi - asin(Arg);
    else
        Phi1(i) = NaN;
        Phi2(i) = NaN;
    end  
end

       Phi1 = mod(Phi1, 2 * pi);
       Phi2 = mod(Phi2, 2 * pi);
        for i = 1:length(Phi1) - 1
           if (abs(Phi1(i) - Phi1(i + 1)) > pi)
               Phi1(i) = NaN;
           end
           if (abs(Phi2(i) - Phi2(i + 1)) > pi)
               Phi2(i) = NaN;
           end
        end
        
plot(ax, Phi, Phi1, 'm', 'LineWidth', 2);
plot(ax, Phi, Phi2, 'm', 'LineWidth', 2);

end

function MainTrajectories(ax, SP, FP, NP)

Step = 30;
[X,Y] = meshgrid(0 : pi/Step : 2*pi, 0 : pi/Step : 2*pi);

Sigma = SP.Sigma;

% Center = SP.Gamma(1) - SP.d / 2;
% if Center > 1
%     Center = pi / 2;
% elseif Center < -1
%         Center = 3 * pi / 2;
% else
%     Center = asin(Center);
% end

Center = pi/2;

D = [Center - Sigma, Center + Sigma];
k = -500;
Alpha = D(1);
Delta = 2 * Sigma;

U = SP.Gamma(1) - sin(X) - SP.d .* 1 ./ (1 + exp(k .* (cos(Delta / 2) - cos(Y - Alpha - Delta ./ 2))));
V = SP.Gamma(2) - sin(Y) - SP.d .* 1 ./ (1 + exp(k .* (cos(Delta / 2) - cos(X - Alpha - Delta ./ 2))));

density = 0.5;
h = streamslice(ax, X, Y, U, V, density, 'method', 'cubic');
for i = 1 : length(h)
    h(i).Color = [1 0.5 0.5];
    h(i).LineWidth = 1;
end

end

function SupportTrajectory(ax, SP, FP, NP)

LineWidth = 2;
    Tspan = [0 NP.Tspan];
    X_0 = NP.InitialCondition;
    %opts = odeset('RelTol',1e-6,'AbsTol',1e-7, 'MaxStep', 1e-1);
   [T, X] = ode23(@(T, X) AdlerSystem(X, SP), Tspan, X_0, NP.Options);

   X = mod(X, 2 * pi);
    for i = 1:length(X) - 1
       if (abs(X(i, 1) - X(i + 1, 1)) > pi)
           X(i, 1) = NaN;
       end
       if (abs(X(i, 2) - X(i + 1, 2)) > pi)
           X(i, 2) = NaN; 
       end
    end

   plot(ax, X(:, 1), X(:, 2), 'Color', 'blue', 'LineWidth', LineWidth);
   plot(ax, X_0(1), X_0(2), '.', 'MarkerSize', 30, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', 'blue', 'LineWidth', 2);

end

function Separatrix(ax, SP, FP, NP)

k = -500;
Sigma = SP.Sigma;
d = SP.d;
gamma = SP.Gamma;

%Кол-во седловых СР

syms x y;
syms a b;

Lambda = zeros(1, 2);

opts = odeset('RelTol',1e-6,'AbsTol',1e-7, 'MaxStep', 1e-1);
LineWidth = 5;

%Приближенные СР
X_0 = mod(real([pi/2 - Sigma, pi/2 - Sigma; ...
       pi/2 + Sigma, pi/2 - Sigma; ...
       pi/2 - Sigma, pi/2 + Sigma; ...
       pi/2 + Sigma, pi/2 + Sigma;
       asin(gamma(1) - d), asin(gamma(1) - d);
       pi - asin(gamma(1) - d), asin(gamma(1) - d);
       asin(gamma(1) - d), pi - asin(gamma(1) - d);
       pi - asin(gamma(1) - d), pi - asin(gamma(1) - d);]), 2*pi);

%Вспомогательная функция - производная от ф-ии связи
dI = @(x, k, Sigma) k .* cos(x) .* exp(k .* (cos(Sigma) - sin(x))) ./ (1 + exp(k .* (cos(Sigma) - sin(x)))).^2;
   
for i = [1, 4]
   S = vpasolve([1.01 - sin(x) - d / (1 + exp(k * (cos(Sigma) - sin(y)))), 1.01 - sin(y) - d / (1 + exp(k * (cos(Sigma) - sin(x))))], [x, y], X_0(i, :));
    if (~isempty(S.x))  && (abs(S.x) < 100) && (abs(S.y) < 100)
        
        Sx = double(S.x);
        Sy = double(S.y);
        
        I1 = dI(Sx, k, Sigma);
        I2 = dI(Sy, k, Sigma);
        
        MatrIX = [cos(Sx), d * I2; d * I1, cos(Sy)];

        [V,D] = eig(MatrIX); 

        
         if(EquilibriaType([D(1,1), D(2,2)]) == 1)

            for j = 1:2
                
                Epsilon = 1e-5;
                
                X_0_1 = [Sx + Epsilon * V(1, j), Sy + Epsilon * V(2, j)];
                X_0_2 = [Sx - Epsilon * V(1, j), Sy - Epsilon * V(2, j)];
               
                
                if i == 1
                   Color = 'cyan';
                else if i == 4
                        Color = 'red';
                    end
                end
                
                if D(j, j) > 0
                    Tspan = [200 0];
                    
                else
                    Tspan = [0 200];
                    
                end

                [T, X] = ode23(@(T, X) AdlerSystem(X, SP), Tspan, X_0_1, opts);
                X = mod(X, 2 * pi);
                for p = 1:length(X) - 1
                   if (abs(X(p, 1) - X(p + 1, 1)) > pi)
                       X(p, 1) = NaN;
                   end
                   if (abs(X(p, 2) - X(p + 1, 2)) > pi)
                       X(p, 2) = NaN; 
                   end
                end
                plot(ax, X(:, 1), X(:, 2), 'Color', Color, 'LineWidth', LineWidth);
                
                [T, X] = ode23(@(T, X) AdlerSystem(X, SP), Tspan, X_0_2, opts);
                X = mod(X, 2 * pi);
                for p = 1:length(X) - 1
                   if (abs(X(p, 1) - X(p + 1, 1)) > pi)
                       X(p, 1) = NaN;
                   end
                   if (abs(X(p, 2) - X(p + 1, 2)) > pi)
                       X(p, 2) = NaN; 
                   end
                end
                plot(ax, X(:, 1), X(:, 2), 'Color', Color, 'LineWidth', LineWidth);
            end
            
        end
    end
 
end
   

end

function Equilibria(ax, SP)

MarkerSize = 12;

Sigma = SP.Sigma;

% Center = SP.Gamma(1) - SP.d / 2;
% if Center > 1
%     Center = pi / 2;
% elseif Center < -1
%         Center = 3 * pi / 2;
% else
%     Center = asin(Center);
% end

Center = pi/2;

D = [Center - Sigma, Center + Sigma];

Alpha = D(1);
Delta = 2 * Sigma;

k = -500;
Sigma = SP.Sigma;
d = SP.d;
gamma = SP.Gamma;

%Кол-во седловых СР

syms x y;

LineWidth = 3;

%Приближенные СР
X_0 = mod(real([Center - Sigma, Center - Sigma; ...
       Center + Sigma, Center - Sigma; ...
       Center - Sigma, Center + Sigma; ...
       Center + Sigma, Center + Sigma;
       asin(gamma(1) - d), asin(gamma(1) - d);
       pi - asin(gamma(1) - d), asin(gamma(1) - d);
       asin(gamma(1) - d), pi - asin(gamma(1) - d);
       pi - asin(gamma(1) - d), pi - asin(gamma(1) - d);]), 2*pi);

%Вспомогательная функция - производная от ф-ии связи
dI = @(x, k, Alpha, Delta) -k .* sin(x - Alpha - Delta / 2) .* exp(k .* (cos(Delta / 2) - cos(x - Alpha - Delta / 2))) ./ (1 + exp(k .* (cos(Delta / 2) - cos(x - Alpha - Delta / 2)))).^2;
   
for i = 1:8

   S = vpasolve([gamma(1) - sin(x) - d / (1 + exp(k * (cos(Delta / 2) - cos(y - Alpha - Delta / 2)))), ...
                 gamma(2) - sin(y) - d / (1 + exp(k * (cos(Delta / 2) - cos(x - Alpha - Delta / 2))))], [x, y], X_0(i, :));
    if (~isempty(S.x)) && (abs(S.x) < 100) && (abs(S.y) < 100)
        
        Sx = mod(double(S.x), 2*pi);
        Sy = mod(double(S.y), 2*pi);
        
        I1 = dI(Sx, k, Alpha, Delta);
        I2 = dI(Sy, k, Alpha, Delta);
        
        MatrIX = [cos(Sx), d * I2; d * I1, cos(Sy)];

        [V,D] = eig(MatrIX); 
        N = EquilibriaType([D(1, 1), D(2, 2)]);
       
        switch N
            case 0 % Сложное СР
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'black', 'LineWidth', LineWidth);
            case 1 % Седло
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0.6350 0.0780 0.1840], 'LineWidth', LineWidth);
            case 2 % Устойчивый узел
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', [0.6350 0.0780 0.1840], 'LineWidth', LineWidth);
            case 3 % Неустойчивый узел
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', [0 0.4470 0.7410], 'LineWidth', LineWidth);
            case 4 % Центр
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'magenta', 'MarkerEdgeColor', [0.4940 0.1840 0.5560], 'LineWidth', LineWidth);
            case 5 % Устойчивый фокус
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [1 1 0], 'MarkerEdgeColor', [0 0.4470 0.7410], 'LineWidth', LineWidth);
            case 6 % Неустойчивый фокус
                plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [1 1 0], 'MarkerEdgeColor', [0.6350 0.0780 0.1840], 'LineWidth', LineWidth);
        end
        
%         switch N
%             case 0 % Сложное СР
%                 plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0 0 0], 'LineWidth', LineWidth);
%             case 1 % Седло
%                 plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0 1 0], 'LineWidth', LineWidth);
%             case 2 % Устойчивый узел
%                 plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [0 0 1], 'LineWidth', LineWidth);
%             case 3 % Неустойчивый узел
%                 plot(ax, Sx, Sy, 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', [1 0 0], 'LineWidth', LineWidth);
%         end

    end
 
end

end

function VectorField(ax, SP)

[X,Y] = meshgrid(0 : pi/8 : 2*pi, 0 : pi/8 : 2*pi);

D = [pi/2 - SP.Sigma, pi/2 + SP.Sigma];
k = -500;
Alpha = D(1);
Delta = 2 * SP.Sigma;

U = SP.Gamma(1) - sin(X) - SP.d .* 1 ./ (1 + exp(k .* (cos(Delta / 2) - cos(Y - Alpha - Delta ./ 2))));
V = SP.Gamma(1) - sin(Y) - SP.d .* 1 ./ (1 + exp(k .* (cos(Delta / 2) - cos(X - Alpha - Delta ./ 2))));

quiver(ax, X, Y, U, V, 'Color', [0.6350 0.0780 0.1840])

end

% k = -500;
% Sigma = pi/2;
% d = 2.01;
% gamma = 1.01;
% 
% syms x y;
% 
% X_0 =  [pi/2 + Sigma, pi/2 + Sigma];
% 
% S = vpasolve([1.01 - sin(x) - d / (1 + exp(k * (cos(Sigma) - sin(y)))), 1.01 - sin(y) - d / (1 + exp(k * (cos(Sigma) - sin(x))))], [x, y], X_0);