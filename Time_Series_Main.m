function Time_Series_Main(Axes, PlotParameters)

SystemParameters = PlotParameters.SystemParameters;
FigureParameters = PlotParameters.FigureParameters;
NumericalParameters = PlotParameters.NumericalParameters;

% --- Задание параметров графика --- %

InitializePlot(Axes, SystemParameters, FigureParameters, NumericalParameters);

% --- Нанесение областей D --- %

AreaOfD(Axes, SystemParameters, NumericalParameters);

% --- Нанесение траекторий --- %
Trajectories(Axes, SystemParameters, NumericalParameters);

end

function InitializePlot(ax, SP, FP, NP)

%figure;
cla(ax);
hold(ax, 'on');
grid(ax, 'on');
ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

ax.YTick = [0, pi / 2, pi, 3*pi/2, 2 * pi];
%ax.YTickLabel = {'0', '$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'};
ax.YTickLabel = {'0', '$\pi / 2$', '$\pi$', '$3\pi / 2$', '$2\pi$'};


ax.Box = 'on';

FontSize = 24;
ax.FontSize = FontSize;

ax.YLim = [0, 2 * pi];
ax.XLim = [0, NP.Tspan];

AxisLabelsFontSize = 24;

ax.YLabel.String = '$\phi_1,~\phi_2$';
ax.YLabel.FontSize = AxisLabelsFontSize;
ax.YAxisLocation = 'left';
%ax.YLabel.Rotation = 0;

ax.XLabel.String = 't';
ax.XLabel.FontSize = AxisLabelsFontSize;



end

function AreaOfD(ax, SP, NP)

Sigma = SP.Sigma;
Tspan = NP.Tspan;

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

fill(ax, [0, 0, Tspan, Tspan], [D(1), D(2), D(2), D(1)] + 2*pi, [159/255, 1, 160/255], 'EdgeColor', 'none');
fill(ax, [0, 0, Tspan, Tspan], [D(1), D(2), D(2), D(1)], [159/255, 1, 160/255], 'EdgeColor', 'none');

plot(ax, [0, Tspan], [Center, Center], 'Color', [0.1 1 0.1]);

end

function Trajectories(ax, SP, NP)

PlotStep = 1;
LineWidth = 3;

[T, X] = ode23(@(T, X) AdlerSystem(X, SP), [0, NP.Tspan], NP.InitialCondition, NP.Options);

X = mod(X, 2 * pi);
for i = 1:length(X) - 1
   if (abs(X(i, 1) - X(i + 1, 1)) > pi)
       X(i, 1) = NaN;
   end
   if (abs(X(i, 2) - X(i + 1, 2)) > pi)
       X(i, 2) = NaN; 
   end
end
plot(ax, T(1:PlotStep:end), X(1:PlotStep:end, 1), 'b', 'LineWidth', LineWidth);
plot(ax, T(1:PlotStep:end), X(1:PlotStep:end, 2), 'r', 'LineWidth', LineWidth);
end
