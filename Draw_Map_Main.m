function  Draw_Map_Main(ax, Map)

% --- Инициализация графика --- %

InitializePlot(ax, Map.MP);

% --- Отрисовка карты --- %

DrawMap(ax, Map.MP, Map.ValueArray)

end

function InitializePlot(ax, MP)

% --- --- %
cla(ax);
hold(ax, 'on');
axis(ax, 'square');

% --- Границы карты --- %
ax.XLim = MP.FirstParamBorder;
ax.YLim = MP.SecondParamBorder;

% --- Подписи осей --- %
AxisLabelsFontSize = 34;
FontSize = 34;
ax.TickLabelInterpreter = 'latex';

ax.XLabel.String = '$\sigma$';
ax.XLabel.FontSize = AxisLabelsFontSize;

ax.YLabel.String = '$d$';
ax.YLabel.FontSize = AxisLabelsFontSize;

set(ax, 'FontSize', FontSize);
end

function DrawMap(ax, MP, ValueArray)

% --- Формирование сетки --- %
FirstParam = linspace(MP.FirstParamBorder(1), MP.FirstParamBorder(2), MP.Width);
SecondParam = linspace(MP.SecondParamBorder(1), MP.SecondParamBorder(2), MP.Height);
[FirstParamMesh, SecondParamMesh] = meshgrid(FirstParam, SecondParam);

% --- Отрисовка вспомогательных пикселей для фиксирования цвета --- %
ColorMin = 1;
ColorMax = 256;
pcolor(ax, [-2, -1;-2, -1], [-2, -2;-1, -1], [ColorMin, ColorMin; ColorMin, ColorMin]);
pcolor(ax, [-3, -2;-3, -2], [-2, -2;-1, -1], [ColorMax, ColorMax; ColorMax, ColorMax]);
colormap(ax, colorcube);

% --- Отрисовка сетки --- %

pc = pcolor(ax, FirstParamMesh, SecondParamMesh, ValueArray*3 + 50);
pc.EdgeColor = 'none';

%shading flat;
%shading interp

end
