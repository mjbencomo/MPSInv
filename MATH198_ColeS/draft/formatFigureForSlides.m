function formatFigureForSlides(fig)
% Format a MATLAB figure that will occupy about half a 16:9 slide.

fig.Color = 'w';
fig.Units = 'inches';
pos = fig.Position;
pos(3:4) = [6.3 4.7];
fig.Position = pos;

axesHandles = findall(fig,'Type','axes');

for k = 1:numel(axesHandles)
    ax = axesHandles(k);

    % Skip colorbar axes.
    if strcmp(ax.Tag,'Colorbar')
        continue
    end

    ax.FontSize = 17;       % tick labels
    ax.LineWidth = 1.2;

    ax.Title.FontSize = 26;
    ax.Title.FontWeight = 'bold';

    ax.XLabel.FontSize = 20;
    ax.YLabel.FontSize = 20;
end

legends = findall(fig,'Type','legend');
set(legends,'FontSize',16);

colorbars = findall(fig,'Type','colorbar');
for k = 1:numel(colorbars)
    colorbars(k).FontSize = 16;
    colorbars(k).Label.FontSize = 18;
end
end