function fig = plotVonKarmanDragSurface(options)
% plotVonKarmanDragSurface 3-D surface of von Karman nose cone drag coefficient over Mach number and L/D.
%
%   plotVonKarmanDragSurface() plots the zero-lift drag coefficient of a von Karman
%   (LD-Haack) nose cone against Mach number and fineness ratio L/D at sea level.
%
%   fig = plotVonKarmanDragSurface(Name=Value) returns the figure and accepts
%       Mach             - free-stream Mach numbers (default linspace(0.1, 3, 291))
%       FinenessRatio    - L/D values (default linspace(2, 8, 121))
%       BaseDiameter     - nose base diameter in m (default 0.1)
%       SurfaceRoughness - surface roughness height in m (default 60e-6)
%       Temperature      - static air temperature in K (default 288.15)
%       Pressure         - static air pressure in Pa (default 101325)
%       ExportPath       - image file to write, "" to skip (default "")
%
%   example
%       plotVonKarmanDragSurface(Mach=linspace(0.5, 2.5, 201), ExportPath="drag.png")
%
%   see also vonKarmanDragCoefficient

    arguments
        options.Mach (1, :) double {mustBeFinite, mustBePositive} = linspace(0.1, 3, 291)
        options.FinenessRatio (1, :) double {mustBeFinite, mustBePositive} = linspace(2, 8, 121)
        options.BaseDiameter (1, 1) double {mustBeFinite, mustBePositive} = 0.1
        options.SurfaceRoughness (1, 1) double {mustBeFinite, mustBeNonnegative} = 60e-6
        options.Temperature (1, 1) double {mustBeFinite, mustBePositive} = 288.15
        options.Pressure (1, 1) double {mustBeFinite, mustBePositive} = 101325
        options.ExportPath (1, 1) string = ""
    end

    mach = unique(options.Mach);
    finenessRatio = unique(options.FinenessRatio);
    minimumPoints = 2;
    if (numel(mach) < minimumPoints) || (numel(finenessRatio) < minimumPoints)
        error("plotVonKarmanDragSurface:gridTooSmall", ...
            "Mach and FinenessRatio each need at least %d distinct values to form a surface.", ...
            minimumPoints);
    end

    % drag over the grid; rows follow L/D, columns follow mach
    dragOptions = rmfield(options, ["Mach", "FinenessRatio", "ExportPath"]);
    dragArguments = namedargs2cell(dragOptions);
    [totalCd, ~, ~, speedOfSound] = vonKarmanDragCoefficient(mach, finenessRatio(:), dragArguments{:});
    [machGrid, finenessGrid] = meshgrid(mach, finenessRatio);
    maxCd = max(totalCd, [], "all");

    % the L/D edge with more drag is drawn at the back so the surface rises away from the camera
    if mean(totalCd(1, :)) > mean(totalCd(end, :))
        finenessDirection = "reverse";
        frontRow = numel(finenessRatio);
    else
        finenessDirection = "normal";
        frontRow = 1;
    end

    % single-hue sequential ramp (light = low drag) and neutral chart chrome
    rampHex = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", ...
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"];
    colormapSize = 256;
    rampColors = validatecolor(rampHex, "multiple");
    dragColormap = interp1(linspace(0, 1, numel(rampHex)), rampColors, linspace(0, 1, colormapSize));
    surfaceColor = validatecolor("#fcfcfb");
    primaryInk = validatecolor("#0b0b0b");
    secondaryInk = validatecolor("#52514e");
    mutedInk = validatecolor("#898781");
    gridColor = validatecolor("#e1e0d9");
    fontSize = 11;
    labelOffsetFraction = 0.04;
    sonicMach = 1;

    fig = figure(Color=surfaceColor, Position=[100, 100, 1000, 700]);
    ax = axes(fig, Color=surfaceColor, FontSize=fontSize, XColor=secondaryInk, YColor=secondaryInk, ...
        ZColor=secondaryInk, GridColor=gridColor, GridAlpha=1, TickDir="out");
    hold(ax, "on")

    dragSurface = surf(ax, mach, finenessRatio, totalCd, EdgeColor="none", FaceColor="interp");

    % iso-lines drawn at full resolution so they follow the surface
    numIsoLines = 13;
    rowIndex = unique(round(linspace(1, numel(finenessRatio), numIsoLines)));
    columnIndex = unique(round(linspace(1, numel(mach), numIsoLines)));
    plot3(ax, machGrid(rowIndex, :).', finenessGrid(rowIndex, :).', totalCd(rowIndex, :).', ...
        Color=mutedInk, LineWidth=0.5);
    plot3(ax, machGrid(:, columnIndex), finenessGrid(:, columnIndex), totalCd(:, columnIndex), ...
        Color=mutedInk, LineWidth=0.5);

    % sonic line at mach 1
    if (sonicMach >= mach(1)) && (sonicMach <= mach(end))
        sonicLineMach = sonicMach*ones(size(finenessRatio));
        sonicCd = interp2(mach, finenessRatio, totalCd, sonicLineMach, finenessRatio);
        plot3(ax, sonicLineMach, finenessRatio, sonicCd, Color=primaryInk, LineWidth=2);
        % label below the front end of the line, where no surface can hide it
        labelOffset = labelOffsetFraction*maxCd;
        text(ax, sonicMach, finenessRatio(frontRow), sonicCd(frontRow) - labelOffset, "Mach 1", ...
            Color=primaryInk, FontSize=fontSize, HorizontalAlignment="center", VerticalAlignment="top");
    end
    hold(ax, "off")

    % data tips report mach number, velocity, L/D and C_D
    dragSurface.DataTipTemplate.DataTipRows = [ ...
        dataTipTextRow("Mach", "XData", "%.2f"), ...
        dataTipTextRow("Velocity (m/s)", machGrid*speedOfSound, "%.0f"), ...
        dataTipTextRow("L/D", "YData", "%.2f"), ...
        dataTipTextRow("C_D", "ZData", "%.4f")];

    colormap(ax, dragColormap)
    clim(ax, [0, maxCd])
    colorBar = colorbar(ax, Color=secondaryInk, TickDirection="out");
    colorBar.Label.String = "Drag coefficient C_D";

    grid(ax, "on")
    xlim(ax, [mach(1), mach(end)])
    ylim(ax, [finenessRatio(1), finenessRatio(end)])
    zlim(ax, [0, maxCd])
    xlabel(ax, "Mach number", Color=primaryInk)
    ylabel(ax, "Fineness ratio L/D", Color=primaryInk)
    zlabel(ax, "Drag coefficient C_D", Color=primaryInk)
    title(ax, "Von Kármán nose cone drag coefficient", Color=primaryInk)
    subtitleFormat = "T = %.2f K, p = %.3f kPa, a = %.1f m/s, D = %g mm, roughness %g µm, " ...
        + "no base drag";
    subtitle(ax, compose(subtitleFormat, options.Temperature, options.Pressure/1e3, speedOfSound, ...
        options.BaseDiameter*1e3, options.SurfaceRoughness*1e6), Color=secondaryInk)
    ax.YDir = finenessDirection;
    view(ax, -37.5, 30)

    % fixed layout keeps the last x tick label clear of the colorbar
    axesPosition = [0.09, 0.12, 0.72, 0.74];
    colorbarPosition = [0.87, axesPosition(2), 0.02, axesPosition(4)];
    ax.InnerPosition = axesPosition;
    colorBar.Position = colorbarPosition;

    if strlength(options.ExportPath) > 0
        % hide the axes toolbar while exporting so it is not captured in the image
        ax.Toolbar.Visible = "off";
        exportgraphics(fig, options.ExportPath, Resolution=200, BackgroundColor=surfaceColor)
        ax.Toolbar.Visible = "on";
    end
end
