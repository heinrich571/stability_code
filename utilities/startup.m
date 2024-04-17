%% -------------------------- User toggles --------------------------- %%
StartupSettings.FontName            = 'Arial';
StartupSettings.FontWeight          = 'Normal';
StartupSettings.AxesFontSize        = 18;
StartupSettings.TitleTextMultiplier = 1.4;
StartupSettings.LabelTextMultiplier = 1.2;
StartupSettings.LineWidth           = 2.5;
StartupSettings.MarkerSize          = 9;
StartupSettings.LegendLocation      = 'northeast';
StartupSettings.LegendFontSize      = 16;
StartupSettings.LegendFontWeight    = 'Normal';
StartupSettings.AxesTickLength      = [0.01 0.01];
StartupSettings.AxesTickDir         = 'in';
StartupSettings.ColorMap            = jet;
StartupSettings.GridMode            = 'on';
StartupSettings.GridMinorMode       = 'off';
StartupSettings.ColorOrder          = [ 000, 090, 255
                                        000, 180, 000
                                        255, 000, 000
                                        255, 000, 255
                                        253, 127, 001
                                        000, 190, 190
                                        176, 080, 255
                                        180, 210, 000
                                        107, 117, 150
                                        170, 170, 170]/255;
% --------------------------------------------------------------------- %

%% ---------------------- Default line settings ---------------------- %%
set(groot, 'defaultLineLineWidth', StartupSettings.LineWidth)
set(groot, 'defaultLineMarkerSize', StartupSettings.MarkerSize)
set(groot, 'defaultErrorBarLineWidth', StartupSettings.LineWidth)
% --------------------------------------------------------------------- %

%% -------------------- Default contour settings --------------------- %%
set(groot, 'defaultContourLineWidth', 1.5)
set(groot, 'defaultContourLabelSpacing', 100)
set(groot, 'defaultTextFontName', StartupSettings.FontName)
set(groot, 'defaultTextFontSize', 16)
% --------------------------------------------------------------------- %

%% --------------------- Default figure settings --------------------- %%
set(groot, 'defaultFigureColor', 0.96*ones(1,3))
set(groot, 'defaultFigureColormap', StartupSettings.ColorMap)
set(groot, 'defaultFigureWindowStyle', 'Docked')
set(groot, 'defaultFigureNextPlot', 'add')
set(groot, 'defaultFigureCreateFcn', @(varargin) figcreateufcn)
set(groot, 'defaultFigureGraphicsSmoothing', 'on')
% --------------------------------------------------------------------- %

%% ----------------------- Default UI settings ----------------------- %%
set(groot, 'defaultUicontrolFontName', StartupSettings.FontName)
set(groot, 'defaultUitableFontName', StartupSettings.FontName)
set(groot, 'defaultUipanelFontName', StartupSettings.FontName)
% --------------------------------------------------------------------- %

%% ---------------------- Default text settings ---------------------- %%
set(groot, 'defaultTextFontName', StartupSettings.FontName)
% --------------------------------------------------------------------- %

%% ---------------------- Default axes settings ---------------------- %%
set(groot, 'defaultAxesBox', 'off')
set(groot, 'defaultAxesCreateFcn', @(ax,~)set(ax.Toolbar,'Visible','off'))
set(groot, 'defaultAxesFontName', StartupSettings.FontName)
set(groot, 'defaultAxesFontWeight', StartupSettings.FontWeight)
set(groot, 'defaultAxesTitleFontWeight', StartupSettings.FontWeight)
set(groot, 'defaultAxesFontSize', StartupSettings.AxesFontSize)
set(groot, 'defaultAxesLabelFontSizeMultiplier', StartupSettings.LabelTextMultiplier)
set(groot, 'defaultAxesTitleFontSize', StartupSettings.TitleTextMultiplier)
set(groot, 'defaultAxesXColor', 0.1*ones(1,3), 'defaultAxesYColor', 0.1*ones(1,3))
set(groot, 'defaultAxesLineWidth', 1.0)
set(groot, 'defaultLineMarkerFaceColor', ones(1, 3))
set(groot, 'defaultAxesXMinorGridMode', 'manual', 'defaultAxesYMinorGridMode', 'manual', 'defaultAxesZMinorGridMode', 'manual')
set(groot, 'defaultAxesXGrid', StartupSettings.GridMode, 'defaultAxesYGrid', StartupSettings.GridMode, 'defaultAxesZGrid', StartupSettings.GridMode)
set(groot, 'defaultAxesXMinorGrid', StartupSettings.GridMinorMode, 'defaultAxesYMinorGrid', StartupSettings.GridMinorMode, 'defaultAxesZMinorGrid', StartupSettings.GridMinorMode)
set(groot, 'defaultAxesGridLineStyle', '--')
set(groot, 'defaultAxesGridAlpha', 1.0)
set(groot, 'defaultAxesMinorGridLineStyle', '-')
set(groot, 'defaultAxesMinorGridAlpha', 0.04)
set(groot, 'defaultAxesColor', 1*ones(1,3))
set(groot, 'defaultAxesTickDirMode', 'manual')
set(groot, 'defaultAxesTickDir', StartupSettings.AxesTickDir)
set(groot, 'defaultAxesTickLength', StartupSettings.AxesTickLength.*ones(1,2))
set(groot, 'defaultAxesXMinorTick', 'on', 'defaultAxesYMinorTick', 'on', 'defaultAxesZMinorTick', 'on')
set(groot, 'defaultAxesNextPlot', 'add')
set(groot, 'defaultAxesBox', 'on')
set(groot, 'defaultAxesUnits', 'Normalized')
set(groot, 'defaultAxesColorOrder', StartupSettings.ColorOrder)

set(groot, 'defaultScatterLineWidth', 1.4)

set(groot, 'defaultHistogramFaceAlpha', 1)

set(groot, 'defaultQuiverAutoScaleFactor', 1.0)
set(groot, 'defaultQuiverLineWidth', StartupSettings.LineWidth)
set(groot, 'defaultQuiverMaxHeadSize', 0.1)
% --------------------------------------------------------------------- %

%% --------------------- Default legend settings --------------------- %%
set(groot, 'defaultLegendFontSizeMode', 'manual')
set(groot, 'defaultLegendFontSize', StartupSettings.LegendFontSize)
set(groot, 'defaultLegendFontWeightMode', 'manual')
set(groot, 'defaultLegendFontWeight', StartupSettings.LegendFontWeight)
set(groot, 'defaultLegendColor', 0.99*ones(1,3))
set(groot, 'defaultLegendEdgeColor', 0.65*ones(1, 3))
set(groot, 'defaultLegendLocation', StartupSettings.LegendLocation)
% --------------------------------------------------------------------- %

%% -------------------------- Miscellaneous -------------------------- %%
set(groot, 'defaultTextarrowshapeFontName', StartupSettings.FontName)
set(groot, 'defaultTextboxshapeFontName', StartupSettings.FontName)
clearvars StartupSettings vec
clc
% --------------------------------------------------------------------- %

%% ------------------------ Required functions ----------------------- %%
function figcreateufcn()

addToolbarExplorationButtons(gcf);
toolbar = findall(gcf,'Type','UIToolbar');

if isempty(toolbar)
    return;
end

% copy figure keyboard shortcut + button
set(findall(gcf,'Type','UIMenu','Tag','figMenuEditCopyFigure'),'Accelerator','X')
copyFigureButton = uipushtool(toolbar);                                 % Adding 'Copy Figure' button to the default toolbar
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','HDF_point.gif'));                       % A nice (representative) image for the new 'Copy Figure' button
copyFigureButtonImage = ind2rgb(img,map);
copyFigureButton.CData = copyFigureButtonImage;                        	% Setting the above icon to the 'Copy Figure' button
copyFigureButton.TooltipString = 'Copy Figure';                        	% Tooltip string for 'Copy Figure'
copyFigureButton.ClickedCallback = @copyFigure;                         % Directing to the actual function that allows the user to copy the current figure

% curve fitting button
CurveFittingButton = uipushtool(toolbar);
[img,~] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','plotpicker-andrewsplot.png'));
CurveFittingButton.CData = img;
CurveFittingButton.TooltipString = 'Curve Fitting';
CurveFittingButton.ClickedCallback = @curveFitting;

% save figure as .fig file in the current directory
SaveFigHereButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','HDF_pointfieldset.gif'));
SaveFigHereButtonImage = ind2rgb(img,map);
SaveFigHereButton.CData = SaveFigHereButtonImage;
SaveFigHereButton.TooltipString = 'Save .fig in the current folder';
SaveFigHereButton.ClickedCallback = @saveAsFigAtCurrentFolder;

% save figure as .png file in the current directory
SavePNGHereButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','HDF_rasterimage.gif'));
SavePNGHereButtonImage = ind2rgb(img,map);
SavePNGHereButton.CData = SavePNGHereButtonImage;
SavePNGHereButton.TooltipString = 'Save .png in the current folder';
SavePNGHereButton.ClickedCallback = @saveAsPNGAtCurrentFolder;

% toggle x axis logarithmic scale
LogXButton = uitoggletool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','plotpicker-semilogx.png'));
LogXButton.CData = img;
LogXButton.TooltipString = 'Set Y axis to logarithmic scale';
LogXButton.ClickedCallback = @logXAxis;

% toggle y axis logarithmic scale
LogYButton = uitoggletool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','plotpicker-semilogy.png'));
LogYButton.CData = img;
LogYButton.TooltipString = 'Set Y axis to logarithmic scale';
LogYButton.ClickedCallback = @logYAxis;

% toggle equal axes
EqualAxesButton = uitoggletool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','pageicon.gif'));
EqualAxesButtonImage = ind2rgb(img,map);
EqualAxesButton.CData = EqualAxesButtonImage;
EqualAxesButton.TooltipString = 'Set axes equal';
EqualAxesButton.ClickedCallback = @axisEqual;

% change axes axes labels and title
setFigureForDocumentButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','tool_text.gif'));
setFigureForDocumentButtonImage = ind2rgb(img,map);
setFigureForDocumentButton.CData = setFigureForDocumentButtonImage;
setFigureForDocumentButton.TooltipString = 'Label settings';
setFigureForDocumentButton.ClickedCallback = @setPlotText;

% set axes limits
AxesLimitsButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','tool_double_arrow.gif'));
AxesLimitsButtonImage = ind2rgb(img,map);
AxesLimitsButton.CData = AxesLimitsButtonImage;
AxesLimitsButton.TooltipString = 'Set axes limits';
AxesLimitsButton.ClickedCallback = @axisLimits;

% quickly create a linear line between two points
computeSlopeButton = uitoggletool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','tool_line.gif'));
computeSlopeButtonImage = ind2rgb(img,map);
computeSlopeButton.CData = computeSlopeButtonImage;
computeSlopeButton.TooltipString = 'Plot slope line';
computeSlopeButton.ClickedCallback = @computeSlope;

% set data all data cursors' design for documents
setFigureForDocumentButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','book_sim.gif'));
setFigureForDocumentButtonImage = ind2rgb(img,map);
setFigureForDocumentButton.CData = setFigureForDocumentButtonImage;
setFigureForDocumentButton.TooltipString = 'Enlarge data cursors and adjust the font to match the figure''s font';
setFigureForDocumentButton.ClickedCallback = @enlargeDataCursors;

% Prepare the entire figure for a document
setFigureForDocumentButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','book_mat.gif'));
setFigureForDocumentButtonImage = ind2rgb(img,map);
setFigureForDocumentButton.CData = setFigureForDocumentButtonImage;
setFigureForDocumentButton.TooltipString = 'Prepare the figure for official documents';
setFigureForDocumentButton.ClickedCallback = @prepFigureForDocument;

% use black lines and symbols instead of colors
setFigureForDocumentButton = uipushtool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','unknownicon.gif'));
setFigureForDocumentButtonImage = ind2rgb(img,map);
setFigureForDocumentButton.CData = setFigureForDocumentButtonImage;
setFigureForDocumentButton.TooltipString = 'Use symbols';
setFigureForDocumentButton.ClickedCallback = @setPlotSymbols;

% toggle minor grid
toggleMinorGridButton = uitoggletool(toolbar);
[img,map] = imread(fullfile(matlabroot,...
    'toolbox','matlab','icons','HDF_grid.gif'));
toggleMinorGridButtonImage = ind2rgb(img,map);
toggleMinorGridButton.CData = toggleMinorGridButtonImage;
toggleMinorGridButton.TooltipString = 'Toggle grid';
toggleMinorGridButton.ClickedCallback = @toggleGrid;

% supporting functions
    function copyFigure(~,~)
        print('-clipboard','-dmeta')
    end

    function curveFitting(~,~)
        cftool
    end

    function saveAsFigAtCurrentFolder(~,~)
        currentFigure = gcf;
        if isempty(currentFigure.Name)
            figureName = 'Nameless_Figure';
        else
            figureName = currentFigure.Name;
        end
        saveas(currentFigure,[figureName,'.fig'])
    end

    function saveAsPNGAtCurrentFolder(~,~)
        currentFigure = gcf;
        if isempty(currentFigure.Name)
            imageName = 'Nameless_Figure';
        else
            imageName = currentFigure.Name;
        end
        saveas(currentFigure,[imageName,'.png'])
    end

    function logXAxis(src,~)
        currentAxes = gca;
        state = src.State;
        if strcmp(state,'on')
            set(currentAxes,'XScale','log')
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','plotpicker-semilogx.png'));
            LogXButton.CData = img;
            LogXButton.TooltipString = 'Set Y axis to linear scale';
        else
            set(currentAxes,'XScale','linear')
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','plotpicker-semilogx.png'));
            LogXButton.CData = img;
            LogXButton.TooltipString = 'Set Y axis to logarithmic scale';
        end
        grid minor
    end

    function logYAxis(src,~)
        currentAxes = gca;
        state = src.State;
        if strcmp(state,'on')
            set(currentAxes,'YScale','log')
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','plotpicker-semilogy.png'));
            LogYButton.CData = img;
            LogYButton.TooltipString = 'Set Y axis to linear scale';
        else
            set(currentAxes,'YScale','linear')
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','plotpicker-semilogy.png'));
            LogYButton.CData = img;
            LogYButton.TooltipString = 'Set Y axis to logarithmic scale';
        end
        grid minor
    end

    function axisEqual(src,~)
        state = src.State;
        if strcmp(state,'on')
            axis equal
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','pageicon.gif'));
            EqualAxesButtonImage = ind2rgb(img,map);
            EqualAxesButton.CData = EqualAxesButtonImage;
            EqualAxesButton.TooltipString = 'Set axes normal';
        else
            axis normal
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','pageicon.gif'));
            EqualAxesButtonImage = ind2rgb(img,map);
            EqualAxesButton.CData = EqualAxesButtonImage;
            EqualAxesButton.TooltipString = 'Set axes equal';
        end
    end

    function axisLimits(~,~)
        axes = gca;
        prompt = {'X axis lower limit:','X axis upper limit:','Y axis lower limit:','Y axis upper limit:','Z axis lower limit:','Z axis upper limit:'};
        dlgtitle = 'Axes limits';
        dims = [1 50];
        currentLimits = [axes.XLim,axes.YLim,axes.ZLim];
        definput = {num2str(currentLimits(1)),num2str(currentLimits(2)),num2str(currentLimits(3)),num2str(currentLimits(4)),num2str(currentLimits(5)),num2str(currentLimits(6))};
        opts = 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        if isempty(answer)
            return
        else
            newLimits = str2double(answer)';
            xlim(newLimits(1:2));
            ylim(newLimits(3:4));
            zlim(newLimits(5:6));
        end
    end

    function computeSlope(src,~)
        state = src.State;
        axes = gca;
        if strcmp(state,'on')
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','tool_line.gif'));
            computeSlopeButtonImage = ind2rgb(img,map);
            computeSlopeButton.CData = computeSlopeButtonImage;
            computeSlopeButton.TooltipString = 'Hide slope line';
            
            d = datacursormode(gcf);
            points = getCursorInfo(d);
            if isempty(points)
                return
            elseif numel(points) == 1
                return
            end
            p0 = points(end-1).Position;
            p1 = points(end).Position;
            x0 = p0(1);
            y0 = p0(2);
            x1 = p1(1);
            y1 = p1(2);
            slope = (y1-y0)/(x1-x0);
            currentAxesLimits = [axes.XLim,axes.YLim];
            x = currentAxesLimits(1):diff(currentAxesLimits(1:2))/100:currentAxesLimits(2);
            y = y0 + slope*(x-x0);
            figure(gcf);
            plot(x,y,'DisplayName','Slope line');
            fprintf('\nLine equation: y(x) = %f*x + %f\n', slope, y0-slope*x0);
            fprintf('\nSlope = %f\n\n', slope);
        else
            [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','tool_line.gif'));
            computeSlopeButtonImage = ind2rgb(img,map);
            computeSlopeButton.CData = computeSlopeButtonImage;
            computeSlopeButton.TooltipString = 'Plot slope line';
            lines = findobj(axes,'Type','Line');
            if strcmp(lines(1).DisplayName,'Slope line')
                delete(lines(1));
            end
        end
    end
    
    function enlargeDataCursors(~,~)
        AllDataCursors = findall(gcf,'type','hggroup');
        set(AllDataCursors,'FontSize',14)
        set(AllDataCursors,'FontName','Arial')
        set(AllDataCursors,'FontWeight','Normal')
    end

    function prepFigureForDocument(~,~)
        AllDataCursors = findall(gcf,'type','hggroup');
        set(AllDataCursors,'FontSize',14)
        set(AllDataCursors,'FontName','Arial')
        set(AllDataCursors,'FontWeight','Normal')
        currax = gca;
        set(currax,'FontName','Arial')
        set(currax,'FontWeight','normal')
        set(currax,'TickLabelInterpreter','latex')
        currax.Title.Interpreter  = 'latex';
        currax.XLabel.Interpreter = 'latex';
        currax.YLabel.Interpreter = 'latex';
        currax.Legend.Interpreter = 'latex';
        XLabel = get(currax,'XLabel');
        YLabel = get(currax,'YLabel');
        XLabel = XLabel.String;
        YLabel = YLabel.String;
        ix = strfind(XLabel,'[');
        iy = strfind(YLabel,'[');
        if ~isempty(ix)
            xlabel(['$\it{' XLabel(1:ix-1) '}$\rm{' XLabel(ix:end) '}'])
        else
            xlabel(['$\it{' XLabel '}$'])
        end
        if ~isempty(iy)
            ylabel(['$\it{' YLabel(1:iy-1) '}$\rm{' YLabel(iy:end) '}'])
        else
            ylabel(['$\it{' YLabel '}$'])
        end
    end

    function setPlotText(~,~)
        axes = gca;
        prompt = {'\fontsize{10} Title text:','\fontsize{10} X label text:','\fontsize{10} Y label text:','\fontsize{10} Z label text:'};
        dlgtitle = 'Axes text';
        dims = [1 50];
        Title = axes.Title;
        XLabel = axes.XLabel;
        YLabel = axes.YLabel;
        ZLabel = axes.ZLabel;
        definput = {Title.String,XLabel.String,YLabel.String,ZLabel.String};
        opts.Resize = 'on'; % So the window is resizable
        opts.WindowStyle = 'normal'; % So that the user can still interact with the axes
        opts.Interpreter = 'tex';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        if isempty(answer)
            return
        else
            Title.String  = answer{1};
            XLabel.String = answer{2};
            YLabel.String = answer{3};
            ZLabel.String = answer{4};
        end
        disp(['title(''' Title.String ''')'])
        disp(['xlabel(''' XLabel.String ''')'])
        disp(['ylabel(''' YLabel.String ''')'])
        disp(['zlabel(''' ZLabel.String ''')'])
    end

    function setPlotSymbols(~,~)
        symbolSet = {'o','s','^','d','x','v','p','>','+','<'};
        lines = findobj(gcf,'type','line');
        N_lines = numel(lines);
        lineSet = {'-','-.','--'};
        for n = N_lines:-1:1
            XLim(1) = lines(n).XData(1);
            XLim(2) = lines(n).XData(end);
            if mod(n,2) ~= 0
                x = [XLim(1),XLim(1):(XLim(2)-XLim(1))/10:XLim(2),XLim(2)];
            else
                x = [XLim(1),XLim(1):(XLim(2)-XLim(1))/13:XLim(2),XLim(2)];
            end
            lineColor = lines(n).Color;
            if N_lines <= 3
                lines(n).LineStyle = lineSet{N_lines-n+1};
            end
            y = interp1(lines(n).XData,lines(n).YData,x);
            p(N_lines-n+1) = plot(x,y,[symbolSet{N_lines-n+1}],'Color',lineColor,'DisplayName',lines(n).DisplayName);
        end
        legend(p)
    end

    function toggleGrid(src,~)
        state = src.State;
        if strcmp(state, 'on')
            grid minor
        else
            grid off
            grid on
        end
    end
end
% --------------------------------------------------------------------- %