function pd1 = fit_max_values(max_values_without_qburst)
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(MAX_VALUES_WITHOUT_QBURST)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 30-Sep-2020 17:38:10

% Output fitted probablility distribution: PD1

% Data from dataset "max_values_without_qburst data":
%    Y = max_values_without_qburst

% Force all inputs to be column vectors
max_values_without_qburst = max_values_without_qburst(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "max_values_without_qburst data"
[CdfF,CdfX] = ecdf(max_values_without_qburst,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 3;
BinInfo.nbins = 40;
[~,BinEdge] = internal.stats.histbins(max_values_without_qburst,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0.666667 0],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'max_values_without_qburst data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('lognormal',[ 13.15930243095, 0.5836819084205])
pd1 = fitdist(max_values_without_qburst, 'lognormal');
YPlot = pdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'fit 1';

% Adjust figure
box on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');
