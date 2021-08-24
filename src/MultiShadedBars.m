function [hobj figure_handle] = MultiShadedBars(plotfile,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Default values for shading, threshold, and edge
shade_lowbound=0.05;
shade_upperbound=0.95;
showthreshold=false;
noedge=false;

%% Get input arguments for shading, threshold, edge, and output file
if ~isempty(varargin)
    for iter = 1:size(varargin,2)
        if ischar(varargin{iter})
            switch(varargin{iter})
                case('LowerBound')
                    shade_lowbound=varargin{iter+1};
                case('UpperBound')
                    shade_upperbound=varargin{iter+1};
                case('Outputfile')
                    outputfile=varargin{iter+1};
                case('ShowThreshold')
                    showthreshold=true;
                case('NoEdge')
                    noedge=true;
            end
        end
    end
end

%% Read plot file
fid = fopen(plotfile);
plot_data = textscan(fid,'%s%f%f%f%f%f','Delimiter',',');
fclose(fid);
nplots = length(plot_data{1});

%% Create plot(s)
for curr_plot=1:nplots
    [xbins,ymeans,ybound] = loadcorr(plot_data{1}{curr_plot},'LowerBound',shade_lowbound,'UpperBound',shade_upperbound);
    if curr_plot==1
        if noedge
             [hobj figure_handle] = shadedbars(xbins,ymeans,ybound,[plot_data{2}(curr_plot) plot_data{3}(curr_plot) plot_data{4}(curr_plot)],plot_data{5}(curr_plot),[],noedge);
        else
             [hobj figure_handle] = shadedbars(xbins,ymeans,ybound,[plot_data{2}(curr_plot) plot_data{3}(curr_plot) plot_data{4}(curr_plot)],plot_data{5}(curr_plot));
        end
    else
        if noedge
            shadedbars(xbins,ymeans,ybound,[plot_data{2}(curr_plot) plot_data{3}(curr_plot) plot_data{4}(curr_plot)],plot_data{5}(curr_plot),hobj,noedge);
        else
            shadedbars(xbins,ymeans,ybound,[plot_data{2}(curr_plot) plot_data{3}(curr_plot) plot_data{4}(curr_plot)],plot_data{5}(curr_plot),hobj);            
        end
    end
    if showthreshold
        hobj.mainLine = plot(hobj.mainLine.XData,zeros(size(hobj.mainLine.YData))+plot_data{6}(curr_plot),'Color',[plot_data{2}(curr_plot) plot_data{3}(curr_plot) plot_data{4}(curr_plot)],'LineStyle','--');
    end
end

%% Save plot(s) to output file
if exist('outputfile','var')
    saveas(figure_handle,outputfile);
end

