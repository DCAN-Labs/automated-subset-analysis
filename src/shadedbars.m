function [h figure_handle] = shadedbars(x,y,errBar,colors,alpha,hobj,noedge)

% Inputs
% x - row vector of x values 
% y - row vector of y values (1 x N) where N is number of bins (i.e.,
%     length(x))
% errBar - if a row vector we draw symmetric errorbars. If it has 2 rows, 
%          error bars are drawn around mean, with row 1 being the upper bar 
%          and row 2 being the lower bar. 
% colors - 1x3 vector containing rgb values (ranging from 0 - 1)
% alpha - transparency value between 0 - 1. 
%
% Output
%   A figure containing shaded error bars around a mean value 
%   h - figure handle 

%%
% Cheking the y data
[C, N] = size(y);
if N==1
    y = y'; % flip y to row vector if input is column vector 
end

% Cheking the x data
if isempty(x)
    error('need a vector of x values of length y')
end

% if sizes of x and y are different error out
if any(size(x) == size(y)) ~= 1  
    error('x and y do not have same dimensions')
end

% If errBar input is a vector, add and subtract it from the mean (y) value.
if size(errBar,1)==1 && size(y,2) == size(errBar,2) 
    uE = y+(errBar./2);
    lE = y-(errBar./2);
    %Calculate the y values at which we will place the error bars
elseif size(errBar,1) == 2 && size(y,2) == size(errBar,2) 
    uE=errBar(1,:);
    lE=errBar(2,:);
elseif size(errBar,1) > 2 
    error('errBar has too many rows')
else
    error('errBar and y have different number of rows')
end

%% Plot it 
if exist('hobj','var')
else
    figure_handle = figure;
end
hold on 

edgeColor=colors+(1-colors)*0.5; % edge color 
patchSaturation=0.15; %How de-saturated or transparent to make the patch
faceAlpha=alpha; % set alpha
patchColor=colors+(1-colors)*(1-patchSaturation); %combine color and saturation
set(gcf,'renderer','painters')

%Make the cordinats for the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

% plot patch 
h.patch=patch(xP,yP,1,'facecolor',patchColor,...
    'edgecolor','none',...
    'facealpha',faceAlpha);

%Make edges around the patch.
if exist('noedge','var')
    if isempty(noedge)
        h.edge(1)=plot(x,lE,'-','color',edgeColor);
        h.edge(2)=plot(x,uE,'-','color',edgeColor);
    end
else        
    h.edge(1)=plot(x,lE,'-','color',edgeColor);
    h.edge(2)=plot(x,uE,'-','color',edgeColor);
end
% plot mean line (y) 
h.mainLine=plot(x,y,'color',colors);


end
