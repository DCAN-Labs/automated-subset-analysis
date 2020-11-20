function [xbins,ymeans,ybound,corrtable] = loadcorr(csvfile,varargin)
shade_lowbound=0.05;
shade_upperbound=0.95;
if isempty(varargin) == 0
    for iter = 1:size(varargin,2)
        if ischar(varargin{iter})
            switch(varargin{iter})
                case('LowerBound')
                    shade_lowbound=varargin{iter+1};
                case('UpperBound')
                    shade_upperbound=varargin{iter+1};
            end
        end
    end
end
disp(['csv file: ' csvfile])
fid = fopen(csvfile);
stuff=textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
if strcmp(stuff{1}{1},'Subjects')
    subjectcol=1;
    corrcol=2;
else
    subjectcol=2;
    corrcol=1;
end
corrtable=zeros(length(stuff{1})-1,2);
corrtable(:,1) = cellfun(@str2num,stuff{subjectcol}(2:end));
corrtable(:,2) = cellfun(@str2num,stuff{corrcol}(2:end));
xbins = transpose(unique(corrtable(:,1))); %needs transpose for shadedbars
nbins = length(xbins);
ymeans = zeros(1,nbins);
ybound = zeros(2,nbins);
for iter = 1:nbins
    ymeans(iter) = mean(corrtable(corrtable(:,1) == xbins(iter),2));
    ybound(:,iter) = quantile(corrtable(corrtable(:,1) == xbins(iter),2),[shade_lowbound shade_upperbound]);
end
