function thLvl = getThreshold(imageIn)
%GETTHRESHOLD  Get a threshold for the image
%
%  T = GETTHRESHOLD(I) gets a greyscale threshold
%  level T for the image I.
%
%  Threshold is determined by looking at image histogram, then
%  looking for the greyscale value where the maximum count
%  drops to at least 10%.

%Get the image intensity histogram
binEdges = linspace(0,double(max(imageIn(:))),200);
[nCnts, binEdges] = histcounts(imageIn(:),binEdges);
binCenters = diff(binEdges) + binEdges(1:end-1);

nCnts = smooth(nCnts,5);

%Find the background peak count
[bgCnt,bgLoc] = findpeaks(nCnts,'Npeaks',1,'SortStr','descend');

%Find where the histogram counts drops to at least 10% of this value
thLoc = find(nCnts(bgLoc:end) <= bgCnt * 0.1, 1, 'first');

if isempty(thLoc)
    thLvl = 0;
    warning('getThreshold:CouldNotGetThreshold',...
        'Auto thresholding failed to find a suitable threshold level. Try specifying one manually.');
end

thLvl = binCenters(thLoc + bgLoc);


end