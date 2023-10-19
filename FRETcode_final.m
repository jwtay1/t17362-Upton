clearvars
close all

%Load YFP and DAPI images taken during bleaching process (for correlation
%coefficient)

%Files are sequences of 3 images for each dataset, but compiled by channel
YFPfile='10_11_BamA_YidC_live_binned_YFP.tif';
DAPIfile='10_11_BamA_YidC_live_binned_DAPI.tif';

%Compute number of frames
YFPinfo = imfinfo(YFPfile);
numFramesYFP = numel(YFPinfo);

DAPIinfo = imfinfo(DAPIfile);
numFramesDAPI = numel(DAPIinfo);

%calibration slides
calDAPI=imread('cal DAPI binned.tif');
calYFP=imread('cal YFP binned.tif');

%de-noise calibration slides
calDAPI_corr=imgaussfilt(calDAPI);
calYFP_corr=imgaussfilt(calYFP);

%normalize calibration 
calDAPI_corr=double(calDAPI_corr);
calDAPI_corr=calDAPI_corr./(max(max(calDAPI_corr)));
calYFP_corr=double(calYFP_corr);
calYFP_corr=calYFP_corr./(max(max(calYFP_corr)));


%For each frame during bleaching, mask DAPI/YFP channels bright points, add
%up cumulative intensities of those points and monitor over bleaching
YFPbleach=[];
DAPIbleach=[];
numCells_total=[];

%input number of first frames
% firstframes = [1 4 7 10 13 16 19 22 25 28]; 
% secondframes = [2 5 8 11 14 17 20 23 26 29];

firstframes = 1:3:numFramesDAPI;
secondframes = firstframes + 1;

for iFrame = 1:numFramesYFP

    %read in images sequentially
    YFPframe = imread(YFPfile,iFrame);
    DAPIframe = imread(DAPIfile,iFrame);

    %convert to double format
    YFPframe = double(YFPframe);
    DAPIframe = double(DAPIframe);

    %divide each frame by calibration image
    YFPframe = YFPframe./calYFP_corr;
    DAPIframe = DAPIframe./calDAPI_corr;

    %Make a mask
    thLvl = getThreshold(DAPIframe);

    maskDAPI = DAPIframe > thLvl;
    maskDAPI = imopen(maskDAPI, strel('diamond', 1));
    maskDAPI = imdilate(maskDAPI, strel('disk', 1));

    maskDAPI = imclearborder(maskDAPI); 
      
    %If first of three frames, set up new tracking object
    if  ismember(iFrame,firstframes)==1

        %Initialize a new tracking object
        L = LAPLinker;

        ctrFrame = 1;   
      
    end

    %Fix this - need to use morphological operations
    bg_DAPI = median(median(DAPIframe));
    bg_YFP = median(median(YFPframe));

    %Compile cell data so we can implement tracking
    data_DAPI = regionprops(maskDAPI, DAPIframe, 'PixelValues', 'Area', 'Centroid');
    data_YFP = regionprops(maskDAPI, YFPframe, 'PixelValues');

    for iCell = 1:numel(data_DAPI)

        celldata(iCell).DAPIpixelvalues = data_DAPI(iCell).PixelValues;
        celldata(iCell).YFPpixelvalues = data_YFP(iCell).PixelValues;
        celldata(iCell).Area = data_DAPI(iCell).Area;
        celldata(iCell).Centroid = data_DAPI(iCell).Centroid;

        %Background correct DAPI and YFP values
        celldata(iCell).DAPIbleach = sum(celldata(iCell).DAPIpixelvalues - bg_DAPI);
        celldata(iCell).meanDAPIbleach = sum(celldata(iCell).DAPIpixelvalues - bg_DAPI)/celldata(iCell).Area;
        
        celldata(iCell).YFPbleach = sum(celldata(iCell).YFPpixelvalues - bg_YFP);
        celldata(iCell).meanYFPbleach = sum(celldata(iCell).YFPpixelvalues - bg_YFP)/celldata(iCell).Area;
    
    end

    %Add data to the the track
    L = assignToTrack(L, ctrFrame, celldata);

    %Every third frame, compile all tracks into a single struct
    if rem(iFrame, 3) == 0

        if exist('combinedCellData', 'var')
            currNumTracks = numel(combinedCellData);
        else
            currNumTracks = 0;
        end

        for iCell = 1:L.NumTracks

            if ~exist('combinedCellData', 'var')
                combinedCellData = getTrack(L, iCell);
            else
                combinedCellData(currNumTracks + iCell) = getTrack(L, iCell);
            end            

        end
    end

    %Write some files to check tracking
    Iout = imfuse(DAPIframe, bwperim(maskDAPI));

    for iCell = 1:L.NumTracks

        ct = getTrack(L, iCell);

        if ismember(ctrFrame, ct.Frames)

            Iout = insertText(Iout, ct.Centroid(end, :), int2str(iCell), ...
                'BoxOpacity', 0, 'TextColor', 'yellow');

        end

    end

    [~, fnOut] = fileparts(YFPfile);
    imwrite(Iout, [fnOut(1:end-3), '_', int2str(iFrame), '.png']);

    ctrFrame = ctrFrame + 1;

end
%%
%Remove any cells that were not tracked the whole way or were tracked
%incorrectly (based on change in size)
idxToDelete = [];
for iCell = 1:numel(combinedCellData)

    if numel(combinedCellData(iCell).Frames) < 3
        idxToDelete = [idxToDelete iCell];

    elseif any(diff(combinedCellData(iCell).Area) > 80)
        idxToDelete = [idxToDelete iCell];
    end

end
combinedCellData(idxToDelete) = [];

figure;
plot(combinedCellData(5).Frames, combinedCellData(5).DAPIbleach, combinedCellData(5).Frames, combinedCellData(5).YFPbleach)


return


%eliminate cells that aren't in both DAPIpre and DAPIpost
%if numel(DAPIpre)>numel(DAPIpost);
   % for i=1:numel(DAPIpre);

numCells_total=(sum(numCells_total))./3;
%sort into individual experiments
r_save=[];
time=[0 2.25 4.5];
for ii=1:(numFramesYFP/3)
    YFPbleach_split=YFPbleach((1+(3*(ii-1))):(3*ii));
    DAPIbleach_split=DAPIbleach((1+(3*(ii-1))):(3*ii));
    
    figure;
    scatter(time,YFPbleach_split);
    title('YFP bleach');
    xlabel('min');
    ylabel('cumulative YFP intensity');
    
    figure;
    scatter(time,DAPIbleach_split);
    title('DAPI signal during YFP bleach');
    xlabel('min');
    ylabel('cumulative DAPI intensity');
    
    r=corrcoef(DAPIbleach_split, YFPbleach_split);
    r_save=[r_save r];
end

r_unique=unique(r_save);
r_unique(r_unique==1)=[];
r_avg=mean(r_unique);

%divide these images
Emax_avg_store=[];
Emax_std_store=[];
Emean_avg_store=[];
Emean_std_store=[];
Emean_store=[];
Emean_numTracks_store=[];
Length_store=[];

mean_pre=mean(DAPIpre);
mean_post=mean(DAPIpost);
FRETaverage=(mean_post-mean_pre)/(mean_post)*100

FRETcell=(DAPIpost-DAPIpre)./DAPIpost*100;
figure;
histogram(FRETcell);
title('Ecell (%)');


