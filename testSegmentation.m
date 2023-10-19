clearvars
clc

DAPIfile='10_11_BamA_YidC_live_binned_DAPI.tif';

I1 = imread(DAPIfile, 1);
I2 = imread(DAPIfile, 2);
I3 = imread(DAPIfile, 3);

I = I2;

thLvl = getThreshold(I);

mask = I > thLvl;
mask = imopen(mask, strel('diamond', 1));
mask = imdilate(mask, strel('disk', 1));

mask = bwareaopen(mask, 10);

imshowpair(I, bwperim(mask))



