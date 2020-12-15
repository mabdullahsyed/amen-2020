close all
clear all
startup
currdir = pwd;

% This line will open a window to ask you for the folder with images
imgdir = uigetdir();


cd(imgdir)
imgfiles = dir('*.tif');
%% Key parameters:
thresh_ki67 = 60; %signal >60 will be considered ki67 positive
overlap_pos_threshold_gfp = 0.5; % cell will be considered gfp positive if 50% of it overlaps with gfp
overlap_pos_threshold_ki67 = 0.3; % cell will be considered ki67 positive if 50% of it overlaps with ki67


%% Less key parameters
blur_dapi = 10; % For contrast normalization. If value is x, dapi channel will be gaussian blurred to radius imagesize/x.
blu_gfp = 50; % For contrast normalization. If value is x, gfp channel will be gaussian blurred to radius imagesize/x.





%%
%Directory where results csvs will be saved
dir_results = strcat(imgdir, '\results');

%directory where segmented images of each channel are saved
dir_segmented = strcat(imgdir, '\segmentedimgs');

percentdoublepos = zeros(size(imgfiles,1),1);


%loops once for each image
for i = 1:size(imgfiles,1)
    
    
    %Read all three channels from the image
    clear currimg dapi gfp ki67
    currimg = imread(imgfiles(i).name);
    dapi = currimg(:,:,3);
    gfp = currimg(:,:,2);
    % ki67 is immediately thresholded using an arbitrary value of 120
    ki67 = currimg(:,:,1)>thresh_ki67; 
    
    %% Threshold nuclei
    
    nuc_blur = gaussf(dapi,2);
    
    %This is local contrast normalization. Divides each pixel by the
    %minimum value of neighboring pixels + 20% of max value (without the
    %20% 
    dapi_local = nuc_blur./(0.2*max(nuc_blur(:))+(im2mat(gaussf(nuc_blur,size(nuc_blur,1)/blur_dapi))));
    
    %Threshold of the local contrast normalized image
    nucthresh = threshold(dapi_local);
    
    
    %% Segment nuclei using a seeded watershed
    
    %choose seeds
    nuc_ids = dilation(maxima(nuc_blur),3);
    
    %create image for watershed. Original intensities except negative areas
    %are converted to +inf so the algorithm ignores them
    nuc2 = nuc_blur;
    nuc2(~nucthresh)=inf;
    
    %seeded watershed command creates boundaries between watershed regions
    nuc_ws = waterseed(nuc_ids,-1*nuc2,1);
    
    %apply watershed to a thresholded image and label features
    nucthresh2 = nucthresh;
    nucthresh2(nuc_ws) = false;
    labelledimg = label(nucthresh2,1);
    
    %% Threshold GFP
    gfp2 = gaussf(gfp,2); %blur gfp
    gfp_local = gfp2./(0.2*max(gfp2)+(im2mat(gaussf(gfp2,size(gfp,1)/blu_gfp)))); %contrast normalize of blurred gfp
    gfpthresh = threshold(gfp_local); %threshold gfp
       
    %% Measure
    gfp_measure = measure(labelledimg, im2mat(gfpthresh),'mean',[], 1, 100, 0);
    ki67_measure = measure(labelledimg, ki67*1,'mean',[], 1, 100, 0);
    
    ID = gfp_measure.ID';
    gfp_data = gfp_measure.mean';
    ki67_data = ki67_measure.mean';
    
    
    %% Post analysis
    clear posgfp poski67 doublepos T_results
    posgfp = gfp_data>overlap_pos_threshold_gfp; 
    poski67 = ki67_data>overlap_pos_threshold_ki67;
    
    doublepos = posgfp & poski67;
    
    T_results = table(ID, gfp_data, ki67_data, posgfp, poski67, doublepos);
    
    cd(dir_results)
    filename_cells = strcat(imgfiles(i).name(1:end-4), '-results_per_cell.csv');
    writetable(T_results,filename_cells);
    
    percentdoublepos(i) = sum(doublepos)/sum(posgfp);
    cd(imgdir)
    
    %% Save segmented images
    finallabelledimg = labelledimg;
    finallabelledimg(~ismember(im2mat(labelledimg),ID))=false;
    
    cd(dir_segmented)
    writeim(finallabelledimg>0,strcat(imgfiles(i).name(1:end-4), '-nucthresh.tif'));
    writeim(gfpthresh,strcat(imgfiles(i).name(1:end-4), '-gfpthresh.tif'));
    writeim(ki67,strcat(imgfiles(i).name(1:end-4), '-ki67thresh.tif'));
    
    writeim(dapi_local>0,strcat(imgfiles(i).name(1:end-4), '-dapi.tif'));
    writeim(gfp_local,strcat(imgfiles(i).name(1:end-4), '-gfp.tif'));
    writeim(ki67,strcat(imgfiles(i).name(1:end-4), '-ki67.tif'));
    
    cd(imgdir)
end

%%
names = cat(1,{imgfiles(:).name})';
T = table(names, percentdoublepos);

filename_results = strcat('combinedresults_',imgdir((find(imgdir=='\', 1, 'last')+1):end),'.csv');
writetable(T,filename_results);
cd(currdir);
