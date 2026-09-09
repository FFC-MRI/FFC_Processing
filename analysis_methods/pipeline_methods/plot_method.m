clear
root = cd;
folders = dir;
datalength = 0;
for fi = 3:numel(folders)
    f = folders(fi);
    cd(root)
    % loading the data
    try
        cd(f.name)
        load fitresult.mat
    catch
        continue
    end
    datalength = datalength+1;
    outfit.patient = f.name;
    data(datalength) = outfit; %#ok<SAGROW>
    % oversampling to increase resolution
    d = [];
    for i = 1:size(data(datalength).R1,3)
        d(:,:,i) = imresize(data(datalength).R1(:,:,i),2);
    end
    data(datalength).R1 = abs(d);
end

for d = 1:length(data)
    % plotting R square maps
    figure(1)
    subplot(4,4,d)
    imagesc(real(data(d).rsquare))
    caxis([0 1])
    axis off
    
%     % plotting R1 maps
%     figure(2)
%     subplot(4,4,d)
%     imagesc(log10(data(d).R1(:,:,1)))
%     caxis(log10([1 8]))
%     axis off
%     figure(3)
%     subplot(4,4,d)
%     imagesc(log10(data(d).R1(:,:,2)))
%     caxis(log10([2 12]))
%     axis off
%     figure(4)
%     subplot(4,4,d)
%     imagesc(log10(data(d).R1(:,:,3)))
%     caxis(log10([3 15]))
%     axis off
%     figure(5)
%     subplot(4,4,d)
%     imagesc(log10(data(d).R1(:,:,4)))
%     caxis(log10([3 15]))
%     axis off
    
    % plotting R1 maps
    figure(2)
    subplot(4,4,d)
    imagesc((data(d).R1(:,:,1)))
    caxis(([1 8]))
    axis off
    figure(3)
    subplot(4,4,d)
    imagesc((data(d).R1(:,:,2)))
    caxis(([2 12]))
    axis off
    figure(4)
    subplot(4,4,d)
    imagesc((data(d).R1(:,:,3)))
    caxis(([3 15]))
    axis off
    figure(5)
    subplot(4,4,d)
    imagesc((data(d).R1(:,:,4)))
    caxis(([3 15]))
    axis off
    
    % Contrast map
    figure(6)
    subplot(4,4,d)
    imagesc(mean(log10(data(d).R1),3))
    caxis([0.3 1])
    title(strrep(data(d).patient,'_',' '))
    axis off
    title('Constrast map')
    
    % slope map
    figure(7)
    subplot(4,4,d)
    imagesc((log10(data(d).R1(:,:,end))-log10(data(d).R1(:,:,2)))/(log(2e-4)-log(0.02)));
    caxis([-0.3 0])
    title(strrep(data(d).patient,'_',' '))
    axis off
    title('Slope map')
end