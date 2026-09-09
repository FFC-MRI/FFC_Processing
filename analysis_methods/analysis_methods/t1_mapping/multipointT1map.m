function [ T1map ] = multipointT1map(imagestack,TI,checkfit,mask,method)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

extra.T1Vec = linspace(20,15000,3000);
extra.T1Init = 20;
extra.tVec = TI;
imagestack = abs(imagestack);
if nargin<5
    method = 'RD-NLS-PR'; % 'RD-NLS' : non-abs
end
if nargin <3
    checkfit = false;
end
if nargin<4
    mask = ones(size(imagestack,1),size(imagestack,2));
end
if isempty(mask)
    mask = ones(size(imagestack,1),size(imagestack,2));
end

 for n=1:size(imagestack,3)
   imagestack(:,:,n) = imgaussfilt(imagestack(:,:,n));
 end

T1map = T1ScanExperiment(imagestack,extra, method,mask,checkfit);


end

