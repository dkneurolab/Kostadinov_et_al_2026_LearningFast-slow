function [xDists, yDists] = fig2_GLMdistCalc(fovParams)
%% Rotate everybody
nRois = numel(fovParams);
medSlope = median(cat(1,fovParams.slope));
rotTheta = pi/2-atan(medSlope);
if rotTheta >= pi/2; rotTheta = rotTheta-pi; end
if rotTheta < -pi/2; rotTheta = rotTheta+pi; end
rotMat = [cos(rotTheta) -sin(rotTheta); sin(rotTheta) cos(rotTheta)];

for iroi = 1:nRois
    fovParams(iroi).ctrs = median(single(fovParams(iroi).rois));
    fovParams(iroi).rois2 = single(fovParams(iroi).rois);
    fovParams(iroi).rois2(:,2) = 512-fovParams(iroi).rois2(:,2);
    fovParams(iroi).rois2 = (rotMat*fovParams(iroi).rois2')';
end

allCtrs = cat(1,fovParams.ctrs); allCtrs(:,2) = 512-allCtrs(:,2);
allCtrsRot = (rotMat*allCtrs')';

% figure; scatter(allCtrs(:,1),allCtrs(:,2));
% hold on;
% scatter(allCtrsRot(:,1),allCtrsRot(:,2));
% 
% for iPlot = 1:20:numel(fovParams)
%     plotCol = [rand rand rand];
% 
%     plot(fovParams(iPlot).rois(:,1),512-fovParams(iPlot).rois(:,2),'Color',plotCol);
%     plot(fovParams(iPlot).rois2(:,1),fovParams(iPlot).rois2(:,2),'Color',plotCol/2);
%     
% end
% axis square;

xDists0 = repmat(allCtrsRot(:,1),1,nRois);
yDists0 = repmat(allCtrsRot(:,2),1,nRois);

xDists = abs(xDists0-xDists0');
yDists = abs(yDists0-yDists0');


end
