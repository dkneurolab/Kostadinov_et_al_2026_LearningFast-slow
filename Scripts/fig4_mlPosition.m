function seshOut = fig4_mlPosition(seshIn, matchData)

seshOut = seshIn;

for iSesh = 1:numel(seshIn)
    seshLocal = seshIn(iSesh);
    for jMatch = 1:numel(matchData)
        if strcmpi(seshLocal.name,matchData(jMatch).name) && strcmpi(seshLocal.fov, matchData(jMatch).fov)
            if strcmpi(seshLocal.date, matchData(jMatch).date1)
                iDate = 1;
            elseif strcmpi(seshLocal.date, matchData(jMatch).date2)                
                iDate = 2;
            else
                continue
            end
            
            % Get median slope of FOV            
            medSlope = median(cat(1,matchData(jMatch).imData(iDate).final.slope));
            medTheta = atan(1/medSlope);
            if medTheta < 0; medTheta = medTheta+2*pi(); end
            % Get original center coordinates
            ctrs = cat(1,matchData(jMatch).imData(iDate).final.ctrs);
            ctrs(:,2) = 512-ctrs(:,2);
            [xctrs, ~] = rotatePoint(ctrs(:,1),ctrs(:,2), double(medTheta));
            seshOut(iSesh).xCtrsPx = xctrs-min(xctrs);

        end
    end
end