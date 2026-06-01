function matchOut = fig1_roiMatch_cleanup(matchFOVs)

matchFOVs(1).iGrp = [];

a = 0;
for ifov = 1:numel(matchFOVs)-1
    for jfov = ifov+1:numel(matchFOVs)
        if strcmpi(matchFOVs(ifov).name,matchFOVs(jfov).name) && ...
                strcmpi(matchFOVs(ifov).fov,matchFOVs(jfov).fov)
            a = a + 1;
            matchFOVs(ifov).iGrp = a; matchFOVs(jfov).iGrp = a;
        end
    end
    if isempty(matchFOVs(ifov).iGrp)
        a = a + 1;
        matchFOVs(ifov).iGrp = a;
    end
end

iGrp = cat(1,matchFOVs.iGrp);
uGrp = unique(iGrp);

% matchOut = struct;
for iOut = 1:numel(uGrp)
    fovLocal = matchFOVs(iGrp == uGrp(iOut));
    if numel(fovLocal) == 1
        matchOut(iOut) = fovLocal; %#ok<*AGROW>
    else
        nFov1 = numel(cat(1,fovLocal(1).iv1sesh,fovLocal(1).iv4_1sesh,fovLocal(1).iv4_Nsesh,fovLocal(1).iv5sesh));
        nFov2 = numel(cat(1,fovLocal(2).iv1sesh,fovLocal(2).iv4_1sesh,fovLocal(2).iv4_Nsesh,fovLocal(2).iv5sesh));
        if nFov1 > nFov2
            matchOut(iOut) = fovLocal(1);
        else
            matchOut(iOut) = fovLocal(2);
        end
    end
end

end