function unitQuat = vec2unitQuat(vec)
    if vecnorm(vec) > 1
        unitQuat = [0; vec / vecnorm(vec)];
    else
        unitQuat = [sqrt(1 - vecnorm(vec)); vec];
    end
end