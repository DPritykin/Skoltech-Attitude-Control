function quat = triad(sunOrb, sunBody, bOrb, bBody)
    eSunOrb = sunOrb / vecnorm(sunOrb);
    eSunBody = sunBody / vecnorm(sunBody);
    eBOrb = bOrb / vecnorm(bOrb);
    eBBody = bBody / vecnorm(bBody);

    exBody = eSunBody;
    exOrb = eSunOrb;
    
    ezBody = crossProduct(exBody, eBBody);
    ezBody = ezBody / vecnorm(ezBody);

    ezOrb = crossProduct(exOrb, eBOrb);
    ezOrb = ezOrb / vecnorm(ezOrb);
    
    eyBody = crossProduct(ezBody, exBody);
    eyOrb = crossProduct(ezOrb, exOrb);

    quat = dcm2quat(([exBody, eyBody, ezBody] * [exOrb, eyOrb, ezOrb]'));
end