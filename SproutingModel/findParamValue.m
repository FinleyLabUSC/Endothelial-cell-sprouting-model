function parameterVal = findParamValue(meanConc,multiplier)
% a = meanConc/multiplier;
% b = meanConc*multiplier;

b = meanConc;
a = meanConc/multiplier;

parameterVal = b + a.*randn(1,1);
return