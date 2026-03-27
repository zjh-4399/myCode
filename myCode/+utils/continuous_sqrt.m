function [rootVal, lastRoot] = continuous_sqrt(z, lastRoot)
%CONTINUOUS_SQRT Continuous branch tracking for complex square roots.
rootVal = sqrt(z);
if isempty(lastRoot)
    lastRoot = rootVal;
    return;
end
if abs(rootVal - lastRoot) > abs(-rootVal - lastRoot)
    rootVal = -rootVal;
end
lastRoot = rootVal;
end
