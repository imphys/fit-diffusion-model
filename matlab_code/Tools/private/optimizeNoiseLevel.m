function currentPoint = optimizeNoiseLevel( currentPoint, opts)


if isempty(currentPoint.noiseLevel)
    currentPoint.noiseLevel = std