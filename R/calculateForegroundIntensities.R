"calculateForegroundIntensities" <-
  function(raw, xs, ys, sharpen=TRUE, k =49777){
    if (sharpen){
      raw = sharpen(raw)
    }
    R = averagePixels(xs, ys, raw, k)
    R
  }
