"calculateForegroundIntensities" <-
  function(raw, xs, ys, sharpen=TRUE, k =length(xs)){
    if (sharpen){
      raw = sharpen(raw)
    }
    R = averagePixels(xs, ys, raw, k)
    R
  }
