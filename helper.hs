module Helper where

removeLeadingSpaces:: String -> String
removeLeadingSpaces = dropWhile (== ' ')

setToZeroIfInfinite :: Double -> Double
setToZeroIfInfinite x 
    | isInfinite x = 0.0
    | otherwise    = x

setToZeroIfNaN :: Double -> Double
setToZeroIfNaN x 
    | isNaN x   = 0.0
    | otherwise = x

setToZeroIfNanOrInfinite:: Double -> Double
setToZeroIfNanOrInfinite = setToZeroIfNaN . setToZeroIfInfinite