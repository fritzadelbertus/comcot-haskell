module Helper where
import Constants
import TypeModule

removeLeadingSpaces:: String -> String
removeLeadingSpaces = dropWhile (== ' ')

setToSomethingIfInfinite:: Double -> Double -> Double
setToSomethingIfInfinite y x
    | isInfinite x = y
    | otherwise = x

setToSomethingIfNan:: Double -> Double -> Double
setToSomethingIfNan y x 
    | isNaN x = y
    | otherwise = x

setToZeroIfInfinite :: Double -> Double
setToZeroIfInfinite = setToSomethingIfInfinite zero

setToZeroIfNaN :: Double -> Double
setToZeroIfNaN = setToSomethingIfNan zero

setToZeroIfNanOrInfinite:: Double -> Double
setToZeroIfNanOrInfinite = setToZeroIfNaN . setToZeroIfInfinite

rxLayer :: LayerConfig -> Double
rxLayer config = dt config / dx config

ryLayer :: LayerConfig -> Double
ryLayer = rxLayer

grxLayer :: LayerConfig -> Double
grxLayer l = grav * rxLayer l

gryLayer :: LayerConfig -> Double
gryLayer l = grav * ryLayer l