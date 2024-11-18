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

generateLayer :: MiniLayer -> (MiniLayer -> Int -> Int -> a) -> [[a]]
generateLayer layer cellFunc = [[cellFunc layer i j | j <- [0..hny layer-1]] | i <- [0..hnx layer-1]]

generateLayerWithConfig :: MiniLayer -> LayerConfig -> (MiniLayer -> LayerConfig -> Int -> Int -> a) -> [[a]]
generateLayerWithConfig layer config cellFunc = [[cellFunc layer config i j | j <- [0..hny layer-1]] | i <- [0..hnx layer-1]]

onResults :: (t1 -> t2 -> t3) -> (t4 -> t5 -> t6 -> t7 -> t1) -> (t4 -> t5 -> t6 -> t7 -> t2) -> t4 -> t5 -> t6 -> t7 -> t3
onResults op f g a b c d = op (f a b c d) (g a b c d)