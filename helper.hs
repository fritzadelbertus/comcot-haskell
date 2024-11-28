module Helper where
import Constants
import TypeModule
import Numeric (showFFloat)

show4:: Double -> String
show4 = flip (showFFloat (Just 4)) ""

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

generateLayerWithConfig :: MiniLayer -> b -> (MiniLayer -> b -> Int -> Int -> a) -> [[a]]
generateLayerWithConfig layer config cellFunc = [[cellFunc layer config i j | j <- [0..hny layer-1]] | i <- [0..hnx layer-1]]

onResults3 :: (t1 -> t2 -> t3) -> (t4 -> t5 -> t6 -> t1) -> (t4 -> t5 -> t6 -> t2) -> t4 -> t5 -> t6 -> t3
onResults3 op f g a b c  = op (f a b c) (g a b c)

onResults2 :: (t1 -> t2 -> t3) -> (t4 -> t5 -> t1) -> (t4 -> t5 -> t2) -> t4 -> t5 -> t3
onResults2 op f g a b = op (f a b) (g a b)

onResults1 :: (t1 -> t2 -> t3) -> (t4 -> t1) -> (t4 -> t2) -> t4 -> t3
onResults1 op f g a = op (f a) (g a)


showLayerId :: LayerConfig -> String
showLayerId = show . layer_id

showCoorSystem :: LayerConfig -> String
showCoorSystem config
    | isCartesian config = "USE CARTESIAN COORDINATE SYSTEM"
    | otherwise = "NOT AVAILABLE YET"

showGoverningEq :: LayerConfig -> String
showGoverningEq config
    | laygov config == 0 = "USE LINEAR SHALLOW WATER EQUATIONS"
    | otherwise = "NOT AVAILABLE YET"








isCartesian:: LayerConfig -> Bool
isCartesian l = laycord l == 1 

isSpherical:: LayerConfig -> Bool
isSpherical l = laycord l == 0


------------------------------------------

showTotalRunTime :: GeneralConfig -> String
showTotalRunTime = show4 . total_time

showTimeInterval :: GeneralConfig -> String
showTimeInterval = show4 . time_interval

showBoundaryCodition :: GeneralConfig -> String 
showBoundaryCodition config
    | boundary_condition config == 0 = "OPEN"
    | otherwise = "NOT AVAILABLE YET"

showStepInterval :: GeneralConfig -> Double -> String
showStepInterval config step = show $ time_interval config/step

showDepthContour:: GeneralConfig -> String
showDepthContour = show4 . height_limit

showTimeHistoryOutput:: GeneralConfig -> String
showTimeHistoryOutput config
    | output_option config == 1 || output_option config == 2 = "ENABLED"
    | otherwise = "DISABLED"

showZMaxOutput:: GeneralConfig -> String
showZMaxOutput config
    | output_option config == 0 || output_option config == 2 = "ENABLED"
    | otherwise = "DISABLED"