module Moment where
import TypeModule
import Constants
import Helper
import Control.Applicative (Applicative(liftA2))
moment:: LayerConfig -> MiniLayer -> MiniLayer
moment layConfig layer
    -- | isSpherical = momtS layer
    | isCartesian = momtC layer layConfig
    | otherwise = layer
    where 
        isSpherical = laycord layConfig == 0
        isCartesian = laycord layConfig == 1

startMoment:: Layer -> Int
startMoment l
    | layer_id (layerConfig l) == 1 = 0
    | otherwise = 1

-- solveMomMS:: Layer -> Int -> Int -> Double
-- solveMomMS l i j = prevMomentM - nextHeight + prevMomentN 
--     where 
--         prevMomentM = m l !! i !! j !! 1
--         nextHeight = (r2 l !! i !! j) * ((z l !! (i+1) !! j !! 2) - (z l !! i !! j !! 2))
--         prevMomentN = (r3 l !! i !! j) * ((n l !! i !! j !! 1) + (n l !! (i+1) !! j !! 1)
--                 + (n l !! i !! prevJ !! 1) + (n l !! (i+1) !! prevJ !! 1))
--         prevJ
--             | j-1 < 1 = 1
--             | otherwise = j-1

-- solveMomNS:: Layer -> Int -> Int -> Double
-- solveMomNS l i j = prevMomentN - nextHeight - prevMomentM
--     where 
--         prevMomentN = n l !! i !! j !! 1
--         nextHeight = (r4 l !! i !! j)*((z l !! i !! (j+1) !! 2) - (z l !! i !! j !! 2))
--         prevMomentM = (r5 l !! i !! j)* ((m l !! prevI !! j !! 1) + (m l !! prevI !! (j+1) !! 1) 
--                     + (m l !! i !! j !! 1) + (m l !! i !! (j+1) !! 1))
--         prevI
--             | i-1 < 1 = 1
--             | otherwise = i-1


-- momtS:: Layer -> Layer
-- momtS layer = (updateLayerN (nextLayerN : n layer) . updateLayerM (nextLayerM : m layer)) layer
--     where
--         nextLayerM :: [[Double]]
--         nextLayerM = [[fluxM layer i j | i <- [start..lastMx]] | j <- [start..lastMy]]
--         lastMy = ny layer + 1
--         lastMx = nx layer
--         start = startMoment layer
--         fluxM l i j
--             | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xm > eps = xm
--             | otherwise = zero
--             where
--                 stillWaterDepth1 = h l !! i !! j
--                 stillWaterDepth2 = h l !! (i + 1) !! j
--                 xm = solveMomMS l i j
--         nextLayerN :: [[Double]]
--         nextLayerN = [[fluxN layer i j | i <- [start..lastNx]] | j <- [start..lastNy]]
--         lastNy = ny layer
--         lastNx = nx layer + 1
--         fluxN l i j
--             | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xn > eps = xn
--             | otherwise = zero
--             where
--                 stillWaterDepth1 = h l !! i !! j
--                 stillWaterDepth2 = h l !! i !! (j + 1)
--                 xn = solveMomNS l i j

solveMomMC:: MiniLayer -> LayerConfig -> Int -> Int -> Double
solveMomMC l c i j = prevMomentM - heightFactor
    where 
        prevMomentM = hmCurr l !! i !! j
        heightFactor = grxLayer c * hm * ((hzNext l !! (i+1) !! j) - (hzNext l !! i !! j)) 
        hm = (hp l !! i !! j) + 0.5 * ((hzNext l !! i !! j) + (hzNext l !! (i+1) !! j))

heightFactorNC:: MiniLayer -> LayerConfig -> Int -> Int -> Double
heightFactorNC l c i j = gryLayer c * hn * ((hzNext l !! i !! (j + 1))-(hzNext l !! i !! j))
    where hn = (hq l !! i !! j) + 0.5 * ((hzNext l !! i !! j) + (hzNext l !! i !! (j+1)))

prevMomentNC:: MiniLayer -> LayerConfig -> Int -> Int -> Double
prevMomentNC l c i j = hnCurr l !! i !! j


solveMomNC:: MiniLayer -> LayerConfig -> Int -> Int -> Double
solveMomNC = onResults (-) prevMomentNC heightFactorNC
        

fluxMC :: MiniLayer -> LayerConfig -> Int -> Int -> Double
fluxMC l c i j
    | i > lastMx = zero
    | stillWaterDepth1 && stillWaterDepth2  && abs xm > eps = xm
    | otherwise = zero
    where
        lastMy = hny l - 1
        lastMx = hnx l - 2
        stillWaterDepth1 = h l !! i !! j > gx
        stillWaterDepth2 = h l !! (i + 1) !! j > gx
        xm = solveMomMC l c i j

fluxNC :: MiniLayer -> LayerConfig -> Int -> Int -> Double
fluxNC l c i j
    | j > lastNy = zero
    | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xn > eps = xn
    | otherwise = zero
    where
        lastNy = hny l - 2
        lastNx = hnx l - 1
        stillWaterDepth1 = h l !! i !! j
        stillWaterDepth2 = h l !! i !! (j + 1)
        xn = solveMomNC l c i j

momtC:: MiniLayer -> LayerConfig -> MiniLayer
momtC layer layConfig = (updateMiniLayerHNNext nextLayerN . updateMiniLayerHMNext nextLayerM) layer
    where
        start = 0
        nextLayerM :: [[Double]]
        nextLayerM = generateLayerWithConfig layer layConfig fluxMC
        nextLayerN :: [[Double]]
        nextLayerN = generateLayerWithConfig layer layConfig fluxNC
