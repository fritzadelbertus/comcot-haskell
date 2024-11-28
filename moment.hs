module Moment where

import TypeModule
    ( updateMiniLayerHMNext,
      updateMiniLayerHNNext,
      LayerConfig(laycord),
      MiniLayer(h, hp, hmCurr, hq, hzNext, hnCurr, hny, hnx) )
import Constants ( eps, gx, zero )
import Helper
    ( generateLayerWithConfig, grxLayer, gryLayer, onResults3 )

--------------------------------------------------
-- Solve moment equation for Cartesian

moment:: LayerConfig -> MiniLayer -> MiniLayer
moment layConfig layer
    -- | isSpherical = momtS layer
    | isCartesian = momtC layer layConfig
    | otherwise = layer
    where 
        isSpherical = laycord layConfig == 0
        isCartesian = laycord layConfig == 1

heightFactorMC:: LayerConfig -> MiniLayer ->  Int -> Int -> Double
heightFactorMC c l i j = grxLayer c * hm * ((hzNext l !! (i+1) !! j) - (hzNext l !! i !! j)) 
    where hm = (hp l !! i !! j) + 0.5 * ((hzNext l !! i !! j) + (hzNext l !! (i+1) !! j))

prevMomentMC:: MiniLayer -> Int -> Int -> Double
prevMomentMC =  ((!!) .) . ((!!) . hmCurr)

solveMomMC:: LayerConfig -> MiniLayer -> Int -> Int -> Double
solveMomMC c = onResults3 (-) prevMomentMC (heightFactorMC c)

heightFactorNC:: LayerConfig -> MiniLayer ->  Int -> Int -> Double
heightFactorNC c l i j = gryLayer c * hn * ((hzNext l !! i !! (j + 1))-(hzNext l !! i !! j))
    where hn = (hq l !! i !! j) + 0.5 * ((hzNext l !! i !! j) + (hzNext l !! i !! (j+1)))

prevMomentNC:: MiniLayer -> Int -> Int -> Double
prevMomentNC =  ((!!) .) . ((!!) . hnCurr)

solveMomNC:: LayerConfig -> MiniLayer -> Int -> Int -> Double
solveMomNC c = onResults3 (-) prevMomentNC (heightFactorNC c)
        

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
        xm = solveMomMC c l i j

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
        xn = solveMomNC c l i j

momtC:: MiniLayer -> LayerConfig -> MiniLayer
momtC layer layConfig = (updateMiniLayerHNNext nextLayerN . updateMiniLayerHMNext nextLayerM) layer
    where
        nextLayerM = generateLayerWithConfig layer layConfig fluxMC
        nextLayerN = generateLayerWithConfig layer layConfig fluxNC

--------------------------------------------------
-- Moment for spherical is not available yet

-- startMoment:: Layer -> Int
-- startMoment l
--     | layer_id (layerConfig l) == 1 = 0
--     | otherwise = 1

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