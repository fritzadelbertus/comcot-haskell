module Mass where
    
import TypeModule
    ( updateMiniLayerHZNext,
      LayerConfig(laycord),
      MiniLayer(hnx, hzCurr, hmCurr, hnCurr, h, hny) )
import Constants ( eps, gx, zero )
import Helper (rxLayer, ryLayer)

--------------------------------------------------
-- Solve mass equation for Cartesian

mass:: LayerConfig -> MiniLayer -> MiniLayer
mass layConfig layer 
    -- | isSpherical = massS layer
    | isCartesian = massC layer layConfig
    | otherwise = layer
    where 
        -- isSpherical = laycord (layerConfig layer) == 0
        isCartesian = laycord layConfig == 1


solveConC :: MiniLayer -> LayerConfig -> Int -> Int -> Double
solveConC l c i j = prevFreeSurface - prevXFlux - prevYFlux
        where 
            prevFreeSurface = hzCurr l !! i !! j
            prevXFlux = rxLayer c * ((hmCurr l !! i !! j) - (hmCurr l !! (i-1) !! j))
            prevYFlux = ryLayer c * ((hnCurr l !! i !! j) - (hnCurr l !! i !! (j-1)))

freeSurfaceC:: MiniLayer -> LayerConfig -> Int -> Int -> Double
freeSurfaceC l c i j
    | stillWaterDepth > gx && abs zz < eps = zero
    | stillWaterDepth > gx && dd < eps = - stillWaterDepth
    | stillWaterDepth > gx = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz = solveConC l c i j
        dd = zz + stillWaterDepth

massC:: MiniLayer -> LayerConfig -> MiniLayer
massC layer layConfig = updateMiniLayerHZNext nextLayer layer
    where
        nextLayer = [[freeSurfaceC layer layConfig i j | i <- [0..lasty]] | j <- [0..lastx]]
        lasty = hny layer - 1
        lastx = hnx layer - 1
            
--------------------------------------------------
-- Mass for spherical is not available yet

-- massS:: Layer -> Layer
-- massS layer = updateLayerZ (nextLayer : z layer) layer
--     where
--         nextLayer :: [[Double]]
--         nextLayer = [[eta layer i j | i <- [1..lastx]] | j <- [1..lasty]]
--         lasty = lastYIndex layer
--         lastx = lastXIndex layer
--         eta l i j
--             | stillWaterDepth > gx && abs zz > eps = zz
--             | otherwise = zero
--             where
--                 stillWaterDepth = h l !! i !! j
--                 zz = solveConS l i j


-- solveConS :: MiniLayer -> Int -> Int -> Double
-- solveConS l i j = prevFreeSurface - prevXFlux - prevYFlux
--         where 
--             prevFreeSurface = z l !! i !! j !! 1
--             prevXFlux = (r1 l !! i !! j) * ((m l !! i !! j !! 1) - (m l !! (i-1) !! j !! 1))
--             prevYFlux = (r11 l !! i !! j) * ((n l !! i !! j !! 1) * (r6 l !! j) 
--                 - (n l !! i !! (j-1) !! 1) * (r6 l !! (j-1)))