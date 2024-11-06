import TypeModule
import Constants

mass:: Layer -> Layer
mass layer 
    | isSpherical = massS layer
    | isCartesian = massC layer
    where 
        isSpherical = laycord (layerConfig layer) == 0
        isCartesian = laycord (layerConfig layer) == 1

solveConS :: Layer -> Int -> Int -> Double
solveConS l i j = prevFreeSurface - prevXFlux - prevYFlux
        where 
            prevFreeSurface = z l !! i !! j !! 1
            prevXFlux = (r1 l !! i !! j) * ((m l !! i !! j !! 1) - (m l !! (i-1) !! j !! 1))
            prevYFlux = (r11 l !! i !! j) * ((n l !! i !! j !! 1) * (r6 l !! j) 
                - (n l !! i !! (j-1) !! 1) * (r6 l !! (j-1)))

solveConC :: Layer -> Int -> Int -> Double
solveConC l i j = prevFreeSurface - prevXFlux - prevYFlux
        where 
            prevFreeSurface = z l !! i !! j !! 1
            prevXFlux = rx l * ((m l !! i !! j !! 1) - (m l !! (i-1) !! j !! 1))
            prevYFlux = ry l * ((n l !! i !! j !! 1) - (n l !! i !! (j-1) !! 1))

lastYIndex :: Layer -> Int
lastYIndex l
    | layer_id (layerConfig l) == 1 = ny l
    | otherwise = ny l + 1

lastXIndex :: Layer -> Int
lastXIndex l
    | layer_id (layerConfig l) == 1 = nx l
    | otherwise = nx l + 1

massS:: Layer -> Layer
massS layer = updateLayerZ (nextLayer : z layer) layer
    where
        nextLayer :: [[Double]]
        nextLayer = [[eta layer i j | i <- [1..lastx]] | j <- [1..lasty]]
        lasty = lastYIndex layer
        lastx = lastXIndex layer
        eta l i j
            | stillWaterDepth > gx && abs zz > eps = zz
            | otherwise = zero
            where
                stillWaterDepth = h l !! i !! j
                zz = solveConS l i j

massC:: Layer -> Layer
massC layer = updateLayerZ (nextLayer : z layer) layer
    where
        nextLayer :: [[Double]]
        nextLayer = [[eta layer i j | i <- [1..lastx]] | j <- [1..lasty]]
        lasty = lastYIndex layer
        lastx = lastXIndex layer
        eta l i j
            | stillWaterDepth > gx && abs zz < eps = zero
            | stillWaterDepth > gx && dd < eps = - stillWaterDepth
            | stillWaterDepth > gx = zz
            | otherwise = zero
            where
                stillWaterDepth = h l !! i !! j
                zz = solveConC l i j
                dd = zz + stillWaterDepth
