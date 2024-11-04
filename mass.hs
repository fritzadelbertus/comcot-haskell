import TypeModule
import Constants

updateLayerZ:: [[[Double]]] -> Layer -> Layer
updateLayerZ newZ layer = layer {z = newZ}

updateLayerDz:: [[[Double]]] -> Layer -> Layer
updateLayerDz newDz layer = layer {dz = newDz}

mass:: Layer -> Layer
mass layer 
    | isSpherical && isLinear       = massS layer
    | isSpherical && isNonlinear    = conmassS layer
    | isCartesian && isLinear       = massC layer
    | isCartesian && isNonlinear    = conmassC layer
    where 
        isSpherical = laycord (layerConfig layer) == 0
        isCartesian = laycord (layerConfig layer) == 1
        isLinear    = laygov (layerConfig layer) == 0
        isNonlinear = laygov (layerConfig layer) == 1

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

conmass:: (Layer -> Int -> Int -> Double) -> Layer -> Layer
conmass continuity layer = (updateLayerDz (nextLayerDz : dz layer) . updateLayerZ (nextLayerZ : z layer)) layer
    where
        lasty = lastYIndex layer
        lastx = lastXIndex layer
        nextLayerZ::[[Double]]
        nextLayerZ = [[eta layer i j | i <- [1..lastx]] | j <- [1..lasty]]
        eta l i j
            | stillWaterDepth > elmax && abs zz < eps = zero
            | stillWaterDepth > elmax && dd > gx = zz
            | otherwise = - stillWaterDepth
            where
                stillWaterDepth = h l !! i !! j
                zz = continuity l i j
                dd = zz + stillWaterDepth
        nextLayerDz::[[Double]]
        nextLayerDz = [[eta2 layer i j | i <- [1..lastx]] | j <- [1..lasty]]
        eta2 l i j
            | stillWaterDepth > elmax && dd > gx = dd
            | otherwise = zero
            where
                stillWaterDepth = h l !! i !! j
                zz = continuity l i j
                dd = zz + stillWaterDepth

conmassC:: Layer -> Layer
conmassC = conmass solveConC

conmassS:: Layer -> Layer
conmassS = conmass solveConS