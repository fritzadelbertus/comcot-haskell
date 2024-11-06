import TypeModule
import Constants

moment:: Layer -> Layer
moment layer 
    | isSpherical = momtS layer
    | isCartesian = momtC layer
    where 
        isSpherical = laycord (layerConfig layer) == 0
        isCartesian = laycord (layerConfig layer) == 1

startMoment:: Layer -> Int
startMoment l
    | layer_id (layerConfig l) == 1 = 0
    | otherwise = 1

solveMomMS:: Layer -> Int -> Int -> Double
solveMomMS l i j = prevMomentM - nextHeight + prevMomentN 
    where 
        prevMomentM = m l !! i !! j !! 1
        nextHeight = (r2 l !! i !! j) * ((z l !! (i+1) !! j !! 2) - (z l !! i !! j !! 2))
        prevMomentN = (r3 l !! i !! j) * ((n l !! i !! j !! 1) + (n l !! (i+1) !! j !! 1)
                + (n l !! i !! prevJ !! 1) + (n l !! (i+1) !! prevJ !! 1))
        prevJ
            | j-1 < 1 = 1
            | otherwise = j-1

solveMomNS:: Layer -> Int -> Int -> Double
solveMomNS l i j = prevMomentN - nextHeight - prevMomentM
    where 
        prevMomentN = n l !! i !! j !! 1
        nextHeight = (r4 l !! i !! j)*((z l !! i !! (j+1) !! 2) - (z l !! i !! j !! 2))
        prevMomentM = (r5 l !! i !! j)* ((m l !! prevI !! j !! 1) + (m l !! prevI !! (j+1) !! 1) 
                    + (m l !! i !! j !! 1) + (m l !! i !! (j+1) !! 1))
        prevI
            | i-1 < 1 = 1
            | otherwise = i-1


momtS:: Layer -> Layer
momtS layer = (updateLayerN (nextLayerN : n layer) . updateLayerM (nextLayerM : m layer)) layer
    where
        nextLayerM :: [[Double]]
        nextLayerM = [[fluxM layer i j | i <- [start..lastMx]] | j <- [start..lastMy]]
        lastMy = ny layer + 1
        lastMx = nx layer
        start = startMoment layer
        fluxM l i j
            | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xm > eps = xm
            | otherwise = zero
            where
                stillWaterDepth1 = h l !! i !! j
                stillWaterDepth2 = h l !! (i + 1) !! j
                xm = solveMomMS l i j
        nextLayerN :: [[Double]]
        nextLayerN = [[fluxN layer i j | i <- [start..lastNx]] | j <- [start..lastNy]]
        lastNy = ny layer
        lastNx = nx layer + 1
        fluxN l i j
            | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xn > eps = xn
            | otherwise = zero
            where
                stillWaterDepth1 = h l !! i !! j
                stillWaterDepth2 = h l !! i !! (j + 1)
                xn = solveMomNS l i j

solveMomMC:: Layer -> Int -> Int -> Double
solveMomMC l i j = prevMomentM - heightFactor
    where 
        prevMomentM = m l !! i !! j !! 1
        heightFactor = grx l * hm * ((z l !! (i+1) !! j !! 2) - (z l !! i !! j !! 2)) 
        hm = (hp l !! i !! j) + 0.5 * ((z l !! i !! j !! 2) + (z l !! (i+1) !! j !! 2))

solveMomNC:: Layer -> Int -> Int -> Double
solveMomNC l i j = prevMomentN - heightFactor
    where 
        prevMomentN = n l !! i !! j !! 1
        heightFactor = gry l * hn * ((z l !! i !! (j + 1) !! 2)-(z l !! i !! j !! 2))
        hn = (hq l !! i !! j) + 0.5 * ((z l !! i !! j !! 2) + (z l !! i !! (j+1) !! 2))

momtC:: Layer -> Layer
momtC layer = (updateLayerN (nextLayerN : n layer) . updateLayerM (nextLayerM : m layer)) layer
    where
        nextLayerM :: [[Double]]
        nextLayerM = [[fluxM layer i j | i <- [start..lastMx]] | j <- [start..lastMy]]
        lastMy = ny layer + 1
        lastMx = nx layer
        start = startMoment layer
        fluxM l i j
            | stillWaterDepth1 && stillWaterDepth2  && abs xm > eps = xm
            | otherwise = zero
            where
                stillWaterDepth1 = h l !! i !! j > gx
                stillWaterDepth2 = h l !! (i + 1) !! j > gx
                xm = solveMomMC l i j
        nextLayerN :: [[Double]]
        nextLayerN = [[fluxN layer i j | i <- [start..lastNx]] | j <- [start..lastNy]]
        lastNy = ny layer
        lastNx = nx layer + 1
        fluxN l i j
            | stillWaterDepth1 > gx && stillWaterDepth2 > gx && abs xn > eps = xn
            | otherwise = zero
            where
                stillWaterDepth1 = h l !! i !! j
                stillWaterDepth2 = h l !! i !! (j + 1)
                xn = solveMomNC l i j
