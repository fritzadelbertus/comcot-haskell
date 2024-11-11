import TypeModule
import Constants

findUpperLeft:: MiniLayer -> Int -> Int -> Double
findUpperLeft l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz = uh/cc
        uh = sqrt(qx**2+qy**2)
        cc = sqrt(grav*(h l !! i !! j))
        qx = hmNext l !! i !! j
        qy = hnNext l !! i !! j

findBottomLeft:: MiniLayer -> Int -> Int -> Double
findBottomLeft l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz = uh/cc
        uh = sqrt(qx**2+qy**2)
        cc = sqrt(grav*(h l !! i !! j))
        qx = hmNext l !! (i-1) !! j
        qy = hnNext l !! i !! j

findUpperRight:: MiniLayer -> Int -> Int -> Double
findUpperRight l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz = uh/cc
        uh = sqrt(qx**2+qy**2)
        cc = sqrt(grav*(h l !! i !! j))
        qx = hmNext l !! i !! j
        qy = hnNext l !! i !! (j-1)

findBottomRight:: MiniLayer -> Int -> Int -> Double
findBottomRight l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz = uh/cc
        uh = sqrt(qx**2+qy**2)
        cc = sqrt(grav*(h l !! i !! j))
        qx = hmNext l !! (i-1) !! j
        qy = hnNext l !! i !! (j-1)

-- findUpper:: MiniLayer -> Int -> Int -> Double
-- findUpper l i j
--     | stillWaterDepth > gx && abs zz <= ub = zz
--     | otherwise = zero
--     where
--         stillWaterDepth = h l !! i !! j
--         zz = uh/cc
--         uh = sqrt(qx**2+qy**2)
--         cc = sqrt(grav*(h l !! i !! j))
--         qx = hmNext l !! (i-1) !! j
--         qy = hnNext l !! i !! (j-1)

open:: MiniLayer -> MiniLayer
open layer = updateMiniLayerHZNext newHz layer
    where 
        newHz :: [[Double]]
        newHz = [[eta layer i j | j <- [0..hny layer - 1]] | i <- [0..hnx layer - 1]]
        eta l i j
            | upBound && leftBound = findUpperLeft l i j
            | upBound && rightBound = findUpperRight l i j
            | botBound && leftBound = findBottomLeft l i j
            | botBound && rightBound = findBottomRight l i j
            | upBound = findUpper l i j
            | botBound = findBottom l i j
            | rightBound = findRight l i j
            | leftBound = findLeft l i j
            | otherwise = hzNext l !! i !! j 
            where 
                upBound = i == 0
                botBound = i == (hnx l - 1)
                rightBound = j == (hny layer - 1)
                leftBound = j == 0