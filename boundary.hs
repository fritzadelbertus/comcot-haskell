module Boundary where
import TypeModule
    ( updateMiniLayerHZNext,
      MiniLayer(hzNext, hnNext, hmNext, h, hnx, hny) )
import Constants ( grav, gx, ub, zero )

--------------------------------------------------
-- Formula for Open Boundary

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

findUpper:: MiniLayer -> Int -> Int -> Double
findUpper l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz
            | hmNext l !! i !! j > zero = -(uu/cc)
            | otherwise = uu/cc
        uu = sqrt(uh**2+(hmNext l !! i !! j)**2)
        uh
            | h l !! i !! (j-1) > gx = 0.5 * ((hnNext l !! i !! j)+(hnNext l !! i !! (j-1)))
            | otherwise = hnNext l !! i !! j
        cc = sqrt(grav*(h l !! i !! j))

findBottom:: MiniLayer -> Int -> Int -> Double
findBottom l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz
            | hmNext l !! (i-1) !! j < zero = -(uu/cc)
            | otherwise = uu/cc
        uu = sqrt(uh**2+(hmNext l !! (i-1) !! j)**2)
        uh = 0.5 * ((hnNext l !! i !! j)+(hnNext l !! i !! (j-1)))
        cc = sqrt(grav*(h l !! i !! j))

findLeft:: MiniLayer -> Int -> Int -> Double
findLeft l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz
            | hnNext l !! i !! j > zero = -(uu/cc)
            | otherwise = uu/cc
        uu = sqrt(uh**2+(hnNext l !! i !! j)**2)
        uh = 0.5 * ((hmNext l !! i !! j)+(hmNext l !! (i-1) !! j))
        cc = sqrt(grav*(h l !! i !! j))

findRight:: MiniLayer -> Int -> Int -> Double
findRight l i j
    | stillWaterDepth > gx && abs zz <= ub = zz
    | otherwise = zero
    where
        stillWaterDepth = h l !! i !! j
        zz
            | hnNext l !! i !! (j-1) < zero = -(uu/cc)
            | otherwise = uu/cc
        uu = sqrt(uh**2+(hnNext l !! i !! (j-1))**2)
        uh = 0.5 * ((hmNext l !! i !! j)+(hmNext l !! (i-1) !! j))
        cc = sqrt(grav*(h l !! i !! j))

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

--------------------------------------------------
-- Other formula not yet implemented