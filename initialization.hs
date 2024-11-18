module Initialization where
    
import Helper ( removeLeadingSpaces )
import Constants (eps, rEarth, radMin, radDeg, grav)
import Data.List (transpose)
import TypeModule
import Data.Maybe ( fromJust, isNothing )
import Deform

--------------------------------------------------------

sectionRange:: [[Int]] -- General, Fault Model, Wave Maker, Landslide
sectionRange = [[10..20],[25..40],[45..49],[54..59]]

getAllParameters:: String -> [[String]]
getAllParameters content = map (extractSection content) sectionRange

getChildLayersParameters:: String -> [[String]]
getChildLayersParameters content = map (extractSection content) layersRange
    where layersRange = [[n+2..n+16] | n <- [84,103..331]]

-------------------------------------------------------

isCartesian:: LayerConfig -> Bool
isCartesian l = laycord l == 1 

----------------------------
bathymetryFileLength:: String -> Int
bathymetryFileLength content = length (lines content)

readBathData:: ([String]-> String) -> String -> [Double]
readBathData index content = map (read . index . words) (lines content) 

readBathDataX:: String -> [Double]
readBathDataX = readBathData head

readBathDataY:: String -> [Double]
readBathDataY = readBathData (!! 1)

readBathDataZ:: String -> [Double]
readBathDataZ = readBathData (!! 2)

isWrittenRowByRow:: [Double] -> [Double] -> Bool
isWrittenRowByRow x y = xDiff > eps && yDiff < eps
    where 
        xDiff = abs (head x - x !! 1) 
        yDiff = abs (head y - y !! 1)

getBathyNXY:: [Double] -> Int -> Int
getBathyNXY (x:xs) acc
    | head xs > x = getBathyNXY xs acc + 1
    | otherwise = acc

getBathN:: ([Double],[Double]) -> Int -> (Int,Int)
getBathN (x,y) nxy
    | isWrittenRowByRow x y = (getBathyNXY x 1, round nxyfromRow)
    | otherwise = (getBathyNXY y 1, round nxyfromCol)
    where 
        nxyfromRow = fromIntegral nxy/fromIntegral (getBathyNXY x 1)
        nxyfromCol = fromIntegral nxy/fromIntegral (getBathyNXY y 1)

generateXCoor:: [Double] -> [Double] -> Int -> Int -> [Double]
generateXCoor bathDataX bathDataY nx ny
    | isWrittenRowByRow bathDataX bathDataY = [bathDataX !! i | i <-[0..nx-1]]
    | otherwise = [bathDataX !! (i*ny) | i <-[0..nx-1]]

generateYCoor:: [Double] -> [Double] -> Int -> Int -> [Double]
generateYCoor bathDataX bathDataY nx ny
    | not (isWrittenRowByRow bathDataX bathDataY) = [bathDataY !! i | i <-[0..ny-1]]
    | otherwise = [bathDataY !! (i*nx) | i <-[0..ny-1]]

generateZCoor:: [Double] -> [Double] -> [Double] -> Int -> Int -> [[Double]]
generateZCoor bathDataX bathDataY bathDataZ nx ny
    | isWrittenRowByRow bathDataX bathDataY = transpose [take nx (drop (i*nx) bathDataZ) | i <-[0..ny-1]]
    | otherwise = [take ny (drop (i*ny) bathDataZ) | i <-[0..nx-1]]

getBathymetryData:: String -> Bathymetry
getBathymetryData bathContent = Bathymetry {bx = x, by = y, bz = z, bnx = nx, bny = ny}
    where
        x = generateXCoor bathDataX bathDataY nx ny
        y = generateYCoor bathDataX bathDataY nx ny
        z = generateZCoor bathDataX bathDataY bathDataZ nx ny
        bathDataX = readBathDataX bathContent
        bathDataY = readBathDataY bathContent
        bathDataZ = readBathDataZ bathContent
        (nx,ny) = getBathN (bathDataX,bathDataY) nxy
        nxy = bathymetryFileLength bathContent

--------------------------------------------------------------------------

findNxy:: LayerConfig -> (LayerConfig -> Double) -> (LayerConfig -> Double) -> Int
findNxy layConfig start end = n
    where
        n 
            | isCartesian layConfig = cartN
            | otherwise = sphereN
        cartN = 1 + round (distance/d)
        sphereN = 1 + round (distance*60/d)
        distance = end layConfig - start layConfig
        d = dx layConfig

dxCalc:: LayerConfig -> MiniLayer
dxCalc layConfig = MiniLayer { 
        hnx = nx, 
        hny = ny, 
        hx = [x_start layConfig + fromIntegral i*d | i <- [0..nx-1]], 
        hy = [y_start layConfig + fromIntegral j*d | j <- [0..ny-1]], 
        h = [[]],
        hp = [[]],
        hq = [[]],
        hzCurr = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hzNext = [[]],
        hmCurr = [[]],
        hmNext = [[]],
        hnCurr = [[]],
        hnNext = [[]]
    }
    where
        nx = findNxy layConfig x_start x_end
        ny = findNxy layConfig y_start y_end
        d = dx layConfig
---------------------------------

findIndex :: (Ord a) => a -> [a] -> Maybe Int
findIndex value xs = go xs 0
  where
    go (x:y:xs) idx
      | value >= x && value < y = Just idx
      | otherwise = go (y:xs) (idx + 1)
    go _ _ = Nothing

isBetweenInterpolation:: Int -> Int -> Int -> Int -> Bool
isBetweenInterpolation idxI idxJ nx ny = interpolateI && interpolateJ
    where
        interpolateI = idxI >= 0 && idxI < (nx-1) 
        interpolateJ = idxJ >= 0 && idxJ < (ny-1)

gridInterpolate:: MiniLayer -> Bathymetry -> MiniLayer
gridInterpolate miniLayer layBath = updateMiniLayerH newH miniLayer
    where
        newH = [[height i j | j <- [0..endJ]] | i <- [0..endI]]
        endI = hnx miniLayer -1
        endJ = hny miniLayer -1
        nx = bnx layBath
        ny = bny layBath
        xArr = bx layBath
        yArr = by layBath
        zArr = bz layBath
        height i j = centerZ + rightZ + downZ + rightDownZ

            where 
                indexI = findIndex (hx miniLayer !! i) xArr
                indexJ = findIndex (hy miniLayer !! j) yArr
                idxI
                    | isNothing indexI = hnx miniLayer -2
                    | otherwise = fromJust indexI
                idxJ
                    | isNothing indexJ = hny miniLayer -2
                    | otherwise = fromJust indexJ
                centerZ = zArr !! idxI !! idxJ * (1.0-cx)*(1.0-cy)
                rightZ = zArr !! (idxI+1) !! idxJ * cx * (1.0-cy)
                downZ = zArr !! idxI !! (idxJ+1) * (1.0-cx) * cy
                rightDownZ = zArr !! (idxI+1) !! (idxJ+1) * cx * cy

                cx = ((hx miniLayer !! i) - xArr !! idxI) / deltaX
                cy = ((hy miniLayer !! j) - yArr !! idxJ) / deltaY
                deltaX = xArr !! (idxI+1) - xArr !! idxI
                deltaY = yArr !! (idxJ+1) - yArr !! idxJ

 ----------------------------------------

parameterIndent:: Int
parameterIndent = 49

getParameterValue:: String -> Int -> String
getParameterValue content line = drop parameterIndent $ lines content !! line


extractSection:: String -> [Int] -> [String]
extractSection content = map (removeLeadingSpaces . getParameterValue content)

getRootLayerParameters:: String -> [String]
getRootLayerParameters content = extractSection content rootLayerRange
    where rootLayerRange = [66..82]

getRootLayer:: String -> LayerConfig
getRootLayer content = LayerConfig {
        layswitch   = read $ head params,
        laycord     = read $ params !! 1,
        laygov      = read $ params !! 2,
        dx          = read $ params !! 3,
        dt          = read $ params !! 4,
        fric_switch = read $ params !! 5,
        fric_coef   = read $ params !! 6,
        fluxswitch  = read $ params !! 7,
        x_start     = read $ params !! 8,
        x_end       = read $ params !! 9,
        y_start     = read $ params !! 10,
        y_end       = read $ params !! 11,
        depth_name  = params !! 12,
        fs          = read $ params !! 13,
        layer_id    = read $ params !! 14,
        level       = read $ params !! 15,
        parent      = read $ params !! 16
    }
    where params = getRootLayerParameters content

getFaultParameters:: String -> [String]
getFaultParameters content = extractSection content faultLayerRange
    where faultLayerRange = [25..40]

getFaultConfig:: String -> FaultConfig
getFaultConfig content = FaultConfig {
        num_fault       = read $ head params,
        rupture_time    = read $ params !! 1,
        switch          = read $ params !! 2,
        focal_depth     = read $ params !! 3,    
        fault_length    = read $ params !! 4,
        fault_width     = read $ params !! 5,
        dislocation     = read $ params !! 6,
        strike_angle    = read $ params !! 7,
        dip_angle       = read $ params !! 8,
        slip_angle      = read $ params !! 9,
        comp_origin_x   = read $ params !! 10,
        comp_origin_y   = read $ params !! 11,
        epicenter_x     = read $ params !! 12,
        epicenter_y     = read $ params !! 13
    }
    where params = getFaultParameters content

getGeneralParameters:: String -> [String]
getGeneralParameters content = extractSection content generalRange
    where generalRange = [10..20]

getGeneralConfig:: String -> GeneralConfig
getGeneralConfig content = GeneralConfig {
        total_time          = read $ head params,
        time_interval       = read $ params !! 1,
        output_option       = read $ params !! 2,
        start_type          = read $ params !! 3,
        start_time          = read $ params !! 4,
        height_limit        = read $ params !! 5,
        initial_surface     = read $ params !! 6,
        boundary_condition  = read $ params !! 7
    }
    where params = getGeneralParameters content

readLayerConfig:: String -> IO LayerConfig
readLayerConfig fileName = do
    content <- readFile fileName
    return $ getRootLayer content

readFaultConfig:: String -> IO FaultConfig
readFaultConfig fileName = do
    content <- readFile fileName
    return $ getFaultConfig content

readGeneralConfig:: String -> IO GeneralConfig
readGeneralConfig fileName = do
    content <- readFile fileName
    return $ getGeneralConfig content

-----------------------------------------
getInitialLayer:: LayerConfig -> IO MiniLayer
getInitialLayer layConfig = do
    content <- readFile (depth_name layConfig)
    let layBath = getBathymetryData content
    let initLayer = dxCalc layConfig
    return $ gridInterpolate initLayer layBath

----------------------------------------------

getFloorDeform:: Double -> MiniLayer -> LayerConfig -> FaultConfig -> MiniLayer
getFloorDeform time layer layConfig fault = (updateMiniLayerHZCurr newHZ . updateBathymetry deform) layer
    where 
        deform
            | rupture_time fault >= time && rupture_time fault < time + dt layConfig = deformOkada layer layConfig fault
            | otherwise = [[0.0 | j <- [0..hny layer -1]] | i <- [0..hnx layer -1]]
        newHZ = [[height i j | j <- [0..hny layer -1]] | i <- [0..hnx layer -1]]
        height i j = (hzCurr layer !! i !! j) + (deform !! i !! j)

getInitialSurface:: MiniLayer -> LayerConfig -> FaultConfig -> MiniLayer
getInitialSurface = getFloorDeform 0.0

----------------------------------------------

updateHPQ:: MiniLayer -> MiniLayer
updateHPQ layer = (updateMiniLayerHP newHP . updateMiniLayerHQ newHQ) layer
    where 
        newHP = [[heightP layer i j | j <- [0..hny layer-1]] | i <- [0..hnx layer-1]]
        newHQ = [[heightQ layer i j | j <- [0..hny layer-1]] | i <- [0..hnx layer-1]]
        heightP l i j = 0.5*((h l !! i !! j) + (h l !! ip1 !! j))  
            where 
                ip1
                    | i+1 == hnx l = i
                    | otherwise = i+1 
        heightQ l i j = 0.5*((h l !! i !! j) + (h l !! i !! jp1))  
            where 
                jp1
                    | j+1 == hny l = j
                    | otherwise = j+1 

adjustBathymetry:: MiniLayer -> MiniLayer
adjustBathymetry = updateHPQ

updateBathymetry:: [[Double]] -> MiniLayer -> MiniLayer
updateBathymetry deform layer = (updateHPQ . updateMiniLayerH newH) layer
    where 
        newH = [[height i j | j <- [0..endJ]] | i <- [0..endI]]
        endI = hnx layer -1
        endJ = hny layer -1
        height i j = (h layer !! i !! j) - (deform !! i !! j)

------------------------------------------


checkCourantCondition:: MiniLayer -> LayerConfig -> LayerConfig
checkCourantCondition layer layConfig = updateLayerConfigDt newDt layConfig 
    where 
        hMax = maximum $ map maximum (h layer)
        latMax = max (abs (head (hy layer))) (abs (last (hy layer))) * radDeg
        diffX
            | isCartesian layConfig = dx layConfig
            | otherwise = rEarth * cos latMax * dx layConfig * radMin
        diffY
            | isCartesian layConfig = dx layConfig
            | otherwise = rEarth * dx layConfig * radMin
        delta = min diffX diffY
        courant = dt layConfig/courantCoef
        courantCoef = delta/(sqrt grav*hMax)
        courantLimit
            | isCartesian layConfig = 0.35
            | otherwise = 0.5
        newDt
            | courant > courantLimit = courantLimit * courantCoef
            | otherwise = dt layConfig
        