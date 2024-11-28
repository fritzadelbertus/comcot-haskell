module Initialization where
    
import Helper ( removeLeadingSpaces, isCartesian, generateLayerWithConfig, generateLayer, onResults1, onResults2 )
import Constants (eps, rEarth, radMin, radDeg, grav)
import Data.List (transpose)
import TypeModule
    ( updateLayerConfigDt,
      updateMiniLayerH,
      updateMiniLayerHP,
      updateMiniLayerHQ,
      updateMiniLayerHZCurr,
      Bathymetry(..),
      FaultConfig(..),
      GeneralConfig(..),
      LayerConfig(..),
      MiniLayer(..) )
import Data.Maybe ( fromJust, isNothing )
import Deform ( deformOkada )


--------------------------------------------------
-- Main Functions

readLayerConfig:: String -> IO LayerConfig
readLayerConfig = readConfig getRootLayer

readFaultConfig:: String -> IO FaultConfig
readFaultConfig = readConfig getFaultConfig

readGeneralConfig:: String -> IO GeneralConfig
readGeneralConfig = readConfig getGeneralConfig

getInitialLayer:: LayerConfig -> IO MiniLayer
getInitialLayer layConfig = do
    content <- readFile (depth_name layConfig)
    let layBath = getBathymetryData content
    let initLayer = dxCalc layConfig
    return $ gridInterpolate initLayer layBath

getInitialSurface:: MiniLayer -> LayerConfig -> FaultConfig -> IO MiniLayer
getInitialSurface = getFloorDeform 0.0

adjustBathymetry:: MiniLayer -> MiniLayer
adjustBathymetry = updateHPQ

checkCourantCondition:: MiniLayer -> LayerConfig -> LayerConfig
checkCourantCondition layer layConfig = updateLayerConfigDt newDt layConfig 
    where 
        hMax = (maximum . (map maximum . h)) layer
        courant = dt layConfig/courantCoef
        courantCoef = deltaCourant layer layConfig/(sqrt grav*hMax)
        newDt
            | courant > courantLimit layConfig = courantLimit layConfig * courantCoef
            | otherwise = dt layConfig

--------------------------------------------------
-- Extract Input from '.ctl'

parameterIndent:: Int
parameterIndent = 49

sectionRange:: [[Int]] -- General, Fault Model, Wave Maker, Landslide
sectionRange = [[10..20],[25..40],[45..49],[54..59]]

getAllParameters:: String -> [[String]]
getAllParameters content = map (extractSection content) sectionRange

getChildLayersParameters:: String -> [[String]]
getChildLayersParameters content = map (extractSection content) layersRange
    where layersRange = [[n+2..n+16] | n <- [84,103..331]]

getParameterValue:: String -> Int -> String
getParameterValue content = drop parameterIndent . (lines content !!)


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
        epicenter_y     = read $ params !! 13,
        fault_src_data  = params !! 14
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

readConfig :: (String -> a) -> String -> IO a
readConfig f fileName = do
    content <- readFile fileName
    return $ f content

--------------------------------------------------
-- Extract Bathymetry data

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

--------------------------------------------------
-- Generate Inital Layer

findNxy:: LayerConfig -> (LayerConfig -> Double) -> (LayerConfig -> Double) -> Int
findNxy layConfig start end
    | isCartesian layConfig = cartN 
    | otherwise = sphereN
    where
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
        h = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hp = [[]],
        hq = [[]],
        hzCurr = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hzNext = [[]],
        hmCurr = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hmNext = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hnCurr = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]],
        hnNext = [[0.0 | j <-[0..ny-1]] | i <- [0..nx-1]]
    }
    where
        nx = findNxy layConfig x_start x_end
        ny = findNxy layConfig y_start y_end
        d = dx layConfig

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

gridHeight :: MiniLayer -> Bathymetry -> Int -> Int -> Double
gridHeight miniLayer layBath i j = centerZ + rightZ + downZ + rightDownZ
    where 
        xArr = bx layBath
        yArr = by layBath
        zArr = bz layBath
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

gridInterpolate:: MiniLayer -> Bathymetry -> MiniLayer
gridInterpolate miniLayer layBath = updateMiniLayerH newH miniLayer
    where newH = generateLayerWithConfig miniLayer layBath gridHeight

--------------------------------------------------
-- Extract Initial Deformation from file

getInitSurfData:: String -> Bathymetry
getInitSurfData bathContent = Bathymetry {bx = x, by = y, bz = z, bnx = nx, bny = ny}
    where
        x = generateXCoor bathDataX bathDataY nx ny
        y = generateYCoor bathDataX bathDataY nx ny
        z = generateZCoor bathDataX bathDataY bathDataZ nx ny
        bathDataX = readBathDataX bathContent
        bathDataY = readBathDataY bathContent
        bathDataZ = readBathDataZ bathContent
        (nx,ny) = getBathN (bathDataX,bathDataY) nxy
        nxy = bathymetryFileLength bathContent

getInitialDeformFromData:: FaultConfig -> IO [[Double]]
getInitialDeformFromData fault = do
    content <- readFile (fault_src_data fault)
    let deformBath = getInitSurfData content
    return $ bz deformBath

createNewHZDeform :: MiniLayer -> [[Double]] -> [[Double]]
createNewHZDeform layer deform = [[height i j | j <- [0..hny layer -1]] | i <- [0..hnx layer -1]]
    where 
        height i j = (hzCurr layer !! i !! j) + (deform !! i !! j)


getFloorDeform:: Double -> MiniLayer -> LayerConfig -> FaultConfig -> IO MiniLayer
getFloorDeform time layer layConfig fault = do
    deform <- getInitialDeformFromData fault 
    let newHZ = createNewHZDeform layer deform
    let ret = (updateMiniLayerHZCurr newHZ . updateBathymetry deform) layer
    return ret
    -- where 
        -- deform
            -- | rupture_time fault >= time && rupture_time fault < time + dt layConfig = deformOkada layer layConfig fault
            -- | otherwise = [[0.0 | j <- [0..hny layer -1]] | i <- [0..hnx layer -1]]
        -- newHZ = [[height i j | j <- [0..hny layer -1]] | i <- [0..hnx layer -1]]
        -- height i j = (hzCurr layer !! i !! j) + (deform !! i !! j)

--------------------------------------------------
-- Adjust bathymetry

updateHeightP :: MiniLayer -> Int -> Int -> Double
updateHeightP l i j = 0.5*((h l !! i !! j) + (h l !! ip1 !! j))  
    where 
        ip1
            | i+1 == hnx l = i
            | otherwise = i+1 

updateHeightQ :: MiniLayer -> Int -> Int -> Double
updateHeightQ l i j = 0.5*((h l !! i !! j) + (h l !! i !! jp1))  
    where 
        jp1
            | j+1 == hny l = j
            | otherwise = j+1 

updateHPQ:: MiniLayer -> MiniLayer
updateHPQ layer = (updateMiniLayerHP newHP . updateMiniLayerHQ newHQ) layer
    where 
        newHP = generateLayer layer updateHeightP
        newHQ = generateLayer layer updateHeightQ

updateBathHeight :: MiniLayer -> [[Double]] -> Int -> Int -> Double
updateBathHeight l deform i j = (h l !! i !! j) - (deform !! i !! j)

updateBathymetry:: [[Double]] -> MiniLayer -> MiniLayer
updateBathymetry deform layer = (updateHPQ . updateMiniLayerH newH) layer
    where newH = generateLayerWithConfig layer deform updateBathHeight


--------------------------------------------------
-- Check Courant Condition for Stability

courantLatMax :: MiniLayer -> Double
courantLatMax = (radDeg *) . onResults1 max (abs . (head . hy)) (abs . (last . hy))
courantDiffX :: MiniLayer -> LayerConfig -> Double
courantDiffX layer layConfig
    | isCartesian layConfig = dx layConfig
    | otherwise = rEarth * cos (courantLatMax layer) * dx layConfig * radMin

courantDiffY :: MiniLayer -> LayerConfig -> Double
courantDiffY layer layConfig
    | isCartesian layConfig = dx layConfig
    | otherwise = rEarth * dx layConfig * radMin

courantLimit :: LayerConfig -> Double
courantLimit lConfig 
    | isCartesian lConfig = 0.35
    | otherwise = 0.5

deltaCourant :: MiniLayer -> LayerConfig -> Double
deltaCourant = onResults2 min courantDiffX courantDiffY
        