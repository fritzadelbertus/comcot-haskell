module Initialization where

import Helper (removeLeadingSpaces)
import TypeModule

-- Constants
parameterIndent:: Int
parameterIndent = 49

getParameterValue:: String -> Int -> String
getParameterValue content line = drop parameterIndent $ lines content !! line

extractSection:: String -> [Int] -> [String]
extractSection content = map (removeLeadingSpaces . getParameterValue content)

--------------------------------------------------------
-- 1. READ CONFIG (FILE READING)

sectionRange:: [[Int]] -- General, Fault Model, Wave Maker, Landslide
sectionRange = [[10..20],[25..40],[45..49],[54..59]]

getAllParameters:: String -> [[String]]
getAllParameters content = map (extractSection content) sectionRange

getFaultParameters:: String -> [String]
getFaultParameters content = extractSection content faultLayerRange
    where faultLayerRange = [25..40]

getRootLayerParameters:: String -> [String]
getRootLayerParameters content = extractSection content rootLayerRange
    where rootLayerRange = [66..82]

getChildLayersParameters:: String -> [[String]]
getChildLayersParameters content = map (extractSection content) layersRange
    where layersRange = [[n+2..n+17] | n <- [85,105..345]]

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

-- getFaultConfig:: String -> FaultConfig
-- getFaultConfig content = FaultConfig {
--         num_flt     :: Int,    
--         hh          :: Double,    
--         l           :: Double,
--         w           :: Double,
--         d           :: Double,
--         th          :: Double,
--         dl          :: Double,
--         rd          :: Double,
--         fault_xo    :: Double,
--         fault_yo    :: Double,
--         x0          :: Double,
--         y0          :: Double,
--         t0          :: Double,
--         switch      :: Int,
        
--         fault_fs    :: Int,
--         deform_name :: String
--     }
    -- where params = getRootLayerParameters content
readComcotCtl:: IO()
readComcotCtl = do
    content <- readFile "comcot.ctl"
    print $ getAllParameters content
    print $ getRootLayer content
    print $ getChildLayersParameters content

-- 2. GET_INI_SURF  (LOGIC)
-- !DESCRIPTION:
-- !	  #. OBTAIN INITIAL FREE SURFACE DISPLACEMENT FROM FAULT MODEL,
-- !	     CUSTOMIZED DATA FILE OR LANDSLIDE MODEL;
-- !	  #. INTERPOLATE FREE SURFACE DISPLACEMENT INTO NESTED GRID LAYERS;
-- !	  #. INI_SURF =
-- !			0: USE OKADA'S FAULT MODEL TO CALCULATE DEFORMATION
-- !			1: USE AN EXTERNAL FILE TO DETERMINE INITIAL SURFACE
-- !			2: USE INCIDENT WAVE MODEL TO DETERMINE INITIAL SURFACE
-- !			3: USE SUBMARINE LANDSLIDE MODEL
-- !			4: USE MULTIPLE FAULTS + LANDSLIDE (REQUIRE FAULT_MULTI.CTL)
-- !			9: USE MANSINHA AND SMYLIES' MODEL TO CALCULATE DEFORMATION
-- !	  #. FAULT MODELS ARE CALLED IN THIS SUBROUTINE
-- !INPUT:
-- !	  #. GRID INFORMATION, FAULT PARAMETERS
-- !OUTPUT:
-- !	  #. INITIAL WATER SURFACE DISPLACEMENTS OF ALL GRID LAYERS
-- !     #. INITIAL SURFACE DISPLACEMENT IS SAVED IN INI_SURFACE.DAT
-- !     #. SEAFLOOR DISPLACEMENTS ARE SAVED IN DEFORM_SEGXX.DAT

-- main :: IO ()
-- main = do

-- 3. GET_MULTIFAULT_PARAMETERS
-- !DESCRIPTION:
-- !	  #. OBTAIN FAULT PARAMETERS FOR SINGLE/MULTIPLE FAULT;
-- !	  #. FAULT_MULTI.CTL IS REQUIRED IF TOTAL NUMBER OF FAULT PLANES
-- !		 IS LARGER THAN 1
-- !INPUT:
-- !	  #. COMCOT.CTL AND FAULT_MULTI.CTL IF REQUIRED;
-- !OUTPUT:
-- !	  #. FAULT PARAMETERS FOR ALL FAULT PLANES INCLUDED;

faultsRange:: [[Int]]
faultsRange = [[n+3..n+17] | n <- [10,29..732]]

getAllFaultModelParameter:: String -> [[String]]
getAllFaultModelParameter content = map (extractSection content) faultsRange

readFaultMultiCtl :: IO ()
readFaultMultiCtl = do
    content <- readFile "fault_multi.ctl"
    print $ getAllFaultModelParameter content

-- 4. READ_MULTIFAULT_DATA (LO,FLT)
-- !DESCRIPTION:
-- !	  #. OBTAIN FAULT PARAMETERS FOR MULTIPLE FAULT SEGMENTS;
-- !	  #. PARAMETERS FOR FAULT SEGMENTS ARE OBTAINED FROM AN EXTERNAL
-- !	     DATA FILE GIVEN AT LINE 40 IN COMCOT.CTL;
-- !	  #. THE DATA FILE CONTAINS ALL THE INFORMATION FOR EACH SEGMENT;
-- !		 EACH ROW FOR ONE SEGMENT: TIME,LON,LAT,L,W,H,THETA,DELTA,LAMBDA,SLIP;
-- !	  #. TO USE THIS FUNCTION, INPUT 999 AT LINE 26 IN COMCOT.CTL;
-- !INPUT:
-- !	  #. PARAMETER DATA FILE;
-- !OUTPUT:
-- !	  #. FAULT PARAMETERS FOR ALL FAULT PLANES INCLUDED;

-- main :: IO ()
-- main = do
--     content <- readFile "comcot.ctl"
--     print $ getAllParameters content
--     print $ getRootLayerParameters content
--     print $ layersRange
--     print $ getChildLayersParameters content

-- main :: IO ()
-- main = do

-- 5. GET_LANDSLIDE_PARAMETERS
-- !DESCRIPTION:
-- !	  #. OBTAIN ADDITIONAL PARAMETERS FOR LANDSLIDE CONFIGURATION;
-- !	  #. THESE ADDITIONAL PARAMETERS ARE USED TO DETERMINE WATER DEPTH
-- !	     VARIATIONS VIA WATTS ET AL (2003)'S LANDSLIDE THEORY;
-- !	  #. LANDSLIDE.CTL IS REQUIRED IF THE OPTION IN LANDSLIDE SECTION
-- !		 IN COMCOT.CTL IS LARGER THAN 1
-- !INPUT:
-- !	  #. COMCOT.CTL AND LANDSLIDE.CTL IF REQUIRED;
-- !OUTPUT:
-- !	  #. ADDITIONAL LANDSLIDE PARAMETERS FOR LANDSLIDE CONFIGURATION;

landslideRange:: [Int]
landslideRange = [13..22]

getLandSlideParameter:: String -> [String]
getLandSlideParameter content = extractSection content landslideRange

readLandslideCtl :: IO ()
readLandslideCtl = do
    content <- readFile "landslide.ctl"
    print $ getLandSlideParameter content