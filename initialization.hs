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
    where layersRange = [[n+2..n+16] | n <- [84,103..331]]

getRootLayer:: String -> LayerConfig
getRootLayer content = LayerConfig {
        layswitch   = read $ head params,
        laycord     = read $ params !! 1,
        dx          = read $ params !! 2,
        dt          = read $ params !! 3,
        fric_switch = read $ params !! 4,
        fric_coef   = read $ params !! 5,
        fluxswitch  = read $ params !! 6,
        x_start     = read $ params !! 7,
        x_end       = read $ params !! 8,
        y_start     = read $ params !! 9,
        y_end       = read $ params !! 10,
        depth_name  = params !! 11,
        fs          = read $ params !! 12,
        layer_id    = read $ params !! 13,
        level       = read $ params !! 14,
        parent      = read $ params !! 15
    }
    where params = getRootLayerParameters content

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
readComcotCtl:: IO()
readComcotCtl = do
    content <- readFile "comcot.ctl"
    print $ getAllParameters content
    print $ getRootLayer content
    print $ getFaultConfig content
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