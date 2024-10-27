-- ADA 25 "BIG" Functions pada initialization.f90

import Helper (removeLeadingSpaces)

-- Constants
parameterIndent:: Int
parameterIndent = 49

getParameterValue:: String -> Int -> String
getParameterValue content line = drop parameterIndent $ lines content !! line

extractSection:: String -> [Int] -> [String]
extractSection content = map (removeLeadingSpaces . getParameterValue content)

-- 1. READ CONFIG (FILE READING)
-- !DESCRIPTION:
-- !	  #. OBTAIN ALL THE PARAMETERS FROM COMCOT.CTL;
-- !	  #. START_TYPE =
-- !				0: COLD START (SIMULATION STARTS FROM T = 0)
-- !				1: HOT START (SIMULATION STARTS FROM RESUMING TIME)
-- !				20: COLD START WITH TIDE LEVEL ADJUSTMENT
-- !				21: HOT START WITH TIDE LEVEL ADJUSTMENT
-- !	  #. INI_SURF =
-- !				0: USE FAULT MODEL TO DETERMINE SEAFLOOR DEFORMATION
-- !				1: USE DATA FILE TO DETERMINE INITIAL WATER SURFACE
-- !				2: USE INCIDENT WAVE MODEL TO GENERATE WAVES
-- !				3: USE TRANSIENT FLOOR MOTION MODEL (LANDSLIDE);
-- !				4: USE FAULT MODEL + LANDSLIDE;
-- !				   FAULT_MULTI.CTL IS REQUIRED FOR MULTI-FAULT SETUP;
-- !				9: USE MANSINHA AND SMYLIES' MODEL TO CALC DEFORMATION
-- !INPUT:
-- !	  #. COMCOT.CTL (AND FAULT_MULTI.CTL FOR MORE THAN ONE FAULT PLANE)
-- !OUTPUT:
-- !	  #. GENERAL INFORMAITON FOR A SIMULATION;
-- !      #. GRID SETUP;

sectionRange:: [[Int]] -- General, Fault Model, Wave Maker, Landslide
sectionRange = [[10..20],[25..40],[45..49],[54..59]]

rootLayerRange :: [Int]
rootLayerRange = [66..82]

layersRange:: [[Int]]
layersRange = [[n+2..n+17] | n <- [85,105..345]]

getAllParameters:: String -> [[String]]
getAllParameters content = map (extractSection content) sectionRange

getRootLayerParameters:: String -> [String]
getRootLayerParameters content = extractSection content rootLayerRange

getChildLayersParameters:: String -> [[String]]
getChildLayersParameters content = map (extractSection content) layersRange

-- main :: IO ()
-- main = do
--     content <- readFile "comcot.ctl"
--     print $ getAllParameters content
--     print $ getRootLayerParameters content
--     print $ layersRange
--     print $ getChildLayersParameters content

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

-- main :: IO ()
-- main = do
--     content <- readFile "fault_multi.ctl"
--     print $ getAllFaultModelParameter content

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

getLandSlideParameter:: String -> [String]
getLandSlideParameter content = extractSection content [13..22]

main :: IO ()
main = do
    content <- readFile "landslide.ctl"
    print $ getLandSlideParameter content