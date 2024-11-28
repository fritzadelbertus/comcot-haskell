module Output where

import TypeModule
    ( FaultConfig(fault_src_data),
      GeneralConfig(initial_surface),
      LayerConfig,
      MiniLayer(hzCurr) )
import Helper
    ( show4,
      showBoundaryCodition,
      showDepthContour,
      showStepInterval,
      showTimeHistoryOutput,
      showTimeInterval,
      showTotalRunTime,
      showZMaxOutput,
      showLayerId,
      showCoorSystem,
      showGoverningEq )
import Numeric (showFFloat)

header :: String
header = "************** COMCOT (haskell)******************\n\
\ *                                          *\n\
\ *            VERSION= 1.7                 *\n\
\ *                                          *\n\
\ ****************************************************"

--------------------------------------------------
-- Simulation Informations

generalInfo:: GeneralConfig -> Double -> Int -> String
generalInfo config step end = "------------------- GENERAL INFORMATION -------------------\n\
\     TOTAL RUN TIME            (SECOND) :  "++ showTotalRunTime config ++"\n\
\     TIME INTERVAL FOR OUTPUT  (SECOND) :  "++ showTimeInterval config ++"\n\
\     TIME STEP SIZE            (SECOND) :  "++ show4 step ++"\n\
\     TOTAL STEPS TO RUN         (STEPS) :  "++ show end ++"\n\
\     STEP INTERVAL FOR OUTPUT   (STEPS) :  "++ showStepInterval config step ++"\n\
\     TIME HISTORY RECORD OUTPUT         :  "++ showTimeHistoryOutput config ++"\n\
\     MAX. SURFACE DISPLACEMENT OUTPUT   :  "++ showZMaxOutput config ++"\n\
\     SHORELINE LOCATED AT DEPTH CONTOUR :  "++ showDepthContour config ++"\n\
\     BOUNDARY CONDITION                 :  "++ showBoundaryCodition config ++"\n\n"

faultInfoFromData:: FaultConfig -> String
faultInfoFromData config = "---------------- INITIAL CONDITION INFORMATION ---------------\n\
\     #USE INITIAL SURFACE DEFORMATION FILE:"++ fault_src_data config  ++"\n\n\n"

faultInfoFromModel:: FaultConfig -> String
faultInfoFromModel config = "NOT YET WRITTEN"

faultInfo:: GeneralConfig -> FaultConfig -> String
faultInfo gConfig fConfig 
    | initial_surface gConfig == 0 = faultInfoFromModel fConfig
    | initial_surface gConfig == 1 = faultInfoFromData fConfig
    | otherwise = "NOT AVAILABLE YET"

layerInfo:: LayerConfig -> String
layerInfo config = "--------------- GRID INFORMATION -----------------\n\
\     #  GRID IDENTIFICATION NUMBER      :  "++ showLayerId config ++"\n\
\        "++ showCoorSystem config ++"\n\
\        "++ showGoverningEq config ++"\n\
\        BOTTOM FRICTION                 : NOT AVAILABLE\n\
\        VOLUME FLUX OUTPUT              : NOT AVAILABLE\n\
\        SURFACE DISPLACEMENT OUTPUT     : ENABLED\n\
\        GRID LAYER POSITIONS            :\n\
\          X_START               (METER) :   0.00000000\n\
\          X_END                 (METER) :   99.0000000\n\
\          Y_START               (METER) :   0.00000000\n\
\          Y_END                 (METER) :   99.0000000\n\
\        GRID DIMENSION          (NX*NY) :         100 *         100\n\
\        GRID DIMENSION          (NX*NY) :         100 *         100\n\  
\        GRID SIZE           (DX, METER) :   1.00000000\n\n\
\        BATHYMETRY DATA FILE NAME       :bathy2.xyz"


simulationInfo :: LayerConfig -> FaultConfig -> GeneralConfig -> Double -> Int -> String
simulationInfo lConfig fConfig gConfig step end = "***********************************************************\n\
\ *           INPUT INFORMATION - COMCOT V1.7               *\n\
\ ***********************************************************\n\
\ "++ generalInfo gConfig step end ++"\
\ "++ faultInfo gConfig fConfig ++"\
\ "++ layerInfo lConfig


--------------------------------------------------
-- Output Free Surface Displacement

printOutput :: MiniLayer -> Int -> IO MiniLayer
printOutput layer int = do 
    let filePath = "output/output" ++ show int ++ ".dat"                   
    writeArrayToFile filePath ((showIn . hzCurr) layer)                      
    return layer

--------------------------------------------------
-- Output Helpers

writeArrayToFile :: FilePath -> [[Double]] -> IO ()
writeArrayToFile filePath array = do
    let formattedRows = unlines (map formatRow array)
    writeFile filePath formattedRows

reshape :: Int -> [a] -> [[a]]
reshape _ [] = []
reshape n xs = take n xs : reshape n (drop n xs)

showIn:: [[Double]] -> [[Double]]
showIn = foldr reduce []
    where 
        reduce item acc = acc ++ reshape 15 item

formatRow :: [Double] -> String
formatRow row = unwords (map (flip (showFFloat (Just 4)) "") row)