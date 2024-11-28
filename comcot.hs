import Initialization
import TypeModule
import Helper
import Mass
import Moment
import Boundary
import Output


change:: MiniLayer -> MiniLayer
change l = (updateMiniLayerHNCurr (hnNext l) . updateMiniLayerHMCurr (hmNext l). updateMiniLayerHZCurr (hzNext l)) l

iterateComcot :: Int -> Int -> Double -> MiniLayer ->  LayerConfig -> FaultConfig -> IO MiniLayer
iterateComcot start end time layer layConfig fault
    | start > end = return layer
    | start `mod` 2 == 0 = do
        layer <- printOutput layer start
        iterateComcot (start+1) end nextTime nextLayer layConfig fault
    | otherwise = do 
        iterateComcot (start+1) end nextTime nextLayer layConfig fault
    where  
        nextTime = time + dt layConfig
        nextLayer = change $ moment layConfig $ open $ mass layConfig layer 


main:: IO()
main = do
    putStrLn header
    putStrLn "READING PARAMETERS FOR SIMULATION..."
    putStrLn "    READING GENERAL INFORMATION......"
    generalConfig <- readGeneralConfig "comcot.ctl"                                 -- Initialization.hs
    putStrLn "    READING PARAMETERS FOR FAULT MODEL......"
    faultConfig <- readFaultConfig "comcot.ctl"                                     -- Initialization.hs
    putStrLn "    READING PARAMETERS FOR GRID LAYER......"
    initlayConfig <- readLayerConfig "comcot.ctl"                                   -- Initialization.hs
    putStrLn $ depth_name initlayConfig
    
    putStrLn "READING BATHYMETRY DATA..."
    initialLayer <- getInitialLayer initlayConfig                                   -- Initialization.hs
    initdeformLayer <- getInitialSurface initialLayer initlayConfig faultConfig  -- Initialization.hs
    putStrLn "ADJUSTING BATHYMETRY DATA..."
    let adjustedBLayer = adjustBathymetry initdeformLayer                           -- Initialization.hs
    

    let layConfig = checkCourantCondition adjustedBLayer initlayConfig
    

    let istart = 1 
    let time = 0.0
    let kEnd = round (total_time generalConfig/dt layConfig)

    putStrLn $ simulationInfo layConfig faultConfig generalConfig (dt layConfig) kEnd

    lastOutput <- iterateComcot istart kEnd time adjustedBLayer layConfig faultConfig
    putStrLn "Simulation Finished."
    
    -- screenDisplay
    -- for k = 1 to end
    -- if k mod 10 == 0 print k,minutes
    -- getFloorDeform
    -- mass
    -- open
    -- if length children layer >= 1 allGrid
    -- moment
    -- if bcType 1 spongeLayer
    -- if bcType 2 bcWall
    -- if k mod iprt == 0 generateFile
    -- maxAmp
    -- change
    --getTs

    

