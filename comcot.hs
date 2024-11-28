import Initialization
import TypeModule
    ( updateMiniLayerHMCurr,
      updateMiniLayerHNCurr,
      updateMiniLayerHZCurr,
      FaultConfig,
      GeneralConfig(total_time),
      LayerConfig(dt, depth_name),
      MiniLayer(hzNext, hnNext, hmNext) )
import Mass ( mass )
import Moment ( moment )
import Boundary ( open )
import Output ( printOutput, header, simulationInfo )


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
        nextLayer = change 
            $ moment layConfig                                          -- Moment.hs 
            $ open                                                      -- Boundary.hs
            $ mass layConfig layer                                      -- Mass.hs


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
    initdeformLayer <- getInitialSurface initialLayer initlayConfig faultConfig     -- Initialization.hs
    putStrLn "ADJUSTING BATHYMETRY DATA..."
    let adjustedBLayer = adjustBathymetry initdeformLayer                           -- Initialization.hs
    

    let layConfig = checkCourantCondition adjustedBLayer initlayConfig              -- Initialization.hs
    

    let istart = 1 
    let time = 0.0
    let kEnd = round (total_time generalConfig/dt layConfig)

    putStrLn $ simulationInfo layConfig faultConfig generalConfig (dt layConfig) kEnd -- Output.hs

    lastOutput <- iterateComcot istart kEnd time adjustedBLayer layConfig faultConfig -- Main Iteration
    putStrLn "Simulation Finished."
    

    

