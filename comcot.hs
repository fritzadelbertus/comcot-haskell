import Initialization
import TypeModule
import Mass
import Moment

change:: MiniLayer -> MiniLayer
change l = (updateMiniLayerHNCurr (hnNext l) . updateMiniLayerHMCurr (hmNext l). updateMiniLayerHZCurr (hzNext l)) l

iterateComcot :: Int -> Int -> Double -> MiniLayer ->  LayerConfig -> FaultConfig -> MiniLayer
iterateComcot start end time layer layConfig fault
    | start > end = layer
    | otherwise = iterateComcot (start+1) end nextTime nextLayer layConfig fault
    where 
        nextTime = time + dt layConfig
        nextLayer = change $ moment layConfig $ open $ mass layConfig $ getFloorDeform nextTime layer layConfig fault

main:: IO()
main = do
    generalConfig <- readGeneralConfig "comcot.ctl"     -- Initialization.hs
    initlayConfig <- readLayerConfig "comcot.ctl"           -- Initialization.hs
    faultConfig <- readFaultConfig "comcot.ctl"         -- Initialization.hs
    initialLayer <- getInitialLayer initlayConfig           -- Initialization.hs
    let initdeformLayer = getInitialSurface initialLayer initlayConfig faultConfig -- Initialization.hs
    let adjustedBLayer = adjustBathymetry initdeformLayer  -- Initialization.hs
    let layConfig = checkCourantCondition adjustedBLayer initlayConfig
    -- sphericalParameters
    let istart = 1
    let time = 0.0
    let kEnd = round (total_time generalConfig/dt layConfig)
    print adjustedBLayer
    print layConfig
    let lastOutput = iterateComcot istart kEnd time adjustedBLayer layConfig faultConfig
    print lastOutput
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

    
