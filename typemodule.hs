module TypeModule where

data MiniLayer = MiniLayer {
    hnx     :: Int,
    hny     :: Int,
    hx      :: [Double],
    hy      :: [Double],
    h       :: [[Double]],
    hp      :: [[Double]],
    hq      :: [[Double]],
    hzCurr  :: [[Double]],
    hzNext  :: [[Double]],
    hmCurr  :: [[Double]],
    hmNext  :: [[Double]],
    hnCurr  :: [[Double]],
    hnNext  :: [[Double]]
} deriving Show

updateMiniLayerH:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerH newH layer = layer {h = newH}

updateMiniLayerHP:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHP newHP layer = layer {hp = newHP}

updateMiniLayerHQ:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHQ newHQ layer = layer {hq = newHQ}

updateMiniLayerHZCurr:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHZCurr newHZ layer = layer {hzCurr = newHZ}

updateMiniLayerHZNext:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHZNext newHZ layer = layer {hzNext = newHZ}

updateMiniLayerHMCurr:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHMCurr newHM layer = layer {hmCurr = newHM}

updateMiniLayerHMNext:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHMNext newHM layer = layer {hmNext = newHM}

updateMiniLayerHNCurr:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHNCurr newHN layer = layer {hnCurr = newHN}

updateMiniLayerHNNext:: [[Double]] -> MiniLayer -> MiniLayer
updateMiniLayerHNNext newHN layer = layer {hnNext = newHN}

data Bathymetry = Bathymetry {
    bx  :: [Double],
    by  :: [Double],
    bz  :: [[Double]],
    bnx :: Int,
    bny :: Int
}

data GeneralConfig = GeneralConfig {
    total_time          :: Double,
    time_interval       :: Double,
    output_option       :: Int,
    start_type          :: Int, 
    start_time          :: Double,
    height_limit        :: Double,
    initial_surface     :: Int,
    boundary_condition  :: Int
} deriving Show


data LayerConfig = LayerConfig {
    layswitch   :: Int,
    laycord     :: Int,
    laygov      :: Int,
    dx          :: Double,
    dt          :: Double,
    fric_switch :: Int,
    fric_coef   :: Double,
    fluxswitch  :: Int,
    x_start     :: Double,
    x_end       :: Double,
    y_start     :: Double,
    y_end       :: Double,
    depth_name  :: String,
    fs          :: Int,
    layer_id    :: Int,
    level       :: Int,
    parent      :: Int
} deriving Show

updateLayerConfigDt:: Double -> LayerConfig -> LayerConfig
updateLayerConfigDt newDt config = config {dt = newDt}

data FaultConfig = FaultConfig {
    num_fault       :: Int,
    rupture_time    :: Double,
    switch          :: Int,
    focal_depth     :: Double,    
    fault_length    :: Double,
    fault_width     :: Double,
    dislocation     :: Double,
    strike_angle    :: Double,
    dip_angle       :: Double,
    slip_angle      :: Double,
    comp_origin_x   :: Double,
    comp_origin_y   :: Double,
    epicenter_x     :: Double,
    epicenter_y     :: Double,
    fault_src_data  :: String
} deriving Show

