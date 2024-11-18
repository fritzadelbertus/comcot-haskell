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



updateMiniLayer:: (a -> MiniLayer -> MiniLayer) -> a -> MiniLayer -> MiniLayer
updateMiniLayer fieldSetter = fieldSetter 

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
    epicenter_y     :: Double
} deriving Show


data Layer = Layer {
    layerConfig     :: LayerConfig,
    name            :: String,
    nx,ny           :: Int,
    deform          :: [[Double]],      -- Seaflor deformation when fault model is implemented
    z               :: [[[Double]]],    -- Free surface elevation
    m               :: [[[Double]]],    -- Volume flux component x 
    n               :: [[[Double]]],    -- Volume flux component y
    ht              :: [[[Double]]],    -- Transient water depth at T = N*DT and T = (N+1)*DT
    dz              :: [[[Double]]],    -- Total water depth, for moving boundary
    r1              :: [[Double]],      -- Coefficients for spherical coor
    r2              :: [[Double]],
    r3              :: [[Double]],
    r4              :: [[Double]],
    r5              :: [[Double]],
    r6              :: [Double],
    r0              :: [[Double]],
    r11             :: [[Double]],
    r21             :: [[Double]],
    r22             :: [[Double]],
    xflux           :: [[Double]],         -- For flux interpolation
    yflux           :: [[Double]],
    x               :: [Double],           -- X coordinate of grid
    y               :: [Double],           -- Y coordinate of grid
    xt              :: [Double],           -- X coordinate of grid (temp)
    yt              :: [Double],           -- Y coordinate of grid (temp)
    dy              :: Double,
    del_x           :: [Double],
    del_y           :: [Double],
    rx              :: Double,
    ry              :: Double,
    grx             :: Double,
    gry             :: Double,
    ini_switch      :: Int,
    bc_type         :: Int,
    dim             :: Int,
    linchk          :: Int,
    modscm          :: Int,
    fric_vcoef      :: [[Double]],
    mask            :: [[Int]],
    alpha           :: [[Double]],
    a1x             :: [[Double]],
    a2x             :: [[Double]],
    a1y             :: [[Double]],
    a2y             :: [[Double]],
    cf              :: [[[Double]]],
    cb              :: [[[Double]]],
    m0              :: [[Double]],
    n0              :: [[Double]],
    sponge_coefx    :: [[Double]],
    sponge_coefy    :: [[Double]],
    z_max           :: [[Double]],
    z_min       :: [[Double]],
    south_lat   :: Double,
    rel_size    :: Int,
    rel_time    :: Int,
    num_child   :: Int,
    corners     :: [Int],
    xo          :: Double,
    yo          :: Double,
    h_limit     :: Double,
    tide_level  :: Double,
    upz         :: Bool,
    sc_option   :: Int,
    pos         :: [[[Int]]],
    cxy         :: [[[Double]]],
    children    :: [Layer]
} deriving Show


updateLayerZ:: [[[Double]]] -> Layer -> Layer
updateLayerZ newZ layer = layer {z = newZ}

updateLayerDz:: [[[Double]]] -> Layer -> Layer
updateLayerDz newDz layer = layer {dz = newDz}

updateLayerM:: [[[Double]]] -> Layer -> Layer
updateLayerM newM layer = layer {m = newM}

updateLayerN:: [[[Double]]] -> Layer -> Layer
updateLayerN newN layer = layer {n = newN}

updateLayerDeform:: [[Double]] -> Layer -> Layer
updateLayerDeform newDeform layer = layer {deform = newDeform}

updateLayerNx:: Int -> Layer -> Layer
updateLayerNx newNx layer = layer {nx = newNx}

updateLayerNy:: Int -> Layer -> Layer
updateLayerNy newNy layer = layer {ny = newNy}

updateLayerX:: [Double] -> Layer -> Layer
updateLayerX newX layer = layer {x = newX}

updateLayerY:: [Double] -> Layer -> Layer
updateLayerY newY layer = layer {y = newY}

updateLayerXO:: Double -> Layer -> Layer
updateLayerXO newXO layer = layer {xo = newXO}

updateLayerYO:: Double -> Layer -> Layer
updateLayerYO newYO layer = layer {yo = newYO}

updateLayerChildren:: [Layer] -> Layer -> Layer
updateLayerChildren newChildren layer = layer {children = newChildren}

