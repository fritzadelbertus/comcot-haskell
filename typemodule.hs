module TypeModule where

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

data FaultConfig = FaultConfig {
    hh          :: Double,    
    l           :: Double,
    w           :: Double,
    d           :: Double,
    th          :: Double,
    dl          :: Double,
    rd          :: Double,
    fault_xo    :: Double,
    fault_yo    :: Double,
    x0          :: Double,
    y0          :: Double,
    t0          :: Double,
    switch      :: Int,
    num_flt     :: Int,
    fault_fs    :: Int,
    deform_name :: String
}


data Layer = Layer {
    layerConfig     :: LayerConfig,
    name            :: String,
    nx,ny           :: Int,
    deform          :: [[Double]],      -- Seaflor deformation when fault model is implemented
    z               :: [[[Double]]],    -- Free surface elevation
    m               :: [[[Double]]],    -- Volume flux component x 
    n               :: [[[Double]]],    -- Volume flux component y
    h               :: [[Double]],      -- Still water depth (no change)
    ht              :: [[[Double]]],    -- Transient water depth at T = N*DT and T = (N+1)*DT
    hp              :: [[Double]],      -- Still water depth
    hq              :: [[Double]],      -- Still water depth
    dh              :: [[[Double]]],    -- Sediment transport
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
    cxy         :: [[[Double]]]
} deriving Show