module Deform where
import TypeModule
import Constants
import Helper

angleL :: FaultConfig -> Double
angleL = (radDeg *) . dip_angle

angleR :: FaultConfig -> Double
angleR = (radDeg *) . slip_angle

angleT :: FaultConfig -> Double
angleT = (radDeg *) . strike_angle

halfFault :: FaultConfig -> Double
halfFault = (0.5 *) . fault_length

hh :: FaultConfig -> Double
hh = onResults1 (+) focal_depth ((0.5 *) . onResults1 (*) fault_width (sin . angleL))

delx :: FaultConfig -> Double
delx = (0.5 *) . onResults1 (*) fault_width (onResults1 (*) (cos . angleL) (cos .angleT))

dely :: FaultConfig -> Double
dely = (0.5 *) . onResults1 (*) fault_width (onResults1 (*) (cos . angleL) (sin . angleT))

projectionX :: Double -> Double -> Double -> Double -> Double
projectionX = okadaProjection fst

projectionY :: Double -> Double -> Double -> Double -> Double
projectionY = okadaProjection snd

okadaProjection :: ((Double,Double) -> Double) -> Double -> Double -> Double -> Double -> Double
okadaProjection f = (((f .) .) .) . stereoProjection

slipDislocation1 :: FaultConfig -> Double
slipDislocation1 = onResults1 (*) dislocation (cos . angleR)

slipDislocation2 :: FaultConfig -> Double
slipDislocation2 = onResults1 (*) dislocation (sin . angleR)


xShift :: FaultConfig -> MiniLayer -> LayerConfig -> Int -> Int -> Double
xShift fault l c i j 
    | isSpherical c = projXS - delx fault
    | otherwise = xx - projXC - delx fault
    where 
        projXS = projectionX (hx l !! i) (hy l !! j) (epicenter_x fault) (epicenter_y fault)
        projXC = projectionX (comp_origin_x fault) (comp_origin_y fault) (epicenter_x fault) (epicenter_y fault)
        xx = dx c * (fromIntegral i - 1.0)

yShift :: FaultConfig -> MiniLayer -> LayerConfig -> Int -> Int -> Double
yShift fault l c i j
    | isSpherical c = projYS - dely fault
    | otherwise = yy - projYC - dely fault
    where 
        projYS = projectionY (hx l !! i) (hy l !! j) (epicenter_x fault) (epicenter_y fault)
        projYC = projectionY (comp_origin_x fault) (comp_origin_y fault) (epicenter_x fault) (epicenter_y fault)
        yy = dx c * (fromIntegral j - 1.0)

okadaDisplacement:: FaultConfig -> MiniLayer -> LayerConfig -> Int -> Int -> Double
okadaDisplacement fault l c i j = us + ud
    where 
        us
            | abs usT2 <= gx = zero
            | otherwise = usT2
        ud
            | abs udT2 <= gx = zero
            | otherwise = udT2
        usT1 = (f1-f2-f3+f4)*slipDislocation1 fault
        udT1 = (g1-g2-g3+g4)*slipDislocation2 fault
        usT2 = setToZeroIfNanOrInfinite usT1
        udT2 = setToZeroIfNanOrInfinite udT1
        x1 = (xShift fault l c i j)*sin (angleT fault) + (yShift fault l c i j)*cos (angleT fault)
        x2 = (yShift fault l c i j)*sin (angleT fault) - (xShift fault l c i j)*cos (angleT fault)
        p = x2*cs+hh fault*sn 
        f1 = strikeSlip x2 (x1+halfFault fault) p (angleL fault) (hh fault)
        f2 = strikeSlip x2 (x1+halfFault fault) (p-fault_width fault) (angleL fault) (hh fault)
        f3 = strikeSlip x2 (x1-halfFault fault) p (angleL fault) (hh fault)
        f4 = strikeSlip x2 (x1-halfFault fault) (p-fault_width fault) (angleL fault) (hh fault)
        g1 = strikeSlip x2 (x1+halfFault fault) p (angleL fault) (hh fault)
        g2 = strikeSlip x2 (x1+halfFault fault) (p-fault_width fault) (angleL fault) (hh fault)
        g3 = strikeSlip x2 (x1-halfFault fault) p (angleL fault) (hh fault)
        g4 = strikeSlip x2 (x1-halfFault fault) (p-fault_width fault) (angleL fault) (hh fault)
        sn = (sin . angleL) fault
        cs = (cos .angleL) fault

deformOkada:: MiniLayer -> LayerConfig -> FaultConfig -> [[Double]]
deformOkada layer layConfig fault = generateLayerWithConfig layer layConfig (okadaDisplacement fault)

faultDBar:: Double -> Double -> Double -> Double -> Double -> Double
faultDBar x2 y1 y2 dd = onResults1 (-) ((*) y2 . sin) (onResults1 (*) (faultQ x2 y1 y2 dd) cos)
faultQ:: Double -> Double -> Double -> Double -> Double -> Double
faultQ x2 y1 y2 dd = onResults1 (-) ((*) x2 . sin) ((*) dd . cos)
faultR:: Double -> Double -> Double -> Double -> Double -> Double
faultR x2 y1 y2 dd = ((sqrt . (sum . map (**2))) . flip (:) [y1,y2]) . faultQ x2 y1 y2 dd

strikeSlip:: Double -> Double -> Double -> Double -> Double -> Double
strikeSlip x2 y1 y2 dp dd = - ((dBar*q/r/(r+y2) + q*sin dp/(r+y2) + a4*sin dp)/2/pi)
    where 
        dBar = faultDBar x2 y1 y2 dd dp
        r = faultR x2 y1 y2 dd dp
        q = faultQ x2 y1 y2 dd dp
        a4 = 0.5/cos dp*(log tmp1 - sin dp*log tmp2)
        tmp1
            | r+dBar < eps = eps
            | otherwise = r+dBar
        tmp2
            | r+y2 < eps = eps
            | otherwise = r+y2

dipSlip:: Double -> Double -> Double -> Double -> Double -> Double
dipSlip x2 y1 y2 dp dd = - ((dBar*q/r/(r+y1) + sin dp*atan (y1*y2/q/r) - a5*sin dp*cos dp)/2/pi)
    where 
        dBar = faultDBar x2 y1 y2 dd dp
        r = faultR x2 y1 y2 dd dp
        q = faultQ x2 y1 y2 dd dp
        a5 = 0.5*2/cos dp*atan ((y2*(xx+q*cos dp)+xx*(r+xx)*sin dp)/y1/(r+xx)/cos dp)
        xx = sqrt (y1**2 + q**2)
        

stereoProjection:: Double -> Double -> Double -> Double -> (Double,Double)
stereoProjection lonin latin lon0 lat0 = (calcX, calcY)
    where
        calcX = xf + 2.0*r*k0*cos xi*sin (lm-lm0)/beta
        calcY = yf + 2.0*r*k0*((sin xi*cos xi0) - (cos xi*sin xi0*cos (lm-lm0)))/beta
        
        beta = 1.0 + (sin xi*sin xi0) + (cos xi*cos xi0*cos (lm-lm0))
        xi = asin ((w-1.0)/(w+1.0))
        lm = n*(lon-lm0)+lm0
        lm0 = ln0
        xi0 = asin ((w2-1.0)/(w2+1.0))
        
        w = c*((sa*(sb**e))**n)
        sa = (1.0+sn)/(1.0-sn)
        sb = (1.0 - e*sn)/(1.0 + e*sn)
        w2 = c*w1

        c = (n+sn0)*(1.0-snxi0)/((n-sn0)*(1.0+snxi0))
        snxi0 = (w1-1.0)/(w1+1.0)
        w1 = (s1*(s2**e))**n
        s2 = (1.0-e*sn0)/(1.0+e*sn0)
        s1 = (1.0+sn0)/(1.0-sn0)

        n = sqrt (1.0+f2*cs0**4/(1.0-f2))
        r = sqrt (rho0*nu0)
        nu0 = a/tmp0
        rho0 = a*(1.0-f2)/tmp0**3
        tmp0 = sqrt (1.0-f2*sn0**2)
        tmp = sqrt (1.0-f2*sn**2)

        xf = 0.0
        yf = 0.0

        a = 6378137.0000          
        b = 6356752.3142
        f = 0.00335281067183    
        e = 0.08181919092891
        f2 = 0.00669438000426
        es = 0.00673949675659

        k0 = 0.9996

        cs = cos lat
        sn = sin lat
        cs0 = cos lt0
        sn0 = sin lt0

        pole = pi/2.0 - eps

        lat
            | latin*radDeg > pole = pole
            | latin*radDeg < (-pole) = -pole
            | otherwise = latin*radDeg
        lt0
            | lat0*radDeg > pole = pole
            | lat0*radDeg < (-pole) = -pole
            | otherwise = lat0*radDeg
        lon = lonin*radDeg
        ln0 = lon0*radDeg
