module Deform where
import TypeModule
import Constants
import Helper

deformOkada:: Layer -> FaultConfig -> Layer
deformOkada layer fault = updateLayerDeform okadaDeform layer
    where 
        okadaDeform = [[eta layer i j | i <- [1..lastx]] | j <- [1..lasty]]
        lastx = nx layer
        lasty = ny layer
        eta l i j = us + ud
            where 
                angL = radDeg*dip_angle fault
                angR = radDeg*slip_angle fault
                angT = radDeg*strike_angle fault
                halfl = 0.5*fault_length fault
                hh = focal_depth fault + 0.5*fault_width fault*sin angL
                delx = 0.5*fault_width fault*cos angL * cos angT
                dely = 0.5*fault_width fault*cos angL * sin angT
                projXS = fst $ stereoProjection (x l !! i) (y l !! j) (epicenter_x fault) (epicenter_y fault)
                projYS = snd $ stereoProjection (x l !! i) (y l !! j) (epicenter_x fault) (epicenter_y fault)
                projXC = fst $ stereoProjection (comp_origin_x fault) (comp_origin_y fault) (epicenter_x fault) (epicenter_y fault)
                projYC = snd $ stereoProjection (comp_origin_x fault) (comp_origin_y fault) (epicenter_x fault) (epicenter_y fault)
                xx = dx (layerConfig layer) * (fromIntegral i - 1.0)
                yy = dy layer * (fromIntegral j - 1.0)
                xShift
                    | laycord (layerConfig l) == 0 = projXS - delx
                    | otherwise = xx - projXC - delx
                yShift
                    | laycord (layerConfig l) == 0 = projYS - dely
                    | otherwise = yy - projYC - dely
                us
                    | abs usT2 <= gx = zero
                    | otherwise = usT2
                ud
                    | abs udT2 <= gx = zero
                    | otherwise = udT2
                usT1 = (f1-f2-f3+f4)*ds
                udT1 = (g1-g2-g3+g4)*dd
                usT2 = setToZeroIfNanOrInfinite usT1
                udT2 = setToZeroIfNanOrInfinite udT1
                x1 = xShift*sin angT + yShift*cos angT
                x2 = yShift*sin angT - xShift*cos angT
                x3 = zero
                p = x2*cs+hh*sn 
                f1 = strikeSlip x1 x2 x3 (x1+halfl) p angL hh
                f2 = strikeSlip x1 x2 x3 (x1+halfl) (p-fault_width fault) angL hh
                f3 = strikeSlip x1 x2 x3 (x1-halfl) p angL hh
                f4 = strikeSlip x1 x2 x3 (x1-halfl) (p-fault_width fault) angL hh
                g1 = strikeSlip x1 x2 x3 (x1+halfl) p angL hh
                g2 = strikeSlip x1 x2 x3 (x1+halfl) (p-fault_width fault) angL hh
                g3 = strikeSlip x1 x2 x3 (x1-halfl) p angL hh
                g4 = strikeSlip x1 x2 x3 (x1-halfl) (p-fault_width fault) angL hh
                ds = dislocation fault * cos angR
                dd = dislocation fault * sin angR
                sn = sin angL
                cs = cos angL
                

strikeSlip:: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
strikeSlip x1 x2 x3 y1 y2 dp dd = - (dBar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/2/pi
    where 
        dBar = y2*sn - q*cs
        sn = sin dp
        cs = cos dp
        r = sqrt (y1**2 + y2**2 + q**2)
        q = x2*sn - dd*cs
        a4 = 0.5/cs*(log tmp1 - sn*log tmp2)
        tmp1
            | r+dBar < eps = eps
            | otherwise = r+dBar
        tmp2
            | r+y2 < eps = eps
            | otherwise = r+y2

dipSlip:: Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
dipSlip x1 x2 x3 y1 y2 dp dd = - (dBar*q/r/(r+y1) + sn*atan (y1*y2/q/r) - a5*sn*cs)/2/pi
    where 
        dBar = y2*sn - q*cs
        sn = sin dp
        cs = cos dp
        r = sqrt (y1**2 + y2**2 + q**2)
        q = x2*sn - dd*cs
        a5 = 0.5*2/cs*atan ((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
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
