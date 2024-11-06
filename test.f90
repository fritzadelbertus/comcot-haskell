program main
    implicit none
    real :: x,y

    ! Call the subroutine
    call STEREO_PROJECTION(x,y,10.5, 75.2, 9.0, 76.5)

    ! Print the result
    print *, "The result is:", x,y
end program main

SUBROUTINE STEREO_PROJECTION (X,Y,LONIN,LATIN,LON0,LAT0)
    !......................................................................
    !DESCRIPTION:
    !     # MAPPING A POINT ON THE ELLIPSOID SURFACE ONTO A PLANE;
    !     # OBLIQUE STEREOGRAPHIC PROJECTION IS ADOPTED
    !INPUT:
    !     LATIN: LATITUDE IN DEGREES
    !     LONIN: LONGITUDE IN DEGREES
    !     LAT0: LATITUDE OF TANGENTIAL POINT IN DEGREES (E.G., EPICENTER)
    !     LON0: LONGITUDE OF TANGENTIAL POINT IN DEGREES (E.G., EPICENTER)
    !OUTPUT:
    !     X: X COORDINATE/EASTING IN METERS RELATIVE TO ORIGIN (I.E., LON0)
    !     Y: Y COORDINATE/NORTHING IN METERS RELATIVE TO ORIGIN (I.E., LAT0)
    !REFERENCES:
    !	  #. SNYDER, J.P. (1987). MAP PROJECTIONS - A WORKING MANUAL.
    !                          USGS PROFESSIONAL PAPER 1395
    !     #. ELLIPSOIDAL DATUM: WGS84
    !WORKING NOTES:
    !     CREATED ON DEC18 2008 (XIAOMING WANG, GNS)
    !     UPDATED ON JAN02 2009 (XIAOMING WANG, GNS)
    !----------------------------------------------------------------------
          REAL XOUT,YOUT,XS,YS,LAT0,LON0,LATIN,LONIN
          REAL COS_X,COS_Y,SIN_X,SIN_Y
          REAL LAT,LON,LT0,LN0,CS,SN,CS0,SN0,TMP,TMP0
          REAL A,B,K0,E,ES,N,C,R,S1,S2,W1,W2,SA,SB,BETA,XI,LM,XI0,LM0
    
          POLE = 3.141592653589793/2.0 - 10**(-10)	  !AVOID SINGULARITY AT POLES
    
          ! CONVERT DEGREE TO RADIAN
          LAT = LATIN*0.01745329252
          LON = LONIN*0.01745329252
          LT0 = LAT0*0.01745329252
          LN0 = LON0*0.01745329252
          IF (LAT .GT. POLE) LAT = POLE
          IF (LAT .LT. -POLE) LAT = -POLE
          IF (LT0 .GT. POLE) LT0 = POLE
          IF (LT0 .LT. -POLE) LT0 = -POLE
    
          CS = COS(LAT)
          SN = SIN(LAT)
          CS0 = COS(LT0)
          SN0 = SIN(LT0)
    
          ! PARAMETERS
          XF = 0.0					! FALSE EASTING
          YF = 0.0					! FALSE NORTHING
    
          A = 6378137.0000          ! ELLIPSOIDAL SEMI-MAJOR AXIS
          B = 6356752.3142			! ELLIPSOIDAL SEMI-MINOR AXIS
          F = 0.00335281067183      ! FLATTENING, F = (A-B)/A
          E = 0.08181919092891		! ECCENTRICITY, E = SQRT(2.0*F-F**2)
          F2 = 0.00669438000426		! F2 = E**2
          ES = 0.00673949675659		! 2ND ECCENTRICITY, ES = E**2/(1-E**2)
    
          K0 = 0.9996				! SCALE FACTOR
    
          TMP = SQRT(1.0-F2*SN**2)
          TMP0 = SQRT(1.0-F2*SN0**2)
          RHO0 = A*(1.0-F2)/TMP0**3
          NU0 = A/TMP0
          R = SQRT(RHO0*NU0)
          N = SQRT(1.0+F2*CS0**4/(1.0-F2))
    
          S1 = (1.0+SN0)/(1.0-SN0)
          S2 = (1.0-E*SN0)/(1.0+E*SN0)
          W1 = (S1*S2**E)**N
          SN_XI0 = (W1-1.0)/(W1+1.0)
          C = (N+SN0)*(1.0-SN_XI0)/(N-SN0)/(1.0+SN_XI0)
    
          W2 = C*W1
          SA = (1.0+SN)/(1.0-SN)
          SB = (1.0-E*SN)/(1.0+E*SN)
          W = C*(SA*SB**E)**N
    
          XI0 = ASIN((W2-1.0)/(W2+1.0))
          LM0 = LN0
    
          LM = N*(LON-LM0)+LM0
          XI = ASIN((W-1.0)/(W+1.0))
    
          BETA = 1.0 + SIN(XI)*SIN(XI0) + COS(XI)*COS(XI0)*COS(LM-LM0)
    
          Y = YF + 2.0*R*K0*(SIN(XI)*COS(XI0) &
                             - COS(XI)*SIN(XI0)*COS(LM-LM0))/BETA
          X = XF + 2.0*R*K0*COS(XI)*SIN(LM-LM0)/BETA
    
    !	  WRITE(*,*) LM,XI,X,Y
    
          RETURN
          END