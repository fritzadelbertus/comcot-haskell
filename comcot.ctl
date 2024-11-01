#################################################
#                                               #
# Control file for COMCOT program (v1.7)        #
# Sample                                        #
#################################################
#--+-----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
#===============================================:===============================
# General Parameters for Simulation             : Value Field                  |
#===============================================:===============================
#Job Description: NZ30sec bathymetry, Spherical Coordinates for code testing
 Total run time (Wall clock, seconds)           :  7200.000
 Time interval to Save Data    ( unit: sec )    :    60.0
 Output Zmax & TS (0-Max Z;1-Timeseries;2-Both) :     2
 Start Type (0-Cold start; 1-Hot start)         :     0
 Resuming Time If hot start (Seconds)           :  8000.00
 Specify Min WaterDepth offshore  (meter)       :     0.00
 Initial Cond. (0:FLT,1:File,2:WM,3:LS,4:FLT+LS):     0
 Specify BC  (0-Open;1-Sponge;2-Wall;3-FACTS)   :     0
 Specify Input Z filename (for BC=3, FACTS)     : mw94_n22_nz_ha.xyt
 Specify Input U filename (for BC=3, FACTS)     : mw94_n22_nz_ua.xyt
 Specify Input V filename (for BC=3, FACTS)     : mw94_n22_nz_va.xyt

#===============================================:===============================
# Parameters for Fault Model (Segment 01)       :Values                        |
#===============================================:===============================
 No. of FLT Planes (With fault_multi.ctl if >1) :   1
 Fault Rupture Time (seconds)                   :   60.0
 Faulting Option (0: Model; 1- Data;)           :   0
 Focal Depth                             (meter):   30000.000
 Length of source area                   (meter):  100000.000
 Width of source area                    (meter):   50000.000
 Dislocation of fault plate              (meter):       5.000
 Strike direction (theta)               (degree):     270.000
 Dip  angle       (delta)               (degree):      30.000
 Slip angle       (lamda)               (degree):      10.000
 Origin of Comp. Domain (Layer 01) (Lat, degree): -2.5
 Origin of Comp. Domain (Layer 01) (Lon, degree): 100.5
 Epicenter: Latitude                    (degree): -3.0
 Epicenter: Longitude                   (degree): 100.0
 File Name of Deformation Data                  :comcot_ini.dep
 Data Format Option (0-COMCOT; 1-MOST; 2-XYZ)   :     2

#===============================================:===============================
#  Parameters for Wave Maker                    :Values                        |
#===============================================:===============================
 Wave type  ( 1:Solit, 2:given, 3:focusing )    :     1
 FileName of Customized Input (for Type=2)      : fse.dat
 Incident direction( 1:top,2:bt,3:lf,4:rt,5:ob ):     2
 Characteristic Wave Amplitude        (meter)   :     0.500
 Typical Water depth                  (meter)   :  2000.000

#===============================================:===============================
#  Parameters for Submarine LS/Transient Motion :ValUes                        |
#===============================================:===============================
 X Coord. of Left/West Edge of Landlide Area    :  177.00
 X Coord. of Right/East Edge of Landlide Area   :  179.00
 Y Coord. of Bottom/South Edge of Landlide Area :  -41.00
 Y Coord. of Top/North Edge of Landlide Area    :  -39.00
 File Name of landslide Data                    : landslide_test.dat
 Data Format Option (0-Old; 1-XYT; 2-Function)  :     2

#===============================================:===============================
# Configurations for all grids                  :Values                        |
#===============================================:===============================
# Parameters for 1st-level grid -- layer 01     :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     0
 Coordinate System    (0:spherical, 1:cartesian):     0
 Governing Equations  (0:linear,    1:nonlinear):     1
 Grid Size  (dx, sph:minute, cart:meter)        :     0.5
 Time step                            ( second ):     1.0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 X_start                                        :  98.5041666666666629
 X_end                                          :  101.995833333333337
 Y_Start                                        :  -5.99583333333333357
 Y_end                                          :  -2.00416666666666687
 File Name of Bathymetry Data                   :mentawaio.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    01
 Grid Level                                     :     1
 Parent Grid's ID Number                        :    -1

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 02    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 GridSize Ratio of Parent layer to current layer:     2
 X_start                                        :   98.5
 X_end                                          :   102.0
 Y_start                                        :   -6.0
 Y_end                                          :   -2.0
 FileName of Water depth data                   :mentawaio.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    02
 Grid Level                                     :     2
 Parent Grid's ID Number                        :    01

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 03    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 GridSize Ratio of Parent layer to current layer:     4
 X_start                                        :  98.5
 X_end                                          :  102.0
 Y_start                                        :   -6.0
 Y_end                                          :   -2.0
 FileName of Water depth data                   :mentawaio.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    03
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 04    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.6354616225919
 X_end                                          :   120.9234496509648
 Y_start                                        :    21.84426599282262
 Y_end                                          :    22.21347088069079
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    04
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 05    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.3092711496687
 X_end                                          :   120.723621321133
 Y_start                                        :    22.1896412622199
 Y_end                                          :    22.53257020096585
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    05
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 06    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.0036509071348
 X_end                                          :   120.5120381818606
 Y_start                                        :    22.42373978748143
 Y_end                                          :    24.24865827343963
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    06
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 07    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.4562042092812
 X_end                                          :   120.9175729323763
 Y_start                                        :    24.15490301533408
 Y_end                                          :    24.80244805853942
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    07
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 08    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.808842317314
 X_end                                          :   121.3231069967601
 Y_start                                        :    24.69987279481362
 Y_end                                          :    25.19808824930069
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    08
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 09    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   121.2320085958588
 X_end                                          :   122.0460161131184
 Y_start                                        :    24.96080357036498
 Y_end                                          :    25.3508512831777
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    09
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 10    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   121.7227643425554
 X_end                                          :   122.0489548154786
 Y_start                                        :    24.22894457248688
 Y_end                                          :    25.0404101584892
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    10
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 11    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   121.3906957787801
 X_end                                          :   121.8079853387364
 Y_start                                        :    23.41747516260113
 Y_end                                          :    24.38422152963531
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    11
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 12    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   121.1203399645682
 X_end                                          :   121.4906102867619
 Y_start                                        :    22.74369554382476
 Y_end                                          :    23.51418785598371
 FileName of Water depth data                   :  etopo_indo.nc
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    12
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 13    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.808842317314
 X_end                                          :   121.2084996631093
 Y_start                                        :    22.14603554034114
 Y_end                                          :    22.78475484366331
 FileName of Water depth data                   :  4_10.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    13
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 14    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.808842317314
 X_end                                          :   121.2084996631093
 Y_start                                        :    22.14603554034114
 Y_end                                          :    22.78475484366331
 FileName of Water depth data                   :  4_10.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    14
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02

#===============================================:===============================
#  Parameters for Sub-level grid -- layer 15    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     0
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     16
 X_start                                        :   120.808842317314
 X_end                                          :   121.2084996631093
 Y_start                                        :    22.14603554034114
 Y_end                                          :    22.78475484366331
 FileName of Water depth data                   :  4_10.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     3
 Grid Identification Number                     :    15
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02
