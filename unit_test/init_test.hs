import Initialization (
        readComcotCtl,
        readFaultMultiCtl,
        readLandslideCtl
    )


main :: IO ()
main = do
    -- print "READING COMCOT.CTL..."
    -- readComcotCtl
    -- print "READING FAULT_MULTI.CTL"
    -- readFaultMultiCtl
    print "READING LANDSLIDE.CTL"
    readLandslideCtl