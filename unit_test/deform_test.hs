import Deform


main :: IO ()
main = do
    print $ stereoProjection 10.5 75.2 9.0 76.5 -- Returns 42772.2109, 144520.000 in Fortran

