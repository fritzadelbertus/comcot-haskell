import Text.Parsec
import Text.Parsec.String
import Data.List (intercalate)
import Control.Monad (liftM)

-- Earth's radius (in kilometers)
radius :: Double
radius = 6371.0

-- Convert degrees to radians
degToRad :: Double -> Double
degToRad degrees = degrees * pi / 180.0

-- Convert spherical coordinates (Longitude, Latitude, Depth) to Cartesian (X, Y, Z)
sphericalToCartesian :: Double -> Double -> Double -> (Double, Double, Double)
sphericalToCartesian lon lat depth = 
    ( x, y, depth )
  where
    -- Convert lon and lat from degrees to radians
    lonRad = degToRad lon
    latRad = degToRad lat
    
    -- Cartesian coordinates
    x = radius * cos(latRad) * cos(lonRad)
    y = radius * cos(latRad) * sin(lonRad)

-- Parse a line in the XYZ file
parseLine :: Parser (Double, Double, Double)
parseLine = do
    lon <- liftM read (many1 (digit <|> char '.' <|> char '-'))
    space
    lat <- liftM read (many1 (digit <|> char '.' <|> char '-'))
    space
    depth <- liftM read (many1 (digit <|> char '.' <|> char '-'))
    return (lon, lat, depth)

-- Read XYZ file and parse all lines
parseXYZFile :: String -> Either ParseError [(Double, Double, Double)]
parseXYZFile input = parse (many parseLine) "" input

-- Convert the content of XYZ file from spherical to Cartesian
convertFile :: String -> Either ParseError String
convertFile input = do
    coordinates <- parseXYZFile input
    let cartesianCoords = map (\(lon, lat, depth) -> sphericalToCartesian lon lat depth) coordinates
    return (unlines (map showCartesian cartesianCoords))

-- Show the Cartesian coordinates as a string
showCartesian :: (Double, Double, Double) -> String
showCartesian (x, y, z) = intercalate " " (map show [x, y, z])

-- Main function to read the file and print the converted Cartesian coordinates
main :: IO ()
main = do
    -- Read XYZ file as a string
    let filePath = "selatsunda.xyz"
    content <- readFile filePath
    
    -- Convert the content from spherical to Cartesian
    case convertFile content of
        Left err -> print err
        Right result -> do
            putStrLn "Converted Cartesian coordinates:"
            putStrLn result