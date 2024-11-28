import Numeric (showFFloat)


writeListsToFile :: FilePath -> [Double] -> [Double] -> [Double] -> IO ()
writeListsToFile filePath listA listB listC = do
    let combined = zip3 listA listB listC  -- Combine the lists by index
        formatLine (a, b, c) = unwords [show a, show b, showFFloat (Just 4) (c/100) ""]  -- Convert to space-separated string
        linesToWrite = map formatLine combined
    writeFile filePath (unlines linesToWrite)
readBathData:: ([String]-> String) -> String -> [Double]
readBathData index content = map (read . index . words) (lines content) 

readBathDataZ:: String -> [Double]
readBathDataZ = readBathData (!! 2)

main:: IO ()
main = do
    content <- readFile "depthy.xyz"
    let topoZ = readBathDataZ content
    let xCoor = concat [ [0..99] | x <-[0..99]]
    let yCoor = concat [replicate 100 x | x <-[0..99]]
    writeListsToFile "depthy2.xyz" xCoor yCoor $ map abs topoZ
    print $ showFFloat (Just 4) 0.23141441 ""