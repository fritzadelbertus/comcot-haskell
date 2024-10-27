
nonCommentLines::  String -> Int
nonCommentLines c = length $ filter notAComment $ lines c
    where 
        notAComment line = not (null line) && head line /= '!'

fortranAdresses:: [String]
fortranAdresses = map addPath [
        "all_grids",
        "boundaries",
        "comcot",
        "deform",
        "dispersion",
        "hotstart",
        "initialization",
        "landslide",
        "mass",
        "moment",
        "output",
        "type_module",
        "wavemaker"
    ]
    where 
        addPath adr = "comcot-fortran/" ++ adr ++ ".f90"

getNonCommentLines:: Int
getNonCommentLines = sum $ map (nonCommentLines . read) fortranAdresses
            

main :: IO ()
main = do
    content <- readFile $ fortranAdresses !! 0
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 0
    content <- readFile $ fortranAdresses !! 1
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 1
    content <- readFile $ fortranAdresses !! 2
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 2
    content <- readFile $ fortranAdresses !! 3
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 3
    content <- readFile $ fortranAdresses !! 4
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 4
    content <- readFile $ fortranAdresses !! 5
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 5
    content <- readFile $ fortranAdresses !! 6
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 6
    content <- readFile $ fortranAdresses !! 7
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 7
    content <- readFile $ fortranAdresses !! 8
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 8
    content <- readFile $ fortranAdresses !! 9
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 9
    content <- readFile $ fortranAdresses !! 10
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 10
    content <- readFile $ fortranAdresses !! 11
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 11
    content <- readFile $ fortranAdresses !! 12
    print $ "There are " ++ show (nonCommentLines content) ++ " lines in " ++ fortranAdresses !! 12