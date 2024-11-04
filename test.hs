data Paper = Paper {
    size:: Int,
    count:: Int
} deriving Show

a4Paper :: Paper
a4Paper = Paper {size = 25, count = 2}

a5Paper :: Paper
a5Paper = a4Paper {size = 20}

main:: IO()
main = do
    print a4Paper
    print a5Paper