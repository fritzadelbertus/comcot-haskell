
data Test = Test {
    hnCurr:: [[Double]]
}

test = Test {hnCurr = [[1,2,3,4,5],[6,7,8,9,0]]}

prevMomentNC =  ((!!) .) . ((!!) . hnCurr)

main:: IO()
main = do
    print $ prevMomentNC test 1 2