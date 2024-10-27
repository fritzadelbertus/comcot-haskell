module Helper (removeLeadingSpaces) where

removeLeadingSpaces:: String -> String
removeLeadingSpaces = dropWhile (== ' ')