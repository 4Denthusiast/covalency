module Util(
    merge,
    swap,
    traceTime,
) where

import Data.Function
import Data.Time.Clock.System
import Data.Time.Clock.TAI
import Debug.Trace
import System.IO.Unsafe

merge :: Ord a => [a] -> [a] -> [a]
merge [] ys = ys
merge xs [] = xs
merge (x:xs) (y:ys) = if x <= y then x : merge xs (y:ys) else y : merge (x:xs) ys

swap :: (a,b) -> (b,a)
swap (x,y) = (y,x)

-- Print to console the amount of time it takes to evaluate x to HNF.
traceTime :: String -> a -> a
traceTime s x = unsafePerformIO $ do
    t0 <- getSystemTime
    t1 <- seq x $ getSystemTime
    let ts = s ++ show (on diffAbsoluteTime systemToTAITime t1 t0)
    traceIO ts
    return x
