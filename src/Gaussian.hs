module Gaussian(
    Gaussian(..),
    evaluate,
    integral,
    convolve,
    multiply,
    shiftGauss,
    scaleGauss
) where

import Data.Complex

data Gaussian = Gaussian [Float] Float (Complex Float)

evaluate :: Gaussian -> [Float] -> Complex Float
evaluate (Gaussian xs c a) ys = (a*) $ pure $ exp (-norm2 (zipWith (-) xs ys) / c^2) / (sqrt pi ^ length xs)

integral :: Gaussian -> Complex Float
integral (Gaussian xs c a) = a * pure (c ^ length xs)

convolve :: Gaussian -> Gaussian -> Gaussian
convolve (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian (zipWith (+) xs xs') (c + c') (a * a')

multiply :: Gaussian -> Gaussian -> Gaussian
multiply (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian ys d (a * a' * pure (m * sqrt pi ^ length xs))
    where d  = 1/(1/c + 1/c')
          ys = zipWith (\x x' -> d * (x/c + x'/c')) xs xs'
          m  = exp (norm2 xs / c + norm2 xs' / c' - norm2 ys / d)

shiftGauss :: [Float] -> Gaussian -> Gaussian
shiftGauss dxs (Gaussian xs c a) = Gaussian (zipWith (+) dxs xs) c a

scaleGauss :: Complex Float -> Gaussian -> Gaussian
scaleGauss a' (Gaussian xs c a) = Gaussian xs c (a*a')

norm2 :: Num n => [n] -> n
norm2 = sum . map (^2)
