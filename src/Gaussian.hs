{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving  #-}
module Gaussian(
    Gaussian(..),
    centralGaussian,
    evaluate,
    integral,
    convolve,
    multiply,
    shiftGauss,
    scaleGauss
) where

import Data.Complex
import Data.List
import GHC.TypeLits --(natVal, Nat, KnownNat)

data Gaussian (n::Nat) = Gaussian [Float] Float (Complex Float) deriving (Eq, Ord, Show)

centralGaussian :: forall n . KnownNat n => Float -> Complex Float -> Gaussian n
centralGaussian = Gaussian (genericReplicate (natVal @n Proxy) 0)

evaluate :: Gaussian n -> [Float] -> Complex Float
evaluate (Gaussian xs c a) ys = (a*) $ real $ exp (-norm2 (zipWith' (-) xs ys) / c)

integral :: forall n . KnownNat n => Gaussian n -> Complex Float
integral (Gaussian xs c a) = a * real (sqrt (c*pi) ^ natVal @n Proxy)

convolve :: forall n . KnownNat n => Gaussian n -> Gaussian n -> Gaussian n
convolve (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian (zipWith' (+) xs xs') (c + c') (a * a' * real (sqrt (pi/(1/c+1/c')) ^ natVal @n Proxy))

multiply :: Gaussian n -> Gaussian n -> Gaussian n
multiply (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian ys d (a * a' * real m)
    where d  = 1/(1/c + 1/c')
          ys = zipWith' (\x x' -> d * (x/c + x'/c')) xs xs'
          m  = exp (norm2 ys / d - norm2 xs / c - norm2 xs' / c')

shiftGauss :: [Float] -> Gaussian n -> Gaussian n
shiftGauss dxs (Gaussian xs c a) = Gaussian (zipWith' (+) dxs xs) c a

scaleGauss :: Complex Float -> Gaussian n -> Gaussian n
scaleGauss a' (Gaussian xs c a) = Gaussian xs c (a*a')

norm2 :: Num n => [n] -> n
norm2 = sum . map (^2)

real :: Num a => a -> Complex a
real = (:+0)

data Proxy (n :: Nat) = Proxy

--Zips lists and returns an error if the lengths don't match.
zipWith' :: (a->b->c) -> [a] -> [b] -> [c]
zipWith' f (x:xs) (y:ys) = f x y : zipWith' f xs ys
zipWith' _ [] [] = []
zipWith' _ _ _ = error "Non-matching list lengths."

deriving instance Ord a => Ord (Complex a)
