{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving  #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Gaussian(
    Gaussian(..),
    Gaussians,
    centralGaussian,
    centralSphereGaussian,
    evaluate,
    evaluates,
    integral,
    convolve,
    multiply,
    shiftGauss,
    scaleGauss,
    laplacian,
    norm2,
    Proxy(..),
) where

import Linear
import Polynomial (Polynomial)
import qualified Polynomial as P

import Data.Complex
import Data.List
import GHC.TypeLits --(natVal, Nat, KnownNat)

data Gaussian (n::Nat) = Gaussian [Float] Float (Polynomial n Cplx) deriving (Eq, Show)
type Gaussians (n::Nat) = Linear Cplx (Gaussian n)

centralGaussian :: forall n . KnownNat n => Float -> (Polynomial n Cplx) -> Gaussian n
centralGaussian = Gaussian (genericReplicate (natVal @n Proxy) 0)

centralSphereGaussian :: KnownNat n => Float -> Cplx -> Gaussian n
centralSphereGaussian c a = centralGaussian c (P.constant a)

evaluate :: Gaussian n -> [Float] -> Cplx
evaluate (Gaussian xs c a) ys = exp (-norm2 ys' / c) *~ P.evaluate' ys' a
    where ys' = zipWith' (-) ys xs

evaluates :: Gaussians n -> [Float] -> Cplx
evaluates gs xs = flatten $ flip evaluate xs <$> gs

integral :: forall n . KnownNat n => Gaussian n -> Cplx
integral (Gaussian xs c a) = (sqrt (c*pi) ^ natVal @n Proxy) *~ ai
    where ai = P.monomialSum' monCoeff a
          monCoeff es = if any odd es then 0 else product $ map powerCoeff es
          powerCoeff 0 = 1
          powerCoeff n = c / 2 * (fromIntegral n - 1) * powerCoeff (n-2) --Integrate by parts

convolve :: forall n . KnownNat n => Gaussian n -> Gaussian n -> Gaussian n
convolve (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian (zipWith' (+) xs xs') (c + c') (((sqrt (pi/(1/c+1/c'))) ^ natVal @n Proxy) *~ a'')
    where aShift :: Polynomial n (Polynomial n Cplx)
          aShift = P.evaluate (xpky (c/(c+c'))) a * P.evaluate (map negate $ xpky (-c'/(c+c'))) a'
          xpky :: Float -> [Polynomial n (Polynomial n Cplx)]
          xpky k = map (\i -> P.variable i + k *~ (P.constant (P.variable i))) [0..]
          a''    = P.monomialSum' (\es -> if any odd es then 0 else product $ map powerCoeff es) aShift
          powerCoeff 0 = 1
          powerCoeff n = c*c'/(c+c') / 2 * (fromIntegral n - 1) * powerCoeff (n-2)

multiply :: KnownNat n => Gaussian n -> Gaussian n -> Gaussian n
multiply (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian ys d a''
    where d  = 1/(1/c + 1/c')
          ys = zipWith' (\x x' -> d * (x/c + x'/c')) xs xs'
          m  = exp (norm2 ys / d - norm2 xs / c - norm2 xs' / c')
          diff y x = real (y-x)
          a'' = P.constant (real m) * shiftPoly (zipWith diff ys xs) a * shiftPoly (zipWith diff ys xs') a'

shiftGauss :: KnownNat n => [Float] -> Gaussian n -> Gaussian n
shiftGauss dxs (Gaussian xs c a) = Gaussian (zipWith' (+) dxs xs) c a

shiftPoly :: (Num a, Eq a, KnownNat n) => [a] -> Polynomial n a -> Polynomial n a
shiftPoly as p = P.evaluate (zipWith (\a i -> P.variable i + P.constant a) as [0..]) p

scaleGauss :: KnownNat n => Cplx -> Gaussian n -> Gaussian n
scaleGauss a' (Gaussian xs c a) = Gaussian xs c (a * P.constant a')

laplacian :: forall n. KnownNat n => Gaussian n -> Gaussian n
{-
laplacian (Gaussian xs c a) = Gaussian xs c $ P.monomialSum monLap a
    where monLap es = sum $ zipWith3 dir (inits es) es (tail $ tails es)
          dir :: [Int] -> Int -> [Int] -> Polynomial n Cplx
          dir h e t =
              e*(e-1)                  *~ mon (h++(e-2):t)
              - (2+fromIntegral e*4)/c *~ mon (h++ e   :t)
              + 4/(c*c)                *~ mon (h++(e+2):t)
          mon es = mon' 0 es
          mon' n (e:es) = if e < 0 then 0 else (P.variable n ^ e) * mon' (n+1) es
          mon' _ [] = 1
-}
laplacian (Gaussian xs c a) = Gaussian xs c $ sum $ map (\i -> d i (d i a)) [0..fromIntegral $ natVal @n Proxy - 1]
    where d :: Int -> Polynomial n Cplx -> Polynomial n Cplx
          d i p = P.monomialSum (dMon . splitAt i) p - (2 / c) *~ p * P.variable i
          dMon (h,e:t) = fromIntegral e * mon (h ++ (e-1) : t)
          mon es = mon' 0 es
          mon' n (e:es) = if e < 0 then 0 else (P.variable n ^ e) * mon' (n+1) es
          mon' _ [] = 1

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

instance Semilinear (Gaussian n) where
    conj (Gaussian xs c a) = Gaussian xs c (conj a)

instance KnownNat n => InnerProduct Cplx (Gaussian n) where
    dot g g' = integral $ multiply (conj g) g'
