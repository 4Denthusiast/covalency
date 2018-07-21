{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving  #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TupleSections #-}
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
    reverseGauss,
    scaleGauss,
    laplacian,
    differentiate,
    reduceGaussians,
    norm2,
    Proxy(..),
) where

import Linear
import {-# SOURCE #-} Polynomial (Polynomial)
import {-# SOURCE #-} qualified Polynomial as P

import Data.List
import Data.Function
import GHC.TypeLits --(natVal, Nat, KnownNat)

data Gaussian (n::Nat) = Gaussian [Rl] Rl (Polynomial n Rl) deriving (Eq, Show, Ord)
type Gaussians (n::Nat) = Linear Rl (Gaussian n)

centralGaussian :: forall n . KnownNat n => Rl -> (Polynomial n Rl) -> Gaussian n
centralGaussian = Gaussian (genericReplicate (natVal @n Proxy) 0)

centralSphereGaussian :: KnownNat n => Rl -> Rl -> Gaussian n
centralSphereGaussian c a = centralGaussian c (P.constant a)

evaluate :: Gaussian n -> [Rl] -> Rl
evaluate (Gaussian xs c a) ys = exp (-norm2 ys' / c) * P.evaluate' ys' a
    where ys' = zipWith' (-) ys xs

evaluates :: Gaussians n -> [Rl] -> Rl
evaluates gs xs = flatten $ flip evaluate xs <$> gs

integral :: forall n . KnownNat n => Gaussian n -> Rl
integral (Gaussian xs c a) = (sqrt (c*pi) ^ natVal @n Proxy) * ai
    where ai = P.monomialSum' monCoeff a
          monCoeff es = if any odd es then 0 else product $ map powerCoeff es
          powerCoeff 0 = 1
          powerCoeff n = c / 2 * (fromIntegral n - 1) * powerCoeff (n-2) --Integrate by parts

convolve :: forall n . KnownNat n => Gaussian n -> Gaussian n -> Gaussian n
convolve (Gaussian xs c a) (Gaussian xs' c' a') = Gaussian (zipWith' (+) xs xs') (c + c') (((sqrt (pi/(1/c+1/c'))) ^ natVal @n Proxy) *~ a'')
    where aShift :: Polynomial n (Polynomial n Rl)
          aShift = P.evaluate (xpky (c/(c+c'))) a * P.evaluate (map negate $ xpky (-c'/(c+c'))) a'
          xpky :: Rl -> [Polynomial n (Polynomial n Rl)]
          xpky k = map (\i -> P.variable i + k *~ (P.constant (P.variable i))) [0..]
          a''    = P.monomialSum' (\es -> if any odd es then 0 else product $ map powerCoeff es) aShift
          powerCoeff 0 = 1
          powerCoeff n = c*c'/(c+c') / 2 * (fromIntegral n - 1) * powerCoeff (n-2)

multiply :: KnownNat n => Gaussian n -> Gaussian n -> Gaussian n
multiply (Gaussian xs c a) (Gaussian xs' c' a') = if xs == xs' then Gaussian xs d (a*a') else Gaussian ys d a''
    where d  = 1/(1/c + 1/c')
          ys = zipWith' (\x x' -> d * (x/c + x'/c')) xs xs'
          m  = exp (norm2 ys / d - norm2 xs / c - norm2 xs' / c')
          a'' = P.constant m * shiftPoly (zipWith (-) ys xs) a * shiftPoly (zipWith (-) ys xs') a'

shiftGauss :: KnownNat n => [Rl] -> Gaussian n -> Gaussian n
shiftGauss dxs (Gaussian xs c a) = Gaussian (zipWith' (+) dxs xs) c a

reverseGauss :: KnownNat n => Gaussian n -> Gaussian n
reverseGauss (Gaussian xs c a) = Gaussian (map negate xs) c (P.evaluate (map (negate . P.variable) [0..]) a)

shiftPoly :: (Num a, Eq a, KnownNat n) => [a] -> Polynomial n a -> Polynomial n a
shiftPoly as p = P.evaluate (zipWith (\a i -> P.variable i + P.constant a) as [0..]) p

scaleGauss :: KnownNat n => Rl -> Gaussian n -> Gaussian n
scaleGauss a' (Gaussian xs c a) = Gaussian xs c (a' *~ a)

laplacian :: forall n. KnownNat n => Gaussian n -> Gaussian n
laplacian (Gaussian xs c a) = Gaussian xs c $ sum $ map (\i -> diffPoly c i (diffPoly c i a)) [0..fromIntegral $ natVal @n Proxy - 1]

differentiate :: KnownNat n => Int -> Gaussian n -> Gaussian n
differentiate i (Gaussian xs c a) = Gaussian xs c $ diffPoly c i a

diffPoly :: KnownNat n => Rl -> Int -> Polynomial n Rl -> Polynomial n Rl
diffPoly c i p = P.monomialSum (dMon . splitAt i) p - (2 / c) *~ p * P.variable i
    where dMon (h,e:t) = fromIntegral e * P.monomial (h ++ (e-1) : t)

reduceGaussians :: forall n. KnownNat n => Gaussians n -> Gaussians n
reduceGaussians = Linear . map mergeGroup . toGroups . reduce
    where toGroups (Linear xs) = groupBy (on (==) (\(_,Gaussian p c _) -> (p,c))) xs
          mergeGroup :: [(Rl,Gaussian n)] -> (Rl,Gaussian n)
          mergeGroup [g] = g
          mergeGroup gs@((_,Gaussian xs c _):_) = (1,) $ Gaussian xs c $ sum $ map (\(a,Gaussian _ _ p) -> a *~ p) gs

norm2 :: Num n => [n] -> n
norm2 = sum . map (^2)

data Proxy (n :: Nat) = Proxy

--Zips lists and returns an error if the lengths don't match.
zipWith' :: (a->b->c) -> [a] -> [b] -> [c]
zipWith' f (x:xs) (y:ys) = f x y : zipWith' f xs ys
zipWith' _ [] [] = []
zipWith' _ _ _ = error "Non-matching list lengths."

instance KnownNat n => InnerProduct Rl (Gaussian n) where
    dot g g' = integral $ multiply g g'
