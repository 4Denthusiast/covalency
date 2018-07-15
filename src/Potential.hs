{-# LANGUAGE TupleSections #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE FlexibleInstances #-}
module Potential(
    UsableDimension,
    Potential,
    coulumbFunction,
    coulumbPotential,
) where

import Linear
import Gaussian
import {-#SOURCE#-} Polynomial as P

import Data.Complex
import GHC.TypeLits
import Debug.Trace

-- a distribution
type Potential (n::Nat) = Gaussians n -> Cplx

-- The range of scales over which this function may have a significant value and variation.
sizeRange :: Gaussians n -> (Rl, Rl)
sizeRange (Linear gs') = (minimum $ map inner gs, maximum $ map outer gs)
    where gs = map snd gs'
          inner (Gaussian xs c _) = max (0.01*sqrt c) (sqrt (norm2 xs) - 3*sqrt c)
          outer (Gaussian xs c _) = sqrt (norm2 xs) + 3*sqrt c

-- A linear combination of gaussians that approximates x^(2-n), the usual spherically-symmetric harmonic function, within the given range.
-- TODO: add dimension-2 stuff. (that's complicated though)
coulumbFunction :: forall n. UsableDimension n => (Rl, Rl) -> Gaussians n
coulumbFunction (l', h') = Linear (map ((1,).g) [0..k])
    where d = 2 - natVal @n Proxy
          d' = div d 2
          l = l'*l' / 2
          h = h'*h' * coulumbApproximationThreshold @n
          k = ceiling $ log (h/l)
          g i = let x = l * exp (fromIntegral i) in centralSphereGaussian (1/(1/x - 1/(l*exp(fromIntegral k+0.5)))) ((c i * s x) :+ 0)
          s :: Rl -> Rl
          s x = if even d then x ^^ d' / fac (-1-fromIntegral d') else x ^^ d' * sqrt x / fac (-1.5-fromIntegral d')
          c i = if i == 0 then 1 / (1 - exp (-1)) else 1 -- Correction to give the region around 0 the correct integral, even if it has the wrong shape.
          fac 0 = 1
          fac (-0.5) = sqrt pi
          fac x = x * fac (x-1)

defaultCoulumbPotential :: UsableDimension n => Potential n
defaultCoulumbPotential gs = if central then flatten (centralCoulumb <$> gs) else dot gs (coulumbFunction $ sizeRange gs)
    where (Linear gsl) = gs
          central = all ((\(Gaussian xs _ _) -> all (==0) xs) . snd) gsl

centralCoulumb :: KnownNat n => Gaussian n -> Cplx
centralCoulumb (Gaussian _ c a) = monomialSum' (\m -> if any odd m then 0 else sphereIntegral m * radialIntegral c (sum m)) a
    where sphereIntegral m = 2* product (map (gammaHalf . (1+)) m) / gammaHalf (sum (map (1+) m))
          gammaHalf 1 = sqrt pi
          gammaHalf 2 = 1
          gammaHalf n = gammaHalf (n-2) * fromIntegral (n-2) / 2
          radialIntegral c n = (c ^ (1+div n 2))/2 * gammaHalf (n+2)

{- The value for (Gaussian p c 1) is exactly (pi c)^3/2 erf(|p|/c)/|p|. Values for higher L are computed
   by expressing the gaussian by derivatives of the basic gaussian, and taking the corresponding derivatives
   of the potential function. -}
coulumb3 :: Gaussian 3 -> Cplx
coulumb3 g@(Gaussian xs c a) = if norm2 xs < c*4e-6 then centralCoulumb g else monomialSum' (P.evaluate xs' . derivatives) a'
    where a' = hermitify c a
          xs' = [1/xsl, exp (- xsl*xsl / c), erf (xsl/sqrt c)] ++ xs --Add auxilliary variables for convenient manipulation.
          xsl = sqrt $ norm2 xs
          derivatives [x,y,z] = diff x 0 $ diff y 1 $ diff z 2 $ ((pi * c) ** 1.5) *~ variable 0 * variable 2
          diff 0 _ = id
          diff n i = diff (n-1) i . monomialSum (\es -> d0 i es + d1 i es + d2 i es + d3 i es)
          d0, d1, d2, d3 :: Int -> [Int] -> Polynomial 6 Rl
          d0 i es@(e:_:_:_) = -e *~ monomial (zipWith (+) es ([2,0,0] ++ δ i))
          d1 i es@(_:e:_:_) = e *~ (-2/c) *~ monomial (zipWith (+) es ([0,0,0] ++ δ i))
          d2 i es@(_:_:e:_) = e *~ (2/sqrt (pi * c)) *~ monomial (zipWith (+) es ([1,1,-1] ++ δ i))
          d3 i es@(_:_:_:es') = (es' !! i) *~ monomial (zipWith ({-ba'e-} -) es ([0,0,0] ++ δ i))
          δ i = case i of {0 -> [1,0,0]; 1 -> [0,1,0]; 2 -> [0,0,1]}
instance NamedDimensions 6 where
    dimName _ 0 = "r-"
    dimName _ 1 = "g"
    dimName _ 2 = "e"
    dimName _ 3 = "x"
    dimName _ 4 = "y"
    dimName _ 5 = "z"

-- Standard error function (approximation)
erf :: Rl -> Rl
erf x = sqrt $ (1-) $ exp $ -x^2 * (4/pi + a*x^2) / (1 + a*x^2)
    where a = 8*(pi-3) / (3*pi*(4-pi))
--erf x = if x > 3 then 1 - exp(-x*x) / (sqrt pi * x) else sum (zipWith (\k a -> x*(-x*x)^k * a) [0..35] erfCoeffs)

erfCoeffs :: [Rl]
erfCoeffs = map (\k -> 2 / sqrt pi / (fac k * (2*k+1))) [0..]
    where fac 0 = 1
          fac n = n * fac (n-1)

-- Change variables from x,y,... to d/dx,d/dy,... as applied to exp(-xs^2/c)
hermitify :: forall n. KnownNat n => Rl -> Polynomial n Cplx -> Polynomial n Cplx
hermitify c = hermitify' 0
    where hermitify' xs' 0 = xs'
          hermitify' xs' xs = let (dxs',dxs) = lead xs in hermitify' (xs' + dxs') (trimDeg (degree xs-1) $ xs - dxs)
          lead :: Polynomial n Cplx -> (Polynomial n Cplx, Polynomial n Cplx)
          lead xs = (c/2)^degree xs *~ monomialSum (\es -> if sum es == degree xs then (monomial es, herms (degree xs) es) else (0,0)) xs
          herms :: Int -> [Int] -> Polynomial n Cplx
          herms d = (1/sqrt c ^ d *~) . product . zipWith (\i e -> P.evaluate [1/sqrt c *~ variable i] (hermitianPoly e)) [0..]
          trimDeg :: Int -> Polynomial n Cplx -> Polynomial n Cplx
          trimDeg d = monomialSum (\es -> if sum es <= d then monomial es else 0) -- Remove the tiny bits left over due to numerical error.

hermitianPoly :: Int -> Polynomial 1 Rl
hermitianPoly 0 = 1
hermitianPoly n = monomialSum (\[e] -> (2::Rl) *~ monomial [e+1] - e *~ monomial [e-1]) $ hermitianPoly (n-1)

class KnownNat n => UsableDimension (n::Nat) where
    coulumbApproximationThreshold :: Rl
    coulumbPotential :: Potential n
    coulumbPotential = defaultCoulumbPotential

instance UsableDimension 3 where
    coulumbApproximationThreshold = 100
    coulumbPotential = flatten . fmap coulumb3

instance UsableDimension 4 where
    coulumbApproximationThreshold = 1
