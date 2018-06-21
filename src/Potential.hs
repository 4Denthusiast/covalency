{-# LANGUAGE TupleSections #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}
module Potential(
    Potential,
    sphericalHarmonicFunction,
    sphericalHarmonicPotential,
) where

import Linear
import Gaussian

import Data.Complex
import GHC.TypeLits

-- a distribution
type Potential (n::Nat) = Gaussians n -> Cplx

-- The range of scales over which this function may have a significant value and variation.
sizeRange :: Gaussians n -> (Float, Float)
sizeRange (Linear gs') = (minimum $ map inner gs, maximum $ map outer gs)
    where gs = map snd gs'
          inner (Gaussian xs c _) = max (0.01*sqrt c) (sqrt (norm2 xs) - 3*sqrt c)
          outer (Gaussian xs c _) = sqrt (norm2 xs) + 3*sqrt c

-- A linear combination of gaussians that approximates x^(2-n), the usual spherically-symmetric harmonic function, within the given range.
-- TODO: add dimension-2 stuff. (that's complicated though)
sphericalHarmonicFunction :: forall n. KnownNat n => (Float, Float) -> Gaussians n
sphericalHarmonicFunction (l', h') = Linear (map ((1,).g) [0..k])
    where d = 2 - natVal @n Proxy
          d' = div d 2
          l = l'*l' / 2
          h = h'*h' * case d of { -1 -> 20.0; -2 -> 1.0; -3 -> 3.0; _ -> error "Unsupported dimension"}
          k = ceiling $ log (h/l)
          g i = let x = l * exp (fromIntegral i) in centralGaussian (1/(1/x - 1/(l*exp(fromIntegral k+0.5)))) ((c i * s x) :+ 0)
          s :: Float -> Float
          s x = if even d then x ^^ d' / fac (-1-fromIntegral d') else x ^^ d' * sqrt x / fac (-1.5-fromIntegral d')
          c i = if i == 0 then 1 / (1 - exp (-1)) else 1 -- Correction to give the region around 0 the correct integral, even if it has the wrong shape.
          fac 0 = 1
          fac (-0.5) = sqrt pi
          fac x = x * fac (x-1)

sphericalHarmonicPotential :: KnownNat n => Potential n
sphericalHarmonicPotential gs = dot gs (sphericalHarmonicFunction $ sizeRange gs)
