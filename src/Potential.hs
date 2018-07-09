{-# LANGUAGE TupleSections #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}
module Potential(
    Potential,
    coulumbFunction,
    coulumbPotential,
) where

import Linear
import Gaussian
import {-#SOURCE#-} Polynomial

import Data.Complex
import GHC.TypeLits

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
coulumbFunction :: forall n. KnownNat n => (Rl, Rl) -> Gaussians n
coulumbFunction (l', h') = Linear (map ((1,).g) [0..k])
    where d = 2 - natVal @n Proxy
          d' = div d 2
          l = l'*l' / 2
          h = h'*h' * case d of { -1 -> 100.0; -2 -> 1.0; -3 -> 3.0; _ -> error "Unsupported dimension"}
          k = ceiling $ log (h/l)
          g i = let x = l * exp (fromIntegral i) in centralSphereGaussian (1/(1/x - 1/(l*exp(fromIntegral k+0.5)))) ((c i * s x) :+ 0)
          s :: Rl -> Rl
          s x = if even d then x ^^ d' / fac (-1-fromIntegral d') else x ^^ d' * sqrt x / fac (-1.5-fromIntegral d')
          c i = if i == 0 then 1 / (1 - exp (-1)) else 1 -- Correction to give the region around 0 the correct integral, even if it has the wrong shape.
          fac 0 = 1
          fac (-0.5) = sqrt pi
          fac x = x * fac (x-1)

coulumbPotential :: KnownNat n => Potential n
coulumbPotential gs = if central then flatten (centralCoulumb <$> gs) else dot gs (coulumbFunction $ sizeRange gs)
    where (Linear gsl) = gs
          central = all ((\(Gaussian xs _ _) -> all (==0) xs) . snd) gsl
          centralCoulumb (Gaussian _ c a) = monomialSum' (\m -> if any odd m then 0 else sphereIntegral m * radialIntegral c (sum m)) a
          sphereIntegral m = 2* product (map (gammaHalf . (1+)) m) / gammaHalf (sum (map (1+) m))
          gammaHalf 1 = sqrt pi
          gammaHalf 2 = 1
          gammaHalf n = gammaHalf (n-2) * fromIntegral (n-2) / 2
          radialIntegral c n = (c ^ (1+div n 2))/2 * gammaHalf (n+2)
