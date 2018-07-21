{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Geometry(
    nuclearForces,
    bondEnergy,
) where

import Linear
import Gaussian
import Atom
import Orbital

import Data.Map (Map)
import qualified Data.Map as M
import Data.List
import GHC.TypeLits
import Debug.Trace

nuclearForces :: forall n. UsableDimension n => Atoms n -> [Orbital] -> Map AtomLabel [Rl]
nuclearForces ats orbs = M.mapWithKey (\al at -> zipWith (+) (electronForce at) (internuclearForces al at)) ats
    where electronDensity = reduceGaussians $ mconcat $ map (gsSquare . evalOrbital ats . snd . snd) orbs
          gsSquare gs = reduceGaussians $ multiply <$> gs <*> gs
          derivatives = map (flip (fmap . differentiate) electronDensity) [0..fromIntegral $ natVal @n Proxy-1]
          electronForce at = map (negate . atomPotentialGlobal at) derivatives
          internuclearForces al at = foldr (zipWith (+)) (genericReplicate (natVal @n Proxy) 0) $ M.mapWithKey (inForce al at) ats
          inForce al at al' at' = if al == al' then repeat 0 else (fromIntegral (atomicNumber at * atomicNumber at') *) <$> coulumbDiff (zipWith (-) (atomPos at) (atomPos at'))
          coulumbDiff xs = map (/ sqrt (norm2 xs) ^ natVal @n Proxy) xs

bondEnergy :: forall n. KnownNat n => Atoms n -> Integrals -> [Orbital] -> Rl
bondEnergy ats ints@(ov,nh,fei) orbs = sum (individualEnergy <$> M.keys ats) - totalEnergy ats ints orbs
    where individualEnergy al = energyWithNh (atomicNumber $ ats M.! al) $ nuclearHamiltonian (trimOrbitals ats) (Just al)
          energyWithNh z nh' = totalElectronEnergy ats (ov,nh',fei) $ hartreeFockSolution ats (ov,nh',fei) z
