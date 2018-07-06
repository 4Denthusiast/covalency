{-# LANGUAGE DataKinds      #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Atom(
    Atom(..),
    OrbitalLabel,
    AtomLabel,
    Atoms,
    Spin(..),
    newAtom,
    changeZ,
    atomOrbitalsGlobal,
    atomPotentialGlobal,
) where

import Linear
import {-# SOURCE #-} qualified Polynomial as P
import Gaussian
import Potential
import {-# SOURCE #-} Orbital

import Data.Complex
import qualified Data.Map as M
import Control.Applicative
import GHC.TypeLits
import Debug.Trace

type AtomLabel = String
type L = Int
type M = Int
data Spin = Up | Down deriving (Eq, Ord)
instance Show Spin where{show Up = "↑"; show Down = "↓"}
type OrbitalLabel = (Int,L,M)
data Atom (n::Nat) = Atom{
    atomPos :: [Rl],
    atomicNumber :: Int,
    atomOrbitals :: M.Map OrbitalLabel (Gaussians n),
    atomPotential :: Potential n
}
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

newAtom :: forall n. KnownNat n => Int -> [Rl] -> Atom n
newAtom z xs = changeZ z $ Atom {
        atomPos = xs,
        atomOrbitals = (M.fromList $ liftA2 (\i (l,m,p) -> ((i,l,m),normalize @Cplx $ return $ centralGaussian (e i) p)) [-4..4] (concatMap (\l -> zipWith (l,,) [0..] $ P.sphericalHarmonicPolys l) [0,1]))
    }
    where e :: Floating a => Int -> a
          e = (exp . (1.5*) . fromIntegral)

changeZ :: KnownNat n => Int -> Atom n -> Atom n
changeZ z at = at{
        atomicNumber = z,
        atomPotential = ((fromIntegral z*) . negate . coulumbPotential)
    }

atomOrbitalsGlobal :: KnownNat n => Atom n -> M.Map OrbitalLabel (Gaussians n)
atomOrbitalsGlobal at = (shiftGauss (atomPos at) <$>) <$> (atomOrbitals at)

atomPotentialGlobal :: KnownNat n => Atom n -> Gaussians n -> Cplx
atomPotentialGlobal at = atomPotential at . (shiftGauss (negate <$> atomPos at) <$>)
