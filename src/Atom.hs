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
    emptyAtom,
    evalAtomOrb,
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
type OrbitalLabel = (Int,L,M)
data Atom (n::Nat) = Atom{
    atomPos :: [Rl],
    atomOrbitals :: M.Map OrbitalLabel (Gaussians n),
    atomPotential :: Potential n,
    atomKineticTerm :: M.Map OrbitalLabel (Linear Cplx OrbitalLabel)
}
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

emptyAtom :: forall n. KnownNat n => [Rl] -> Atom n
emptyAtom xs = addKinetic $ Atom
        xs
        (M.fromList $ liftA2 (\i (l,m,p) -> ((i,l,m),normalize @Cplx $ return $ centralGaussian (e i) p)) [-3..2] (concatMap (\l -> zipWith (l,,) [0..] $ P.sphericalHarmonicPolys l) [0,1]))
        ((1*) . negate . coulumbPotential)
        undefined
    where e :: Floating a => Int -> a
          e = (exp . (1.5*) . fromIntegral)

evalAtomOrb :: Atom n -> OrbitalLabel -> [Rl] -> Cplx
evalAtomOrb (Atom ax os _ _) ol xs = evaluates o xs'
    where o   = os M.! ol
          xs' = zipWith (-) xs ax

atomOrbitalsGlobal :: KnownNat n => Atom n -> M.Map OrbitalLabel (Gaussians n)
atomOrbitalsGlobal at = (shiftGauss (atomPos at) <$>) <$> (atomOrbitals at)

atomPotentialGlobal :: KnownNat n => Atom n -> Gaussians n -> Cplx
atomPotentialGlobal at = atomPotential at . (shiftGauss (negate <$> atomPos at) <$>)

addKinetic :: KnownNat n => Atom n -> Atom n
addKinetic at = at{atomKineticTerm = (decompose . ((-0.5::Cplx)*~) . fmap laplacian) <$> atomOrbitals at}
    where decompose gs = reduce $ (doInvert allOverlaps M.!) =<< overlaps gs
          overlaps gs  = Linear (map swap $ M.toList $ dot gs <$> (atomOrbitals at))
          allOverlaps  = overlaps <$> atomOrbitals at
