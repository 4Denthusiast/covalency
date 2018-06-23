{-# LANGUAGE DataKinds      #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleContexts #-}
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
import Gaussian
import Potential
import {-# SOURCE #-} Orbital

import Data.Complex
import qualified Data.Map as M
import GHC.TypeLits
import Debug.Trace
import qualified Polynomial as P

type AtomLabel = String
type OrbitalLabel = Int
data Atom (n::Nat) = Atom{
    atomPos :: [Float],
    atomOrbitals :: M.Map OrbitalLabel (Gaussians n),
    atomPotential :: Potential n,
    atomKineticTerm :: M.Map OrbitalLabel (Linear Cplx OrbitalLabel)
}
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

emptyAtom :: KnownNat n => [Float] -> Atom n
emptyAtom xs = addKinetic $ Atom
        xs
        (M.fromList $ map (\i -> (i,normalize @Cplx $ return $ centralGaussian (e i) 1)) [-4..4])
        (((-2)*) . sphericalHarmonicPotential)
        undefined
    where e :: Floating a => Int -> a
          e = (exp . (1.5*) . fromIntegral)

evalAtomOrb :: Atom n -> OrbitalLabel -> [Float] -> Cplx
evalAtomOrb (Atom ax os _ _) ol xs = evaluates o xs'
    where o   = os M.! ol
          xs' = zipWith (-) xs ax

atomOrbitalsGlobal :: KnownNat n => Atom n -> M.Map OrbitalLabel (Gaussians n)
atomOrbitalsGlobal at = (shiftGauss (atomPos at) <$>) <$> (atomOrbitals at)

atomPotentialGlobal :: KnownNat n => Atom n -> Gaussians n -> Cplx
atomPotentialGlobal at = atomPotential at . (shiftGauss (negate <$> atomPos at) <$>)

addKinetic :: KnownNat n => Atom n -> Atom n
addKinetic at = seq (traceShowMatId $ kineticTest at) at{atomKineticTerm = (decompose . ((-0.5::Cplx)*~) . fmap laplacian) <$> atomOrbitals at}
    where decompose gs = reduce $ (invert allOverlaps M.!) =<< overlaps gs
          overlaps gs  = Linear (map swap $ M.toList $ dot gs <$> (atomOrbitals at))
          allOverlaps  = overlaps <$> atomOrbitals at

kineticTest :: KnownNat n => Atom n -> Matrix OrbitalLabel
kineticTest at = (\o -> Linear $ map swap $ M.toList $ (dot o . fmap laplacian) <$> atomOrbitals at) <$> atomOrbitals at

cheatPotential :: Potential n
cheatPotential gs = flatten $ cp <$> gs
    where cp (Gaussian xs c a)
              | norm2 xs < 0.001*c = (2*pi*c :+ 0)* (P.evaluate [0,0,0] a)^2
              | otherwise          = error (show xs)
