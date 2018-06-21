{-# LANGUAGE DataKinds      #-}
{-# LANGUAGE KindSignatures #-}
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

import Data.Complex
import qualified Data.Map as M
import GHC.TypeLits

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
emptyAtom xs = Atom
        xs
        (M.fromList $ map (\i -> (i,normalize $ return $ centralGaussian (e i) 1)) [-4..4])
        (negate . sphericalHarmonicPotential)
        (M.fromList $ (4,mempty) : map (\i -> (i,Linear [(e (-2*i),i),(-e (-2*i),i+1)])) [-4..3])
    where e :: Floating a => Int -> a
          e = (exp . (0.9*) . fromIntegral)

evalAtomOrb :: Atom n -> OrbitalLabel -> [Float] -> Cplx
evalAtomOrb (Atom ax os _ _) ol xs = evaluates o xs'
    where o   = os M.! ol
          xs' = zipWith (-) xs ax

atomOrbitalsGlobal :: Atom n -> M.Map OrbitalLabel (Gaussians n)
atomOrbitalsGlobal at = (shiftGauss (atomPos at) <$>) <$> (atomOrbitals at)

atomPotentialGlobal :: Atom n -> Gaussians n -> Cplx
atomPotentialGlobal at = atomPotential at . (shiftGauss (negate <$> atomPos at) <$>)
