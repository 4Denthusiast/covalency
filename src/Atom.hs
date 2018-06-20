{-# LANGUAGE DataKinds      #-}
{-# LANGUAGE KindSignatures #-}
module Atom(
    Atom(..),
    OrbitalLabel,
    AtomLabel,
    Atoms,
    emptyAtom,
    evalAtomOrb,
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
    atomPotential :: Potential n
    atomKineticTerm :: M.Map OrbitalLabel (Linear (Complex Float) OrbitalLabel)
}
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

emptyAtom :: KnownNat n => [Float] -> Atom n
emptyAtom xs = Atom
    xs
    (M.fromList $ map (\i -> (i,normalize $ single $ centralGaussian (exp $ fromIntegral i) 1)) [-4..4])
    sphericalHarmonicPotential
    (M.fromList $ (4,Linear [(-0.5,3),(1,4)]) : map (\i -> (i,Linear [(-1,i),(1,i+1)])) [-4..3])

evalAtomOrb :: Atom n -> OrbitalLabel -> [Float] -> Complex Float
evalAtomOrb (Atom ax os _ _) ol xs = evaluates o xs'
    where o   = os M.! ol
          xs' = zipWith (-) xs ax
