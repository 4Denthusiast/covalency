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

import Data.Complex
import qualified Data.Map as M
import GHC.TypeLits

type AtomLabel = String
type OrbitalLabel = ()
data Atom (n::Nat) = Atom [Float] (M.Map OrbitalLabel (Linear (Complex Float) (Gaussian n)))
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

emptyAtom :: KnownNat n => [Float] -> Int -> Atom n
emptyAtom xs ix = Atom xs (M.singleton () (single $ centralGaussian 0.5 (cis $ fromIntegral ix)))

evalAtomOrb :: Atom n -> OrbitalLabel -> [Float] -> Complex Float
evalAtomOrb (Atom ax os) ol xs = flatten $ lmap (flip evaluate xs') o
    where o   = os M.! ol
          xs' = zipWith (-) xs ax
