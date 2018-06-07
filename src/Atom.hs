{-# LANGUAGE DataKinds      #-}
{-# LANGUAGE KindSignatures #-}
module Atom(
    Atom(..),
    OrbitalLabel,
    AtomLabel,
    Atoms,
    emptyAtom
) where

import Gaussian

import Data.Complex
import qualified Data.Map as M
import GHC.TypeLits

type AtomLabel = String
type OrbitalLabel = ()
data Atom (n::Nat) = Atom [Float] (M.Map OrbitalLabel (Gaussian n))
type Atoms (n::Nat) = M.Map AtomLabel (Atom n)

emptyAtom :: KnownNat n => [Float] -> Int -> Atom n
emptyAtom xs ix = Atom xs (M.singleton () (centralGaussian 0.5 (cis $ fromIntegral ix)))
