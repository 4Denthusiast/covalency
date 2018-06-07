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

type AtomLabel = String
type OrbitalLabel = ()
data Atom = Atom Float Float (M.Map OrbitalLabel Gaussian)
type Atoms = M.Map AtomLabel Atom

emptyAtom :: Float -> Float -> Int -> Atom
emptyAtom x y n = Atom x y (M.singleton () (Gaussian [0,0] 0.5 (2*(cis $ fromIntegral n))))
