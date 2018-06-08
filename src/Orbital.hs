module Orbital(
    Orbital(..),
    evalOrbital
) where

import Linear
import Gaussian
import Atom

import qualified Data.Map as M
import Data.Complex

type Orbital = Linear (Complex Float) (AtomLabel, OrbitalLabel)

evalOrbital :: Atoms n -> Orbital -> [Float] -> Complex Float
evalOrbital as o xs = flatten $ lmap (\(al, ol) -> evalAtomOrb (as M.! al) ol xs) o
