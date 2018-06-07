module Orbital(
    Orbital(..)
) where

import Gaussian
import Atom

import qualified Data.Map as M
import Data.Complex

data Orbital = Orbital [(AtomLabel, OrbitalLabel, Complex Float)]
