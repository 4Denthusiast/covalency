{-# LANGUAGE FlexibleContexts #-}
module Orbital(
    Matrix,
    invert,
    doInvert,
    tabulate,
    swap,
) where

import Linear

import qualified Data.Map as M
import Data.Maybe
import Data.Complex
import Data.Bifunctor
import Data.Monoid
import Control.Applicative
import GHC.TypeLits

type Matrix a = M.Map a (Linear Cplx a)

overlaps :: (InnerProduct Cplx v, Ord a) => [(a,v)] -> Matrix a

invert :: (Ord a) => Matrix a -> Maybe (Matrix a)

doInvert :: (Ord a) => Matrix a -> Matrix a

tabulate :: (Ord k) => [k] -> (k -> a) -> M.Map k a

swap :: (x,y) -> (y,x)