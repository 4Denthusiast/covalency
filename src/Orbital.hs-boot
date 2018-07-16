{-# LANGUAGE FlexibleContexts #-}
module Orbital(
    Matrix,
    matTimes,
    invert,
    doInvert,
    tabulate,
) where

import Linear

import qualified Data.Map as M
import Data.Maybe
import Data.Bifunctor
import Data.Monoid
import Control.Applicative
import GHC.TypeLits
import GHC.Stack

type Matrix a = M.Map a (Linear Rl a)

overlaps :: (InnerProduct Rl v, Ord a) => [(a,v)] -> Matrix a

matTimes :: (HasCallStack, Ord a) => Matrix a -> Linear Rl a -> Linear Rl a

invert :: (Ord a) => Matrix a -> Maybe (Matrix a)

doInvert :: (Ord a) => Matrix a -> Matrix a

tabulate :: (Ord k) => [k] -> (k -> a) -> M.Map k a

swap :: (x,y) -> (y,x)
