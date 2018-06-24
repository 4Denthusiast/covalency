{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeApplications #-}
module Eigen(
    inverseIterants,
    rayleighIterate,
    allEigenvecs,
) where

import Linear
import Atom
import Orbital

import qualified Data.Map as M
import Data.Monoid
import Data.Complex
import Debug.Trace

-- This doesn't make much physical sense, but it's handy for an arbitrary normalization factor.
instance InnerProduct Cplx (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0

matVecTimes :: Ord a => Matrix a -> Linear Cplx a -> Linear Cplx a
matVecTimes m v = reduce $ v >>= (m M.!)

rayleighQuotient :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Cplx
rayleighQuotient m v = dot v (matVecTimes m v) / dot v v

-- Can't exclude linear combinations of eigenvectors with almost-equal values.
eigenvectorQuality :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Cplx
eigenvectorQuality m v = norm (v' <> (-(dot v v' / dot v v)::Cplx) *~ v)
    where v' = matVecTimes m v
          norm x = dot x x

-- (m - μ I)^-1
offsetInverse :: (Ord a) => Matrix a -> Cplx -> Matrix a
offsetInverse m μ = invert (M.mapWithKey ((reduce .) . (<>) . ((-μ) *~) . return) m)

inverseIterants :: (InnerProduct Cplx a, Ord a) => Matrix a -> Cplx -> Linear Cplx a -> [Linear Cplx a]
inverseIterants m μ v = iterate (normalize @Cplx . matVecTimes (offsetInverse m μ)) v

rayleighIterate :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Linear Cplx a
rayleighIterate m v = normalize @Cplx $ matVecTimes (offsetInverse m $ rayleighQuotient m v) v

--Assumes all eigenvectors to be orthogonal, which isn't actually the case here.
--In particular, the orthogonalisation step doesn't work correctly otherwise
allEigenvecs :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> [Linear Cplx a]
allEigenvecs m v0 = allEigenvecs' []
    where allEigenvecs' vs = let v = next vs $ search 5 vs v0 in v : allEigenvecs' (v:vs)
          search 0 vs v = trace (seq v "searched.") v
          search n vs v = search (n-1) vs $ orthogonalise vs $ matVecTimes (offsetInverse m (less $ traceShowId $ rayleighQuotient m v)) v
          less x = minimum [x-1, x*2, -x]
          next vs v
              | eigenvectorQuality m v < 0.0001 = normalize @Cplx v
              | otherwise = let v' = rayleighIterate m $ orthogonalise vs v in if wrong v' then v else next vs v'
          orthogonalise [] v = normalize @Cplx v
          orthogonalise (x:xs) v = orthogonalise xs $ reduce $ v <> (-dot v x::Cplx) *~ x
          wrong (Linear xs) = any (\(a:+b,_) -> wrong' a || wrong' b) xs
          wrong' a = isNaN a || isInfinite a
          checkO v v' = if magnitude (dot v v' :: Cplx) > 0.0001 then error "Not orthogonal." else v'
