{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Eigen(
    inverseIterants,
    rayleighIterate,
    negativeEigenvecsFrom,
) where

import Linear
import Atom
import Orbital

import qualified Data.Map as M
import qualified Data.Set as S
import Data.List
import Data.Function
import Data.Maybe
import Data.Monoid
import Data.Complex
import Debug.Trace

-- This doesn't make much physical sense, but it's handy for an arbitrary normalization factor.
instance InnerProduct Cplx (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0

matVecTimes :: Ord a => Matrix a -> Linear Cplx a -> Linear Cplx a
matVecTimes m v = reduce $ v >>= (m M.!)

idMat :: Ord a => [a] -> Matrix a
idMat = M.fromList . map (\x -> (x,return x))

rayleighQuotient :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Cplx
rayleighQuotient m v = dot v (matVecTimes m v) / dot v v

-- Can't exclude linear combinations of eigenvectors with almost-equal values.
eigenvectorQuality :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Cplx
eigenvectorQuality m v = norm (v' <> (-(dot v v' / dot v v)::Cplx) *~ v)
    where v' = matVecTimes m v
          norm x = dot x x

powerConvergent :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a
powerConvergent m = fromJust $ find ((<0.001) . eigenvectorQuality m) $ iterate (normalize @Cplx . matVecTimes m) $ snd $ M.findMin m

offsetMat :: (Ord a) => Matrix a -> Cplx -> Matrix a
offsetMat m μ = M.mapWithKey ((reduce .) . (<>) . ((-μ) *~) . return) m

-- (m - μ I)^-1
offsetInverse :: (Ord a) => Matrix a -> Cplx -> Maybe (Matrix a)
offsetInverse m μ = invert $ offsetMat m μ

inverseIterants :: (InnerProduct Cplx a, Ord a) => Matrix a -> Cplx -> Linear Cplx a -> [Linear Cplx a]
inverseIterants m μ v = iterate (normalize @Cplx . matVecTimes (fromJust $ offsetInverse m μ)) v

rayleighIterate :: (InnerProduct Cplx a, Ord a) => Matrix a -> Linear Cplx a -> Linear Cplx a
rayleighIterate m v = case offsetInverse m $ rayleighQuotient m v of
    Nothing -> v
    Just m' -> let v' = normalize @Cplx $ matVecTimes m' v in
        if eigenvectorQuality m v' > eigenvectorQuality m v then v else rayleighIterate m v'

-- Sometimes misses eigenvectors due to (presumably) numerical instability.
negativeEigenvecsFrom :: (InnerProduct Cplx a, Ord a) => Matrix a -> Cplx -> [Linear Cplx a]
negativeEigenvecsFrom m b = concatMap (fst . removeKernel . offsetMat m) $ takeWhile (<0) $ eigenvalsFrom m b

eigenvalsFrom :: (InnerProduct Cplx a, Ord a) => Matrix a -> Cplx -> [Cplx]
eigenvalsFrom m b = if M.null m then [] else b' : eigenvalsFrom (removeEigenval b' m) b'
    where b' = rayleighQuotient m $ rayleighIterate m $ converged (abs . (b-) . eigenvectorQuality m) $ drop 20 $ inverseIterants m b v
          v = snd (M.findMin m)

converged f (x:x':xs) = if f x' >= f x then x else converged f (x':xs)

-- (removeEigenval μ m) is a matrix similar to the restriction of m to a subspace perpendicular to the μ-eigenspace.
removeEigenval :: (InnerProduct Cplx a, Ord a) => Cplx -> Matrix a -> Matrix a
removeEigenval μ m = flip offsetMat (-μ) $ snd $ removeKernel $ offsetMat m μ

removeKernel :: forall a. (InnerProduct Cplx a, Ord a) => Matrix a -> ([Linear Cplx a], Matrix a)
removeKernel m0 = removeKernel' m0 (idMat xs0) xs0 []
    where xs0 = M.keys m0
          -- invariant: m' * m = m0
          removeKernel' m m' (x:xs) ker = {-trace ("m:\n"++showMatrix m++"\nm':\n"++showMatrix m'++"\nm' * m:\n"++showMatrix (matVecTimes m' <$> m)++"\n") $-} let
                  v = m M.! x
                  (Linear vl) = v
                  (a,y) = maximumBy (on compare (magnitude . fst)) $ dropWhile ((<x).snd) vl ++ [(0,undefined)]
                  (a',_) = maximumBy (on compare (magnitude . fst)) $ takeWhile ((<x).snd) vl ++ [(0,undefined)]
                  f :: a -> a
                  f z = if z == x then y else if z == y then x else z
                  k :: a -> Linear Cplx a
                  k  z = if z == x then (-1/a) *~ (f <$> v) <> (1+1/a) *~ return x else return z
                  k' :: Matrix a
                  k' = flip M.union (idMat xs0) $ M.singleton x (f <$> v)
                  succeed = removeKernel' ((reduce . ((k.f) =<<)) <$> m) (matVecTimes m' <$> fmap (f <$>) k') xs ker
                  fail = removeKernel' m m' xs (x:ker)
              in
                  if magnitude a > 0.001* magnitude a' then succeed else fail
          removeKernel' m m' [] ker = (
                  map (\x -> normalize @Cplx $ m M.! x <> (-1::Cplx) *~ return x) ker,
                  ((\z -> if elem z ker then mempty else return z) =<<) <$> matVecTimes m <$> M.withoutKeys m' (S.fromList ker)
              )

--Testing
instance InnerProduct Cplx Int where
    dot a b = if a == b then 1 else 0

matrixFromList :: [[Cplx]] -> Matrix Int
matrixFromList = M.fromList . zip [0..] . map (Linear . flip zip [0..])

eigenvalNear :: (Ord a, InnerProduct Cplx a) => Matrix a -> Cplx -> Cplx
eigenvalNear m b = rayleighQuotient m $ rayleighIterate m $ (!!40) $ inverseIterants m b $ return $ head $ M.keys m
