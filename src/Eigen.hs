{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module Eigen(
    inverseIterants,
    rayleighIterate,
    negativeEigenvecsFrom,
    negativeEigenvecs,
    removeKernel,
    eigenvecNear,
    rayleighQuotient,
) where

import Linear
import Atom
import {-#SOURCE#-} Orbital

import qualified Data.Map as M
import qualified Data.Set as S
import Data.List
import Data.Function
import Data.Maybe
import Data.Monoid
import Debug.Trace

-- This doesn't make much physical sense, but it's handy for an arbitrary normalization factor.
instance InnerProduct Rl (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0

idMat :: Ord a => [a] -> Matrix a
idMat = M.fromList . map (\x -> (x,return x))

rayleighQuotient :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Rl
rayleighQuotient m v = dot v (matTimes m v) / dot v v

-- A vector in the v. space on which the matrix acts. This is somewhat less likely to be a member of a proper invariant subspace than simply taking one element, though it still isn't guaranteed.
arbitrary :: Ord a => Matrix a -> Linear Rl a
arbitrary = reduce . foldr (<>) mempty

-- Can't exclude linear combinations of eigenvectors with almost-equal values.
eigenvectorQuality :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Rl
eigenvectorQuality m v0 = norm (v' <> (-dot v v' ::Rl) *~ v)
    where v = normalize @Rl v0
          v' = matTimes m v
          norm x = dot x x

powerConvergent :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a
powerConvergent m = fromJust $ find ((<0.001) . eigenvectorQuality m) $ iterate (normalize @Rl . matTimes m) $ snd $ M.findMin m

offsetMat :: (Ord a) => Matrix a -> Rl -> Matrix a
offsetMat m μ = M.mapWithKey ((reduce .) . (<>) . ((-μ) *~) . return) m

-- (m - μ I)^-1
offsetInverse :: (Ord a) => Matrix a -> Rl -> Maybe (Matrix a)
offsetInverse m μ = invert $ offsetMat m μ

inverseIterants :: (InnerProduct Rl a, Ord a, Show a) => Matrix a -> Rl -> Linear Rl a -> [Linear Rl a]
inverseIterants m μ v = maybe (repeat $ head ker) actuallyIterate $ offsetInverse m μ
    where actuallyIterate m' = iterate (normalize @Rl . matTimes m') v
          (ker,_) = removeKernel $ offsetMat m μ

rayleighIterate :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Linear Rl a
rayleighIterate m v = case offsetInverse m $ rayleighQuotient m v of
    Nothing -> v
    Just m' -> let v' = normalize @Rl $ matTimes m' v in
        if eigenvectorQuality m v' > eigenvectorQuality m v then v else rayleighIterate m v'

slowRayleighIterate :: (InnerProduct Rl a, Ord a) => Matrix a -> Rl -> Linear Rl a -> Linear Rl a
slowRayleighIterate m μ v = case offsetInverse m μ of
    Nothing -> v
    Just m' -> let v' = normalize @Rl $ matTimes m' v in
        if min 1e-12 (eigenvectorQuality m v') >= eigenvectorQuality m v then v else slowRayleighIterate m (0.7*μ+0.3*rayleighQuotient m v) v'

eigenvecNear :: (InnerProduct Rl a, Ord a, Show a) => Matrix a -> Rl -> Linear Rl a
eigenvecNear m μ0 = slowRayleighIterate m μ0 $ fromJust $ find ((<0.1).eigenvectorQuality m) $ drop 20 $ inverseIterants m μ0 $ arbitrary m

negativeEigenvecs :: (InnerProduct Rl a, Ord a, Show a) => Matrix a -> [(Rl,Linear Rl a)]
negativeEigenvecs m = negativeEigenvecsFrom m b
    where b = converged id $ iterate (eigenvalNear m . (\a -> minimum [a-1,a*18,-a])) 0

-- Sometimes misses eigenvectors due to (presumably) numerical instability.
negativeEigenvecsFrom :: (InnerProduct Rl a, Ord a, Show a) => Matrix a -> Rl -> [(Rl,Linear Rl a)]
negativeEigenvecsFrom m b = map traceQuality $ concatMap (\μ -> map (μ,) $ fst $ removeKernel $ offsetMat m $ eigenvalNear m μ) $ takeWhile (<0) $ eigenvalsFrom m b
    where traceQuality (μ,v) = trace ("quality: " ++ show (eigenvectorQuality m v)) (μ,v)

eigenvalsFrom :: (InnerProduct Rl a, Ord a, Show a) => Matrix a -> Rl -> [Rl]
eigenvalsFrom m b = if M.null m then [] else seq m' b' : eigenvalsFrom m' b'
    where b' = rayleighQuotient m $ eigenvecNear m b
          m' = removeEigenval b' m

converged f (x:x':xs) = if f x' >= f x then x else converged f (x':xs)
converged f [x] = x

-- (removeEigenval μ m) is a matrix similar to the restriction of m to a subspace complementary to the μ-eigenspace.
removeEigenval :: (InnerProduct Rl a, Ord a, Show a) => Rl -> Matrix a -> Matrix a
removeEigenval μ m = flip offsetMat (-μ) $ uncurry (traceShow . (μ,) . checkSomething . length) $ removeKernel $ offsetMat m μ
    where checkSomething 0 = error ("Missing eigenvec: "++show μ++" in Iμ + \n"++showMatrix (offsetMat m μ))
          checkSomething n = n

removeKernel :: forall a. (InnerProduct Rl a, Ord a, Show a) => Matrix a -> ([Linear Rl a], Matrix a)
removeKernel m0 = removeKernel' m0 (idMat xs0) xs0 []
    where xs0 = M.keys m0
          checkKer k = k--trace ("zero? "++show (norm (matTimes m0 k) / norm k)) k
          norm v = dot @Rl v v
          without ker z = if elem z ker then mempty else return z
          -- invariant: m' * m = m0
          removeKernel' m m' (x:xs) ker = {-trace ("m:\n"++showMatrix m++"\nm':\n"++showMatrix m'++"\nm' * m:\n"++showMatrix (matTimes m' <$> m)++"\n")-} $ let
                  v = m M.! x
                  (Linear vl) = v
                  (a,y) = maximumBy (on compare (abs . fst)) $ dropWhile ((<x).snd) vl ++ [(0,undefined)]
                  (a',_) = maximumBy (on compare (abs . fst)) $ takeWhile ((<x).snd) vl ++ [(0,undefined)]
                  f :: a -> a
                  f z = if z == x then y else if z == y then x else z
                  k :: a -> Linear Rl a
                  k  z = if z == x then (-1/a) *~ (f <$> v) <> (1+1/a) *~ return x else return z
                  k' :: Matrix a
                  k' = flip M.union (idMat xs0) $ M.singleton x (f <$> v)
                  succeed = removeKernel' ((reduce . ((k.f) =<<)) <$> m) (matTimes m' <$> fmap (f <$>) k') xs ker
                  fail = removeKernel' m m' xs (x:ker)
              in
                  if abs a > 1e-6* abs a' + 1e-14 then succeed else fail
          removeKernel' m m' [] ker = (
                  map checkKer $ map (\x -> normalize @Rl $ reduce $ (without ker =<< m M.! x) <> (-1::Rl) *~ return x) ker,
                  (without ker =<<) <$> matTimes m <$> M.withoutKeys m' (S.fromList ker)
              )

--Testing
instance InnerProduct Rl Int where
    dot a b = if a == b then 1 else 0

matrixFromList :: [[Rl]] -> Matrix Int
matrixFromList = M.fromList . zip [0..] . map (Linear . flip zip [0..])

eigenvalNear :: (Ord a, InnerProduct Rl a, Show a) => Matrix a -> Rl -> Rl
eigenvalNear m b = rayleighQuotient m $ eigenvecNear m b
