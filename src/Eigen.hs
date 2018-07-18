{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module Eigen(
    inverseIterants,
    rayleighIterate,
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
import Data.Bifunctor
import qualified Numeric.LinearAlgebra.Data as ND
import qualified Numeric.LinearAlgebra as N
import Data.Complex
import Debug.Trace

-- This doesn't make much physical sense, but it's handy for an arbitrary normalization factor.
instance InnerProduct Rl (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0
instance InnerProduct Cplx (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0
instance Semilinear (AtomLabel, OrbitalLabel) where

idMat :: Ord a => [a] -> Matrix a
idMat = M.fromList . map (\x -> (x,return x))

rayleighQuotient :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Rl
rayleighQuotient m v = dot v (matTimes m v) / dot v v

-- Can't exclude linear combinations of eigenvectors with almost-equal values.
eigenvectorQuality :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Rl
eigenvectorQuality m v0 = norm (v' <> (-dot v v' ::Rl) *~ v)
    where v = normalize @Rl v0
          v' = matTimes m v
          norm x = dot x x

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
eigenvecNear m μ0 = slowRayleighIterate m μ0 $ fromJust $ find ((<0.1).eigenvectorQuality m) $ snd $ minimumBy (on compare fst) $ (\(v:vs) -> (abs (μ0-rayleighQuotient m v),v:vs)) <$> drop 20 <$> inverseIterants m μ0 <$> return <$> M.keys m

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
          removeKernel' m m' (x:xs) ker = {-trace ("m:\n"++showMatrix m++"\nm':\n"++showMatrix m'++"\nm' * m:\n"++showMatrix (matTimes m' <$> m)++"\n") $-} let
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

matToNumericsMat :: (Ord a, Show a) => Matrix a -> ND.Matrix Rl
matToNumericsMat m = ND.tr $ ND.fromLists $ map (colToList . reduce) $ M.elems m
    where colToList (Linear xs) = match xs (M.keys m)
          match ((a,x):xs) (k:ks) = if x <= k then a : match xs ks else 0 : match ((a,x):xs) ks
          match [] ks = map (const 0) ks

negativeEigenvecs :: forall a. (InnerProduct Cplx a, InnerProduct Rl a, Ord a, Semilinear a, Show a) => Matrix a -> Matrix a -> [(Rl,Linear Rl a)]
negativeEigenvecs ov m = concat $ map orthonormalize $ groupBy (on (==) fst) $ sortOn fst $ filter ((<0) . fst) ps
    where (ne,nv) = fmap ND.tr $ N.eig $ matToNumericsMat m
          ps = zip (map toReal $ ND.toList ne) (toLin <$> ND.toLists nv)
          toReal (a:+b) = if abs b > 1e-10 then error ("Complex energy: "++show (a:+b)) else a
          toLin :: (Num f, Eq f) => [f] -> Linear f a
          toLin = reduce . Linear . flip zip (M.keys m)
          orthonormalize ps' = snd $ mapAccumL (
                  \vs (e,v) -> let v' = reduce $ onmlizeOnce vs v in ((complexify $ matTimes ov v', complexify v'):vs,(e,v'))
              ) [] ps'
          onmlizeOnce [] = normalizeWith ov . (\(Linear xs) -> Linear (first realPart <$> xs)) . rescale
          onmlizeOnce ((f,f'):fs) = onmlizeOnce fs . (\v -> reduce $ v <> (-dot f v::Cplx) *~ f') . rescale
          rescale (Linear xs) = (1/maximumBy (on compare magnitude) (map fst xs)) *~ Linear xs
          complexify = flatten . fmap return

--Testing
instance InnerProduct Rl Int where
    dot a b = if a == b then 1 else 0
instance InnerProduct Cplx Int where
    dot a b = if a == b then 1 else 0

matrixFromList :: [[Rl]] -> Matrix Int
matrixFromList = M.fromList . zip [0..] . map (Linear . flip zip [0..])
