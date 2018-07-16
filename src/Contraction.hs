{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
module Contraction(
    loadAtomWithBasis,
    basisify,
) where

import Util
import Linear
import qualified Polynomial as P
import Gaussian
import Atom
import Orbital
import Eigen

import Data.List
import Data.Maybe
import qualified Data.Map as M
import Control.Applicative
import Control.Monad
import Data.Monoid
import Data.Bifunctor
import System.Posix.Files
import GHC.TypeLits
import Debug.Trace

type BasisSet = [(L,Linear Rl Rl)]

minimalBasisSet :: forall n. UsableDimension n => Int -> BasisSet
minimalBasisSet = minimalBasisSetFrom @n (-1) 1 1

minimalBasisSetFrom :: forall n. UsableDimension n => Int -> Int -> L -> Int -> BasisSet
minimalBasisSetFrom inner outer l0 z = if insufficientBounds then expanded else b'
    where insufficientBounds = inner' < inner || outer' > outer || l0' > l0
          expanded = minimalBasisSetFrom @n inner' outer' l0' z 
          b' = fmap (fmap (fmap expSize)) b
          inner' = minimum ns - 1
          outer' = maximum ns + 1
          l0' = maximum (map fst b) + 1
          ns = foldr union [] $ map (linEls . snd) b
          linEls (Linear xs) = map snd xs
          b = minimalBasisOf (testAtom @n inner outer l0 z)

testAtom :: forall n. UsableDimension n => Int -> Int -> L -> Int -> Atom n
testAtom inner outer l0 z = changeZ z $ Atom {
        atomPos = genericReplicate (natVal @n Proxy) 0,
        atomOrbitals = (M.fromList $ liftA2 (\i (l,m,p) -> ((i,l,m),normalize @Rl $ return $ centralGaussian (expSize i) p)) [inner..outer] (concatMap (\l -> zipWith (l,,) [0..] $ P.sphericalHarmonicPolys l) [0..l0]))
    }

expSize :: Floating a => Int -> a
expSize = (exp . (1.5*) . fromIntegral)

minimalBasisOf :: forall n. UsableDimension n => Atom n -> [(L,Linear Rl Int)]
minimalBasisOf at = fmap (fmap trimBasisEl) $ basisConverged $ map snd $ iterate (stepBasis @n z ints) (M.empty, [])
    where ints = calculateIntegrals $ M.singleton "" at
          trimBasisEl (Linear xs) = Linear $ filter ((>1e-2) . abs . fst) xs
          z = atomicNumber at

basisConverged :: [[(L,Linear Rl Int)]] -> [(L,Linear Rl Int)]
basisConverged (x:x':xs) = if similar x x' then x' else basisConverged (x':xs)
    where similar (y:ys) (y':ys') = elSimilar y y' && similar ys ys'
          similar [] [] = True
          similar _ _ = False
          elSimilar (l,b) (l',b') = l == l' && norm (b <> (-1::Rl) *~ b') < 1e-10
          norm :: Linear Rl Int -> Rl
          norm b = dot b b

-- Iterate the basis towards convergence, taking the occupancy to smoothly drop to 0 wrt energy.
stepBasis :: forall n. KnownNat n => Int -> Integrals -> (Matrix Label,[(L,Linear Rl Int)]) -> (Matrix Label,[(L,Linear Rl Int)])
stepBasis z (ov,nh,fei) (peeh,b) = (eeh,b')
    where h = addMat nh peeh
          -- Split a matrix representing a spherically-symmetric operator into components wrt L.
          splitMat :: Matrix Label -> [Matrix Int]
          splitMat =
              fmap (fmap (fmap getN)) .
              map (M.mapKeysMonotonic getN) .
              unfoldr (\(m,l) -> guard (not (M.null m)) >> (Just $ (,l+1) <$> M.partitionWithKey (\(_,(_,l',0)) _ -> l' == l) m)) .
              (,0) . M.filterWithKey (\(_,(_,_,m)) _ -> m == 0)
          getN (_,(n,_,_)) = n
          orbss = zipWith (\h' ov' -> map (fmap $ normalizeWith ov') (negativeEigenvecs h')) (splitMat h) (splitMat ov)
          orbs = foldr merge [] $ zipWith (\l os -> map (fmap (l,)) os) [0..] orbss
          thisMultiplicity :: Num a => L -> a
          thisMultiplicity = multiplicity (natVal @n Proxy)
          weightedSum f e = sum $ map (\(e',(l,_)) -> f e' * 2*thisMultiplicity l) $ takeWhile ((<e) . fst) orbs
          eInit = fmap fst $ listToMaybe $ dropWhile ((<=fromIntegral z).snd) $ scanl (\(e,c) (e',(l,_)) -> (e',c+2*thisMultiplicity l)) (-1/0,0) orbs
          secSearchNear t f x0 = if f x0 > t then secSearchDown t f (2*x0) x0 else secSearchUp t f x0 (x0/2)
          secSearchDown t f xl xh = if f xl < t then secSearch t f xl xh else secSearchDown t f (2*xl) xh
          secSearchUp   t f xl xh = if f xh > t then secSearch t f xl xh else secSearchUp   t f xl (xh/2)
          secSearch t f xl xh = if xh - xl < 1e-9 then xl else let x' = xl + (xh-xl)*(min 0.9 $ max 0.1 $ (t-f xl)/(f xh-f xl)) in if f x' < t then secSearch t f x' xh else secSearch t f xl x'
          e0 = secSearchNear (fromIntegral z) (\e -> weightedSum (weight e) (e*0.6)) (fromJust eInit)
          weight e = if isNothing eInit then const 1 else exp . negate . (^4) . (e/)
          tOrbs = if isNothing eInit then orbs else takeWhile ((<0.6*e0) . fst) orbs
          neeh = foldr addMat M.empty $ map (\(e,(l,o)) -> weight e0 e *~ head (eeHamiltonian fei (map (\m -> (Nothing,("",) <$> (,l,m) <$> o)) [0..thisMultiplicity l - 1]))) tOrbs
          eeh = addMat ((0.6::Rl) *~ peeh) ((0.4::Rl) *~ neeh)
          b' = map snd tOrbs

-- Explicitly memoise because this is used lots.
multiplicity :: (Integral i, Num n) => i -> L -> n
multiplicity d l = fromIntegral $ mults !! (fromIntegral d) !! l

mults = map (\d' -> map (\l' -> derive d' l') [0..]) [0..]
    where derive d' l' = binomial (d'+l'-2) l' + binomial (d'+l'-3) (l'-1)

binomial :: Int -> Int -> Int
binomial a b
    | a < 0 || b < 0  || b > a = 0
    | b > div a 2              = binomial a (a-b)
    | b == 0                   = 1
    | otherwise                = div (binomial a (b-1) * (a-b+1)) b

loadAtomWithBasis :: forall n. UsableDimension n => [Rl] -> Int -> IO (Atom n)
loadAtomWithBasis xs z = atomWithBasis xs z <$> (loadBasis @n) z

basisify :: UsableDimension n => Atom n -> IO (Atom n)
basisify at = loadAtomWithBasis (atomPos at) (atomicNumber at)

atomWithBasis :: UsableDimension n => [Rl] -> Int -> BasisSet -> Atom n
atomWithBasis xs z basis = changeZ z $ Atom{
        atomPos = xs,
        atomOrbitals = fmap (normalize @Rl) $ M.fromList $ concat $ zipWith (\(l,cs) n -> zipWith (\p m -> ((n,l,m),flip centralGaussian p <$> cs)) (P.sphericalHarmonicPolys l) [0..]) basis [0..]
    }

loadBasis :: forall n. UsableDimension n => Int -> IO BasisSet
loadBasis z = putFolder >> putBasis >> readBasis
    where filePath = "bases/"++show d++"D "++show z++".txt"
          d = natVal @n Proxy
          putFolder = return ()
          putBasis = fileExist filePath >>= \e -> if e then return () else writeFile filePath (show $ minimalBasisSet @n z)
          readBasis = read <$> readFile filePath
