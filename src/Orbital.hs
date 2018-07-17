{-# LANGUAGE TupleSections #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
module Orbital(
    Orbital(..),
    Matrix,
    Label,
    Integrals,
    evalOrbital,
    calculateIntegrals,
    nuclearHamiltonian,
    eeHamiltonian,
    hartreeFockIterants,
    hartreeFockStep,
    totalEnergy,
    matTimes,
    normalizeWith,
    addMat,
    invert,
    doInvert,
    tabulate,
    showMatrix,
    traceShowMatId,
) where

import Util
import Linear
import Gaussian
import Potential
import Atom
import Eigen

import qualified Data.Map as M
import Data.Map (Map)
import Data.Maybe
import Data.List
import Data.Function
import Data.Bifunctor
import Data.Monoid
import Control.Applicative
import GHC.TypeLits
import Debug.Trace
import GHC.Stack

type Label = (AtomLabel, OrbitalLabel)
type Orbital = (Maybe Spin, Linear Rl Label)
type Matrix a = Map a (Linear Rl a)
-- Integrals = (overlaps, nuclear hamiltonian, four-electron integrals)
type Integrals = (Matrix Label, Matrix Label, Map Label (Map Label (Matrix Label)))

evalOrbital :: KnownNat n => Atoms n -> Linear Rl Label -> Gaussians n
evalOrbital as o = reduce $ o >>= (\(al,ol) -> atomOrbitalsGlobal (as M.! al) M.! ol)

nuclearHamiltonian :: KnownNat n => Atoms n -> Matrix Label
nuclearHamiltonian ats = M.unions $ map atomH (M.toList ats)
    where --atomH :: (AtomLabel, Atom n) -> Matrix Label
          atomH (al, at) = M.mapKeysMonotonic (al,) $ M.mapWithKey (\ol o -> orbH ol o al) (atomOrbitalsGlobal at)
          orbH ol o al = reduce $ approximate (overlapsH o)
          overlapsH o = Linear $ flip map allOrbs (swap . second (singleH o))
          singleH o o' = totalPotential (liftA2 multiply o o') - 0.5*flatten (liftA2 (dot . laplacian) o o')
          approximate o = trim $ reduce $ o >>= ((trim <$> doInvert (overlaps allOrbs)) M.!)
          allOrbs = M.toList $ allOrbitals ats
          totalPotential o = sum $ map (flip atomPotentialGlobal o) (M.elems ats)

-- this ! a ! b ! c = int(r) (r-r')^(2-d)|c⟩⟨a|δ(r)|b⟩
fourElectronIntegrals :: UsableDimension n => Atoms n -> Map Label (Map Label (Matrix Label))
fourElectronIntegrals ats = strict $ tabulate allLabels (\a -> tabulate allLabels (\b -> tabulate allLabels (col a b)))
    where col :: Label -> Label -> Label -> Linear Rl Label
          col a b c = traceCol a b c $ approximate $ mapToLinear $ tabulate allLabels (getFei a b c)
          getFei a b c d = let (a',b',c',d') = semisort (a,b,c,d) in cache M.! a' M.! b' M.! c' M.! d'
          semisort = (\(a,b,c,d) -> if a > c then (c,d,a,b) else (a,b,c,d)) . (\(a,b,c,d) -> (min a b, max a b, min c d, max c d))
          -- Just compute the whole thing then lazily copy some bits over others, thus reducing the amount of actual integrals performed.
          cache = flip fmap allOrbs (\a -> flip fmap allOrbs (\b -> flip fmap allOrbs (\c -> flip fmap allOrbs (fei a b c))))
          fei a b c d = coulumbPotential (convolve . reverseGauss <$> (multiply <$> a <*> b) <*> (multiply <$> c <*> d))
          allLabels = M.keys allOrbs
          allOrbs = allOrbitals ats
          approximate = trim . reduce . ((doInvert (overlaps (M.toList allOrbs)) M.!) =<<)
          traceCol a b c x = seq (x == x) $ trace ("col" ++ show a ++ show b ++ show c) x
          strict x = seq (x == x) x

-- Produces [up electron hamiltonian, down electron hamiltonian]
eeHamiltonian :: Map Label (Map Label (Matrix Label)) -> [Orbital] -> [Matrix Label]
eeHamiltonian fei orbs = foldr (zipWith addMat) [M.empty,M.empty] $ map orbField orbs
    where orbField (s,o) = case s of
              Just Up   -> [addMat (coulumb o) (exchange o), coulumb o]
              Just Down -> [coulumb o, addMat (coulumb o) (exchange o)]
              Nothing   -> [addMat ((2::Rl) *~ coulumb o) (exchange o)]
          coulumb o = flatten' ((tei o M.!) <$> o)
          exchange o = ((-1::Rl) *~ (flip matTimes o <$> tei o))
          tei o = fmap (flatten' . (<$> o) . (M.!)) fei
          flatten' (Linear ms) = foldr addMat M.empty $ map (uncurry (*~)) ms --Can't just use flatten as Map has the wrong monoid instance.

hartreeFockIterants :: UsableDimension n => Atoms n -> Int -> [[Orbital]]
hartreeFockIterants ats n = map snd $ iterate (hartreeFockStep 0.5 n (calculateIntegrals ats)) ([M.empty,M.empty],[])

hartreeFockStep :: Rl -> Int -> Integrals -> ([Matrix Label],[Orbital]) -> ([Matrix Label],[Orbital])
hartreeFockStep s n (overlaps,nh,fei) (peeh,orbs) = (eeh, map (fmap (normalizeWith overlaps) . snd) $ take n $ foldr1 merge newOrbs)
    where eeh = zipWith (\a b -> addMat ((1-s) *~ a) $ (s *~ b)) peeh (eeHamiltonian fei orbs)
          newOrbs = zipWith (\h s -> map (fmap (Just s,)) $ negativeEigenvecs $ addMat nh h) eeh [Up,Down]

-- Doesn't work with orbitals that don't have a spin.
totalEnergy :: forall n. KnownNat n => Atoms n -> Integrals -> [Orbital] -> Rl
totalEnergy ats (overlaps, nh, fei) orbs = sum (map orbEnergy orbs) + atomEnergy
    where orbEnergy (s,o) = matTimes overlaps o `dot` (matTimes nh o <> (0.5::Rl) *~ matTimes (getEeh s) o) / matTimes overlaps o `dot` o
          getEeh (Just s) = eeh !! fromEnum s
          eeh = eeHamiltonian fei $ map (fmap (normalizeWith overlaps)) orbs
          atomEnergy = sum $ map (uncurry atomLabelPairEnergy) $ filter (uncurry (>)) $ (,) <$> atList <*> atList
          atList = M.keys ats
          atomLabelPairEnergy x y = atomPairEnergy (ats M.! x) (ats M.! y)
          atomPairEnergy x y = fromIntegral (atomicNumber x * atomicNumber y) * dist x y ^^ (2 - (natVal @n) Proxy)
          dist x y = sqrt (norm2 (zipWith (-) (atomPos x) (atomPos y)))

overlaps :: (InnerProduct Rl v, Ord a) => [(a,v)] -> Matrix a
overlaps xs = M.fromList $ flip map xs (second $ \x -> Linear (map (\(l,x') -> (dot x x',l)) xs))

orbitalOverlaps :: KnownNat n => Atoms n -> Matrix Label
orbitalOverlaps = overlaps . M.toList . allOrbitals

allOrbitals :: KnownNat n => Atoms n -> Map Label (Gaussians n)
allOrbitals = mconcat . map (\(al,at) -> M.mapKeysMonotonic (al,) $ atomOrbitalsGlobal at) . M.toList

calculateIntegrals :: UsableDimension n => Atoms n -> Integrals
calculateIntegrals ats = (orbitalOverlaps ats', nuclearHamiltonian ats', fourElectronIntegrals ats')
    where ats' = trimOrbitals ats

-- Remove superfluous orbitals which are almost the same as (linear combinations of) others
trimOrbitals :: KnownNat n => Atoms n -> Atoms n
trimOrbitals ats = M.mapWithKey trimAt ats
    where trimAt al at = at{atomOrbitals = M.filterWithKey (\ol _ -> M.member (al,ol) trimmedOrbs) (atomOrbitals at)}
          trimmedOrbs = trimOrbs $ orbitalOverlaps ats
          trimOrbs ovs = let
                  x = eigenvecNear ovs 0
                  d = rayleighQuotient ovs x
                  w = snd $ maximum $ first abs <$> (\(Linear a) -> a) x
                  next = M.delete w $ (reduce . ((\l -> if l == w then mempty else return l) =<<)) <$> ovs
              in if abs d > 2e-3 then ovs else trimOrbs next

doInvert :: forall a. (Ord a) => Matrix a -> Matrix a
doInvert m = maybe (error "Singular matrix") id $ invert m

matTimes :: (HasCallStack, Ord a) => Matrix a -> Linear Rl a -> Linear Rl a
matTimes m v = reduce $ (m M.!) =<< v

normalizeWith :: (InnerProduct Rl a, Ord a) => Matrix a -> Linear Rl a -> Linear Rl a
normalizeWith m o = positiveMultiple $ (1/sqrt (dot o (matTimes m o))::Rl) *~ o
    where positiveMultiple o' = if flatten (const (1::Rl) <$> o') < 0 then (-1::Rl) *~ o' else o'

invert :: forall a. (Ord a) => Matrix a -> Maybe (Matrix a)
invert m0 = invert' m0 m0' xs0
    where xs0 = M.keys m0
          m0' = tabulate xs0 return :: Matrix a
          invert' :: Matrix a -> Matrix a -> [a] -> Maybe (Matrix a)
          invert' m m' [] = Just m'
          invert' m m' (x:xs) =
              let v = m M.! x
                  (Linear vl) = v
                  (a, y) = maximumBy (on compare (abs . fst)) $ dropWhile ((<x).snd) vl ++ [(0,error "singular matrix")]
                  f z = if z == x then y else if z == y then x else z
                  k :: a -> Linear Rl a
                  k z = if x /= z then return z else (1+1/a) *~ return x <> (-(1/a)) *~ (f <$> v)
                  k' = reduce . (>>= k) . fmap f
              in if a == 0 then Nothing else invert' (k' <$> m) (k' <$> m') xs

trim :: (Fractional f, Ord f) => Linear f a -> Linear f a
trim (Linear xs) = Linear $ xs --filter ((>threshold) . abs . fst) xs
    where threshold = (0.001*) $ maximum $ map (abs . fst) xs

tabulate :: (Ord k) => [k] -> (k -> a) -> Map k a
tabulate ks f = M.fromList $ map (\k -> (k, f k)) ks

mapTranspose :: (Ord k, Ord k') => Map k (Map k' a) -> Map k' (Map k a)
mapTranspose = M.foldl (M.unionWith M.union) M.empty . M.mapWithKey (\k m' -> M.singleton k <$> m')

addMat :: Ord a => Matrix a -> Matrix a -> Matrix a
addMat = M.unionWith (\a b -> reduce (a <> b))

showMatrix :: (Show a, Ord a) => Matrix a -> String
showMatrix m = intercalate "\n" $ zipWith (++) xsl $ map (intercalate ", ") $ transpose $ map registerColumn $ zipWith ((:) . show) xs $ map (map show) cs
    where xs = M.keys m
          cs = map (linearToList xs . reduce) $ M.elems m
          linearToList [] (Linear []) = []
          linearToList xs (Linear []) = map (const 0) xs
          linearToList (x':xs) l@(Linear ((a,x):ls))
              | x == x'   = a:linearToList xs (Linear ls)
              | otherwise = 0:linearToList xs l
          registerColumn c = let n = maximum $ map length c in map (take n . (++ repeat ' ')) c
          xsl = registerColumn $ "" : map ((++": ") . show) xs

traceShowMatId :: (Show a, Ord a) => Matrix a -> Matrix a
traceShowMatId m = trace (showMatrix m) m

debugPotential :: UsableDimension n => Atoms n -> Matrix Label
debugPotential ats = (Linear . map swap . M.toList . flip fmap orbs . (p .) . liftA2 multiply) <$> orbs
    where orbs = allOrbitals ats
          p gs = sum $ flip atomPotentialGlobal gs <$> M.elems ats
