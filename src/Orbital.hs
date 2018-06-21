{-# LANGUAGE TupleSections #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Orbital(
    Orbital(..),
    evalOrbital,
    nuclearHamiltonian,
) where

import Linear
import Gaussian
import Potential
import Atom

import qualified Data.Map as M
import Data.Maybe
import Data.Complex
import Data.Bifunctor
import Data.Monoid
import Control.Applicative
import GHC.TypeLits

type Orbital = Linear Cplx (AtomLabel, OrbitalLabel)
type Matrix a = M.Map a (Linear Cplx a)

evalOrbital :: Atoms n -> Orbital -> [Float] -> Cplx
evalOrbital as o xs = flatten $ (\(al, ol) -> evalAtomOrb (as M.! al) ol xs) <$> o

nuclearHamiltonian :: KnownNat n => Atoms n -> Matrix (AtomLabel, OrbitalLabel)
nuclearHamiltonian ats = M.unions $ map atomH (M.toList ats)
    where --atomH :: (AtomLabel, Atom n) -> Matrix (AtomLabel, OrbitalLabel)
          atomH (al, at) = M.mapKeysMonotonic (al,) $ M.mapWithKey (\ol o -> orbH ol o (atomKineticTerm at M.! ol) al) (atomOrbitalsGlobal at)
          --orbH :: OrbitalLabel -> Gaussians n -> Linear Cplx OrbitalLabel -> AtomLabel -> Linear Cplx (AtomLabel, OrbitalLabel)
          orbH ol o k al = fmap (al,) k <> approximate (overlapsH o)
          overlapsH o = Linear $ flip map allOrbs (swap . second (totalPotential . liftA2 multiply o))
          approximate :: Linear Cplx (AtomLabel, OrbitalLabel) -> Orbital
          approximate o = trim $ reduce $ o >>= ((trim <$> invert (overlaps allOrbs)) M.!)
          allOrbs = concatMap (\(al,at) -> map (first (al,)) $ M.toList $ atomOrbitalsGlobal at) (M.toList ats)
          totalPotential o = sum $ map (flip atomPotentialGlobal o) (M.elems ats)

overlaps :: (InnerProduct v Cplx, Ord a) => [(a,v)] -> Matrix a
overlaps xs = M.fromList $ flip map xs (second $ \x -> Linear (map (\(l,x') -> (dot x x',l)) xs))

invert :: (Ord a) => Matrix a -> Matrix a
invert m0 = invert' m0 m0' xs0
    where xs0 = M.keys m0
          m0' = tabulate xs0 return
          invert' m m' [] = m'
          invert' m m' (x:xs) =
              let k y = if x /= y then return y else invertCol (m M.! x) x
                  k' = reduce . (>>= k)
              in invert' (k' <$> m) (k' <$> m') xs
          invertCol (Linear xs) x = let ax = fromJust $ lookup x $ map swap xs in Linear ((1/ax,x):map (first (/(-ax))) (filter ((/=x).snd) xs))

trim :: (Fractional f, Ord f) => Linear f a -> Linear f a
trim (Linear xs) = Linear $ filter ((>threshold) . abs . fst) xs
    where threshold = (0.001*) $ maximum $ map (abs . fst) xs

tabulate :: (Ord k) => [k] -> (k -> a) -> M.Map k a
tabulate ks f = M.fromList $ map (\k -> (k, f k)) ks

swap (x,y) = (y,x)
