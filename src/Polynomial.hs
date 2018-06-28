{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Polynomial(
    Polynomial,
    constant,
    variable,
    degree,
    monomialSum,
    monomialSum',
    evaluate,
    evaluate',
    laplacian,
    sphericalHarmonicPolys,
) where

import Linear
import Orbital
import Eigen

import Data.List
import Data.Complex
import Data.Monoid hiding ((<>))
import Data.Semigroup
import GHC.TypeLits --(natVal, Nat, KnownNat)

-- Polynomials in n variables.
data Polynomial (n :: Nat) a = Polynomial [Monomial n a]
deriving instance Eq a => Eq (Polynomial n a)
deriving instance Ord a => Ord (Polynomial n a)
data Monomial   (n :: Nat) a = Monomial{
    monCoefficient :: a,
    monExponents   :: [Int]
} deriving (Eq, Ord)

monomialTimes :: Num a => Monomial n a -> Monomial n a -> Monomial n a
monomialTimes (Monomial x ex) (Monomial y ey) = Monomial (x*y) (zipWith (+) ex ey)

instance Functor (Monomial n) where
    fmap f (Monomial a e) = Monomial (f a) e

simplifyPoly :: (Eq a, Num a) => [Monomial n a] -> Polynomial n a
simplifyPoly = Polynomial . filter ((/=0) . monCoefficient) . mergeBy mergeMon . sortBy compareExponents
    where compareExponents (Monomial _ e0) (Monomial _ e1) = mappend (compare (sum e1) (sum e0)) (compare e1 e0)
          mergeMon (Monomial x ex) (Monomial y ey)
              | ex == ey  = Just $ Monomial (x+y) ex
              | otherwise = Nothing

data Proxy (n::Nat) = Proxy
constant :: forall a n. (KnownNat n, Eq a, Num a) => a -> Polynomial n a
constant x = simplifyPoly [Monomial x $ genericReplicate (natVal @n Proxy) 0]

variable :: forall a n. (KnownNat n, Eq a, Num a) => Int -> Polynomial n a
variable i = Polynomial [Monomial 1 es]
    where n  = fromInteger $ natVal @n Proxy
          es = replicate i 0 ++ 1 : replicate (n-i-1) 0

degree :: Polynomial n a -> Int
degree (Polynomial []) = -1 -- negative infinity is the usual convention, but making this an Int seems more appropriate.
degree (Polynomial (Monomial _ es : _)) = sum es

instance (Eq a, Num a, KnownNat n) => Num (Polynomial n a) where
    (Polynomial xs) + (Polynomial ys) = simplifyPoly (xs ++ ys)
    negate (Polynomial xs)            = Polynomial (map (fmap negate) xs)
    (Polynomial xs) * (Polynomial ys) = simplifyPoly (monomialTimes <$> xs <*> ys)
    abs _ = error "The typeclass Num doesn't contain a method: conjugate"
    signum _ = error "Signum doesn't make much sense for polynomials."
    fromInteger = constant . fromInteger

instance (Eq a, Num a, KnownNat n) => Semigroup (Polynomial n a) where (<>) = (+)
instance (Eq a, Num a, KnownNat n) => Monoid    (Polynomial n a) where {mappend = (<>); mempty = 0}
instance (Eq a, Num a, KnownNat n) => Vector a  (Polynomial n a) where (*~) = (*) . constant
instance {-# OVERLAPS #-} (Vector a b, Num b, Eq b, KnownNat n) => Vector a (Polynomial n b) where
    (*~) = towerScale @b

monomialSum  :: (Monoid a, Vector b a) => ([Int] -> a) -> Polynomial n b -> a
monomialSum  f (Polynomial ms) = mconcat $ map (\(Monomial b es) -> b *~ f es) ms
monomialSum' :: (Monoid b, Vector a b) => ([Int] -> a) -> Polynomial n b -> b
monomialSum' f (Polynomial ms) = mconcat $ map (\(Monomial b es) -> f es *~ b) ms

evaluate  :: (Monoid a, Num a, Vector b a) => [a] -> Polynomial n b -> a
evaluate  xs = monomialSum  (product . zipWith (^) xs)
evaluate' :: (Monoid b, Num a, Vector a b) => [a] -> Polynomial n b -> b
evaluate' xs = monomialSum' (product . zipWith (^) xs)

laplacian :: (KnownNat n, Eq a, Num a) => Polynomial n a -> Polynomial n a
laplacian (Polynomial ms) = sum $ map monLap ms
    where monLap (Monomial a es) = Polynomial $ zipWith3 (\h t e -> Monomial (a*fromIntegral (e*(e-1))) (h++(e-2):t)) (inits es) (tail $ tails es) es

sphericalHarmonicPolys :: forall n. KnownNat n => Int -> [Polynomial n Cplx]
sphericalHarmonicPolys = map toPoly . fst . removeKernel . flip tabulate linearLaplacian . allMons (fromIntegral $ natVal @n Proxy)
    where toPoly = flatten . fmap (Polynomial . (:[]) . Monomial 1)
          allMons n l = filter ((==l) . sum) $ sequence $ replicate n [0..l]
          linearLaplacian :: [Int] -> Linear Cplx [Int]
          linearLaplacian m = mconcat $ zipWith (\h (e:t) -> if e < 2 then mempty else (e*(e-1))*~return (shift$h++(e-2):t)) (inits m) (init $ tails m)
          shift (x:xs) = x+2 : xs

instance InnerProduct Cplx [Int] where dot a b = if a == b then 1 else 0


assumePolyReal :: (Eq a, Num a) => Polynomial n (Complex a) -> Polynomial n a
assumePolyReal (Polynomial ms) = Polynomial $ map (fmap assumeReal) ms
    where assumeReal (x :+ 0) = x
          assumeReal _        = error "Real number required."

instance Semilinear a => Semilinear (Polynomial n a) where
    conj (Polynomial ms) = Polynomial $ map (fmap conj) ms


class NamedDimensions (n :: Nat) where
    dimName :: forall proxy. proxy n -> Int -> String

instance NamedDimensions 4 where
    dimName _ 0 = "x"
    dimName _ 1 = "y"
    dimName _ 2 = "z"
    dimName _ 3 = "w"

instance NamedDimensions 3 where
    dimName _ 0 = "x"
    dimName _ 1 = "y"
    dimName _ 2 = "z"

instance NamedDimensions 1 where
    dimName _ 0 = "r" --"r" is more likely to get used, but "x" is more general.

instance (Show a, Num a, Eq a, NamedDimensions n, KnownNat n) => Show (Polynomial n a) where
    show 0 = "0"
    show (Polynomial ms) = intercalate " + " $ map show ms

instance (Show a, Num a, Eq a, NamedDimensions n, KnownNat n) => Show (Monomial n a) where
    show (Monomial a es) = showCoeff ++ concatMap (\(i, e) -> dimName @n Proxy i ++ showPower e) (filter ((>0) . snd) $ zip [0..] es)
        where showCoeff   = if a == 1 && any (/=0) es then "" else show a
              showPower e = if e == 1 then "" else '^' : show e

mergeBy :: (a -> a -> Maybe a) -> [a] -> [a]
mergeBy f = mergeBy' []
    where mergeBy' ys [] = ys
          mergeBy' ys (x:xs) = mergeBy' (mergeInto ys x) xs
          mergeInto [] x = [x]
          mergeInto (y:ys) x = case f y x of
              Just y' -> y':ys
              Nothing -> y: mergeInto ys x
