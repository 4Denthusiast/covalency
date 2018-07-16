{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE StandaloneDeriving #-}
module Polynomial(
    Polynomial,
    constant,
    variable,
    degree,
    monomial,
    monomialSum,
    monomialSum',
    evaluate,
    evaluate',
    laplacian,
    sphericalHarmonicPolys,
    NamedDimensions(..),
) where

import Linear

import Data.List
import Data.Monoid hiding ((<>))
import Data.Semigroup
import GHC.TypeLits --(natVal, Nat, KnownNat)

data Polynomial (n :: Nat) a = Polynomial [Monomial n a]
data Monomial   (n :: Nat) a = Monomial{
    monCoefficient :: a,
    monExponents   :: [Int]
}
instance Eq a => Eq (Polynomial n a)
instance Ord a => Ord (Polynomial n a)

data Proxy (n::Nat) = Proxy
constant :: forall a n. (KnownNat n, Eq a, Num a) => a -> Polynomial n a

variable :: forall a n. (KnownNat n, Eq a, Num a) => Int -> Polynomial n a

degree :: Polynomial n a -> Int

monomial :: forall n a. (Num a, KnownNat n) => [Int] -> Polynomial n a

instance (Eq a, Num a, KnownNat n) => Num (Polynomial n a) where

instance (Eq a, Num a, KnownNat n) => Semigroup (Polynomial n a)
instance (Eq a, Num a, KnownNat n) => Monoid    (Polynomial n a)
instance (Eq a, Num a, KnownNat n) => Vector a  (Polynomial n a)
instance {-# OVERLAPS #-} (Vector a b, Num b, Eq b, KnownNat n) => Vector a (Polynomial n b)

monomialSum  :: (Monoid a, Vector b a) => ([Int] -> a) -> Polynomial n b -> a
monomialSum' :: (Monoid b, Vector a b) => ([Int] -> a) -> Polynomial n b -> b

evaluate  :: (Monoid a, Num a, Vector b a) => [a] -> Polynomial n b -> a
evaluate' :: (Monoid b, Num a, Vector a b) => [a] -> Polynomial n b -> b

laplacian :: (KnownNat n, Eq a, Num a) => Polynomial n a -> Polynomial n a

sphericalHarmonicPolys :: (KnownNat n) => Int -> [Polynomial n Rl]

class NamedDimensions (n :: Nat) where
    dimName :: forall proxy. proxy n -> Int -> String
instance NamedDimensions 4
instance NamedDimensions 3
instance NamedDimensions 1
instance (Show a, Num a, Eq a, NamedDimensions n, KnownNat n) => Show (Polynomial n a)
