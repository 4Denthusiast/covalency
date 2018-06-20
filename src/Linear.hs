{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE StandaloneDeriving #-}
module Linear(
    Linear(..),
    lmap,
    reduce,
    scale,
    join,
    flatten,
    single,
    Semilinear(..),
    InnerProduct(..),
    normalize,
) where

import Data.List
import Data.Monoid hiding ((<>))
import Data.Semigroup
import Data.Complex

data Linear f a where
    Linear :: (Num f, Eq f, Ord a, Eq a) => [(f,a)] -> Linear f a

deriving instance (Ord f, Ord a) => Ord (Linear f a)
deriving instance (Eq  f, Eq  a) => Eq  (Linear f a)
deriving instance (Show f, Show a) => Show (Linear f a)

lmap :: (Ord b, Eq b) => (a -> b) -> Linear f a -> Linear f b
lmap f (Linear xs) = reduce $ Linear $ map (\(n,x) -> (n,f x)) xs

instance Semigroup (Linear f a) where
    (Linear a) <> (Linear b) = reduce (Linear (a++b))

instance (Num f, Eq f, Ord a, Eq a) => Monoid (Linear f a) where
    mempty = Linear []
    mappend = (<>)

reduce :: Linear f a -> Linear f a
reduce (Linear xs) = Linear $ foldr merge [] $ sortOn snd $ xs
    where merge (0,_) xs = xs
          merge x []     = [x]
          merge (m,y) ((n,x):xs)
              | x /= y    = (m,y):(n,x):xs
              | m+n == 0  = xs
              | otherwise = (m+n,x):xs

scale :: f -> Linear f a -> Linear f a
--scale 0 (Linear xs) = Linear []
scale a (Linear xs) = Linear $ map (\(n,x) -> (a*n, x)) xs

join :: (Ord a, Eq a) => Linear f (Linear f a) -> Linear f a
join (Linear xs) = mconcat $ map (uncurry scale) xs

flatten :: Linear f f -> f
flatten (Linear xs) = sum $ map (uncurry (*)) xs

single :: (Num f, Eq f, Ord a, Eq a) => a -> Linear f a
single x = Linear [(1,x)]

class Semilinear n where
    conj :: n -> n

instance {- OVERLAPS #-} Num n => Semilinear (Complex n) where
    conj = conjugate

instance {-# OVERLAPS #-} (Semilinear n, Semilinear a) => Semilinear (Linear n a) where
    conj (Linear xs) = Linear (map (\(n,x) -> (conj n, conj x)) xs)

instance {-# OVERLAPS #-} Semilinear n where
    conj = id

almap :: (Ord b, Eq b, Semilinear f) => (a -> b) -> Linear f a -> Linear f b
almap f (Linear xs) = reduce $ Linear $ map (\(n,x) -> (conj n, f x)) xs

class InnerProduct a n where
    dot :: a -> a -> n

instance Num n => InnerProduct n n where
    dot = (*) . conj

instance (InnerProduct a n, Num n, Ord n, Ord a, Eq a) => InnerProduct (Linear n a) n where
    dot x y = flatten $ join $ almap (flip lmap y . dot) x

normalize :: InnerProduct (Linear n a) n, Floating n => Linear n a -> Linear n a
normalize l = scale (1/sqrt (dot l l)) l
