{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE StandaloneDeriving #-}
module Linear(
    Linear(..),
    reduce,
    scale,
    flatten,
    Semilinear(..),
    InnerProduct(..),
    normalize,
    Cplx,
) where

import Data.List
import Data.Monoid hiding ((<>))
import Data.Semigroup
import Control.Applicative
import Control.Monad
import Data.Bifunctor
import Data.Complex

data Linear f a = Linear [(f,a)] deriving (Ord, Eq, Show)

instance Semigroup (Linear f a) where
    (Linear a) <> (Linear b) = Linear (a++b)

instance Monoid (Linear f a) where
    mempty = Linear []
    mappend = (<>)

instance Functor (Linear f) where
    fmap f (Linear xs) = Linear $ map (second f) xs

instance Num f => Applicative (Linear f) where
    pure x = Linear [(1,x)]
    liftA2 f (Linear xs) (Linear ys) = Linear [(a*b, f x y) | (a,x) <- xs, (b,y) <- ys]

instance Num f => Monad (Linear f) where
    (Linear xs) >>= f = mconcat $ map (\(a,x) -> scale a (f x)) xs

reduce :: (Num f, Eq f, Ord a, Eq a) => Linear f a -> Linear f a
reduce (Linear xs) = Linear $ foldr merge [] $ sortOn snd $ xs
    where merge (0,_) xs = xs
          merge x []     = [x]
          merge (m,y) ((n,x):xs)
              | x /= y    = (m,y):(n,x):xs
              | m+n == 0  = xs
              | otherwise = (m+n,x):xs

scale :: Num f => f -> Linear f a -> Linear f a
scale a (Linear xs) = Linear $ map (\(n,x) -> (a*n, x)) xs

flatten :: Num f => Linear f f -> f
flatten (Linear xs) = sum $ map (uncurry (*)) xs

class Semilinear n where
    conj :: n -> n

instance {- OVERLAPS #-} Num n => Semilinear (Complex n) where
    conj = conjugate

instance {-# OVERLAPS #-} (Semilinear n, Semilinear a) => Semilinear (Linear n a) where
    conj (Linear xs) = Linear (map (\(n,x) -> (conj n, conj x)) xs)

instance {-# OVERLAPS #-} Semilinear n where
    conj = id

conj' :: Semilinear f => Linear f a -> Linear f a
conj' (Linear xs) = Linear $ map (\(n,x) -> (conj n, x)) xs

class InnerProduct a n where
    dot :: a -> a -> n

instance Num n => InnerProduct n n where
    dot = (*) . conj

instance (InnerProduct a n, Num n, Ord n, Ord a, Eq a) => InnerProduct (Linear n a) n where
    dot x y = flatten $ dot <$> conj' x <*> y

normalize :: (InnerProduct a n, Num n, Ord n, Ord a, Eq a, Floating n) => Linear n a -> Linear n a
normalize l = scale (1/sqrt (dot l l)) l

type Cplx = Complex Float
