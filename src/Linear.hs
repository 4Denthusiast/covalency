{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
module Linear(
    Linear(..),
    lmap,
    reduce,
    scale,
    join,
    flatten,
    single,
) where

import Data.List
import Data.Monoid hiding ((<>))
import Data.Semigroup

data Linear f a where
    Linear :: (Num f, Eq f, Ord a, Eq a) => [(f,a)] -> Linear f a

lmap :: (Ord b, Eq b) => (a -> b) -> Linear f a -> Linear f b
lmap f (Linear xs) = Linear $ map (\(n,x) -> (n,f x)) xs

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
