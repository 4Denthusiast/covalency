{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TupleSections #-}
module Linear(
    Linear(..),
    reduce,
    scale,
    flatten,
    flatten',
    Vector(..),
    towerScale,
    InnerProduct(..),
    normalize,
    mapToLinear,
    Rl,
) where

import Data.List
import qualified Data.Map as M
import Data.Monoid hiding ((<>))
import Data.Semigroup
import Control.Applicative
import Control.Monad
import Data.Bifunctor

data Linear f a = Linear [(f,a)] deriving (Ord, Eq, Show, Read)

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

flatten  :: forall f a. (Vector f a, Monoid a) => Linear f a -> a
flatten  (Linear xs) = mconcat $ map (uncurry (*~)) xs
flatten' :: forall f a. (Vector a f, Num f) => Linear f a -> f
flatten' (Linear xs) = sum $ map (uncurry $ flip (*~)) xs

class Num f => Vector f a where
    infixl 7 *~
    (*~) :: f -> a -> a

instance Num f => Vector f f where (*~) = (*)
instance Num f => Vector f (Linear f a) where (*~) = scale
instance Vector Int Rl   where n *~ x = fromIntegral n * x
instance {-# OVERLAPS #-} (Vector a b, Num b) => Vector a (Linear b x) where
    (*~) = towerScale @b
instance Vector f a => Vector f (M.Map k a) where (*~) = fmap . (*~)
instance (Vector f a, Vector f b) => Vector f (a,b) where x *~ (a,b) = (x *~ a, x *~ b)

towerScale :: forall a f b. (Num a, Vector f a, Vector a b) => f -> b -> b
towerScale f b = (f *~ (1 :: a)) *~ b

instance Semigroup Rl where (<>)   = (+)
instance Monoid    Rl where {mappend = (<>); mempty = 0}

class InnerProduct f a where
    dot :: a -> a -> f

instance Num f => InnerProduct f f where
    dot = (*)

instance (InnerProduct f a, Num f, Ord f, Monoid f) => InnerProduct f (Linear f a) where
    dot x y = flatten $ dot <$> x <*> y

normalize :: forall f a. (InnerProduct f a, Vector f a, Floating f) => a -> a
normalize l = (1/sqrt (dot @f l l)) *~ l

mapToLinear :: M.Map a f -> Linear f a
mapToLinear = Linear . map (uncurry $ flip (,)) . M.toList

type Rl = Double --In case I want to change my representation again. "Real" is taken by a Prelude typeclass.
