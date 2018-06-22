{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleInstances #-}
import Gaussian
import Polynomial (Polynomial, constant, variable)
import qualified Polynomial
import Linear

import System.Random
import Data.Complex
import Control.Monad
import GHC.TypeLits

main :: IO ()
main = testGaussians

testGaussians :: IO ()
testGaussians = sequence_ $ replicate 5 $ randTestConvolve >> testMultiply

randTestConvolve :: IO ()
randTestConvolve = testConvolve <$> randomGaussian <*> randomGaussian <>> randomPoint

testConvolve :: Gaussian 4 -> Gaussian 4 -> [Float] -> IO ()
testConvolve a b x = let
        a' = shiftGauss x (reverseGauss a)
        multipliedX = integral $ multiply a' b
        convolvedX  = evaluate (convolve a b) x
    in
        (if abs (multipliedX - convolvedX) > 0.001 then
            putStr "ERR"
        else
            putStr "   "
        )
        >> putStrLn (" convolve: "++show multipliedX++", "++show convolvedX)

testMultiply :: IO ()
testMultiply = do
    a <- randomGaussian
    b <- randomGaussian
    x <- randomPoint
    let distX = evaluate a x * evaluate b x
        multX = evaluate (multiply a b) x
    if abs (distX - multX) > 0.001 then
        putStr "ERR"
    else
        putStr "   "
    putStrLn (" multiply: "++show distX++", "++show multX)

reverseGauss :: KnownNat n => Gaussian n -> Gaussian n
reverseGauss (Gaussian xs c a) = Gaussian (map negate xs) c (Polynomial.evaluate (map (negate . variable) [0..]) a)

randomGaussian :: IO (Gaussian 4)
randomGaussian = Gaussian <$> randomPoint <*> ((1+) <$> randomIO) <*> randomPoly

randomPoint = sequence $ replicate 4 randomIO

randomPoly :: IO (Polynomial 4 Cplx)
randomPoly = product <$> sequence (replicate 5 (randomFactor 4))
    where randomFactor :: Int -> IO (Polynomial 4 Cplx)
          randomFactor 0 = constant <$> randomIO
          randomFactor n = (+) <$> randomFactor (n-1) <*> (flip (*~) (variable n) <$> (randomIO :: IO Cplx))

instance Random (Complex Float) where
    randomR ((r0:+i0),(r1:+i1)) g = let (r,g') = randomR (r0,r1) g in let (i,g'') = randomR (i0,i1) g' in (r:+i,g'')
    random = randomR (-(1:+1),1:+1)

infixl 4 <>>
(<>>) :: Monad m => m (a -> m b) -> m a -> m b
f <>> x = join $ f <*> x
