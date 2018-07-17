{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Render(
    bitmapOfOrbital
) where

import Linear
import Gaussian
import Atom
import Orbital

import Graphics.Gloss.Data.Bitmap
import Graphics.Gloss

import qualified Data.Map as M
import qualified Data.ByteString.Lazy as BS
import Data.ByteString.Builder
import Data.List
import GHC.TypeLits

bitmapFormat :: BitmapFormat
bitmapFormat = BitmapFormat BottomToTop PxRGBA

bitmapOfOrbital :: KnownNat n => Int -> Int -> Rl -> Rl -> Linear Rl Label -> Atoms n -> Picture
bitmapOfOrbital x y scale vScale orbital atoms = bitmapOfByteString x y bitmapFormat (BS.toStrict $ bytestringOfOrbital x y scale vScale (evalOrbital atoms orbital)) True

bytestringOfOrbital :: forall n. KnownNat n => Int -> Int -> Rl -> Rl -> Gaussians n -> BS.ByteString
bytestringOfOrbital x0 y0 scale vScale orb = toLazyByteString rows
    where rows  = mconcat $ map (row  . (/scale) . fromIntegral . (+(-div y0 2))) [0..y0-1]
          row y = mconcat $ map (px y . (/scale) . fromIntegral . (+(-div x0 2))) [0..x0-1]
          px y x = colour $ evaluates orb (genericTake (natVal @n Proxy) (x:y:repeat 0))
          colour z = mconcat $ (map (word8 . floor . (*255)) (colourCode (vScale * z))) ++ [word8 255]

colourCode :: Rl -> [Rl]
colourCode z = map ((1-).(ri*).(1-).(r0*) . (\c -> if c then 1 else 0) . (==(z>0))) [True,False,False]
    where r         = abs z
          clamp     = max 0 . min 1
          softClamp = (1-) . (/4) . (^2) . (2-) . max 0 . min 2
          r0        = softClamp r
          ri        = softClamp (1/r)
