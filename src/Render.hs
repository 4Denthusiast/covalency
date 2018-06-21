{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Render(
    bitmapOfOrbital
) where

import Linear (Cplx)
import Gaussian
import Atom
import Orbital

import Graphics.Gloss.Data.Bitmap
import Graphics.Gloss

import qualified Data.Map as M
import qualified Data.ByteString.Lazy as BS
import Data.ByteString.Builder
import Data.List
import Data.Complex
import GHC.TypeLits

bitmapFormat :: BitmapFormat
bitmapFormat = BitmapFormat BottomToTop PxRGBA

bitmapOfOrbital :: KnownNat n => Int -> Int -> Float -> Orbital -> Atoms n -> Picture
bitmapOfOrbital x y scale orbital atoms = bitmapOfByteString x y bitmapFormat (BS.toStrict $ bytestringOfOrbital x y scale orbital atoms) True

bytestringOfOrbital :: forall n. KnownNat n => Int -> Int -> Float -> Orbital -> Atoms n -> BS.ByteString
bytestringOfOrbital x0 y0 scale orb atoms = toLazyByteString rows
    where rows  = mconcat $ map (row  . (/scale) . fromIntegral . (+(-div y0 2))) [0..y0-1]
          row y = mconcat $ map (px y . (/scale) . fromIntegral . (+(-div x0 2))) [0..x0-1]
          px y x = colour $ evalOrbital atoms orb (genericTake (natVal @n Proxy) (x:y:repeat 0))
          colour z = mconcat $ (map (word8 . floor . (*255)) (colourCode z)) ++ [word8 255]

colourCode :: Cplx -> [Float]
colourCode z = map ((1-).(ri*).(1-).(r0*) . softClamp . (2-) . abs . (/(pi/3)) . rerange . (θ+)) [0, 2*pi/3, 4*pi/3]
    where (r, θ)    = polar z
          rerange t = if t > pi then t - 2*pi else t
          clamp     = max 0 . min 1
          softClamp = (1-) . (/4) . (^2) . (2-) . max 0 . min 2
          r0        = softClamp r
          ri        = softClamp (1/r)
