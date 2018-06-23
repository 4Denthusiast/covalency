{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeApplications #-}
module Lib(
    someFunc
) where

import Linear
import Gaussian
import Atom
import Orbital
import Render

import Graphics.Gloss
import Graphics.Gloss.Interface.IO.Interact

import qualified Data.Map as M
import Data.Monoid
import Data.Complex
import GHC.TypeLits (KnownNat)
import Debug.Trace

someFunc :: IO ()
someFunc = interactIO
    (InWindow "Covalency" (800,800) (0,0))
    white
    emptyWorld
    (return . renderWorld)
    ((return .) . handleEvent)
    (const (pure ()))

data World = World InputState (Atoms 3) [Orbital]

data InputState = Typing String | Editing String

viewScale :: Float
viewScale = 200.0

viewLOD :: Int
viewLOD = 4

emptyWorld :: World
emptyWorld = World (Typing []) M.empty [mempty]

renderWorld :: World -> Picture
renderWorld (World input atoms orbs) = pictures $ renderOrbitals atoms orbs : renderInput input : map (uncurry renderAtom) (M.toList atoms)

renderAtom :: AtomLabel -> Atom n -> Picture
renderAtom label at = let (x:y:_) = atomPos at in Color white $ Translate (x*viewScale) (y*viewScale) $ Pictures [Scale 0.1 0.1 $ Text label, Circle 20]

renderInput :: InputState -> Picture
renderInput (Typing  s) = Color red   $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('>':s)
renderInput (Editing s) = Color white $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('#':s)

renderOrbitals :: KnownNat n => Atoms n -> [Orbital] -> Picture
renderOrbitals atoms [] = Blank
renderOrbitals atoms (o:_) = Scale lodF lodF $ bitmapOfOrbital px px (viewScale / lodF) o atoms
    where px   = div 800 viewLOD
          lodF = fromIntegral viewLOD

handleEvent :: Event -> World -> World
handleEvent event w@(World inputState atoms orbs) = case inputState of
    (Typing s) -> case event of
        (EventKey (Char '\b') Down _ _) -> World (Typing (if null s then s else init s)) atoms orbs
        (EventKey (Char c) Down _ _) -> World (Typing (s++[c])) atoms orbs
        (EventKey (SpecialKey KeyEnter) Down _ _) -> World (Editing s) atoms orbs
        _ -> w
    (Editing s) -> case event of
        (EventKey (SpecialKey KeyEnter) Down _ _) -> World (Typing s) atoms orbs
        (EventKey (MouseButton LeftButton) Down _ (x,y)) -> World inputState (M.insert s (emptyAtom [x/viewScale, y/viewScale, 0]) atoms) orbs
        (EventKey (Char 'f') Down _ _) -> World inputState atoms (testOrbs atoms)
        (EventKey (Char '=') Down _ _) -> World inputState atoms (drop 40 orbs)
        (EventKey (Char '[') Down _ _) -> traceShow (head orbs) $ World inputState atoms (tail orbs ++ [head orbs])
        (EventKey (Char ']') Down _ _) -> World inputState atoms (last orbs : init orbs)
        _ -> w

testOrbs :: KnownNat n => Atoms n -> [Orbital]
testOrbs ats = {-map (return . ("",)) [-6..8]-} --map (\(i,o) -> exp (2*fromIntegral i :: Cplx) *~ fmap ("",) o) $ M.toList $ traceShowMatId $ atomKineticTerm $ ats M.! ""
  iterate (\o -> normalize @Cplx $ reduce $ o <> Linear.scale (-(0.001:+0.000)) (o >>= (nuclearHamiltonian ats M.!))) $ return ("",0)

-- Temporary for testing.
instance InnerProduct Cplx (AtomLabel, OrbitalLabel) where
    dot a b = if a == b then 1 else 0
