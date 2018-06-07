{-# LANGUAGE TupleSections #-}

module Lib(
    someFunc
) where

import Gaussian
import Atom
import Orbital
import Render

import Graphics.Gloss
import Graphics.Gloss.Interface.IO.Interact

import qualified Data.Map as M
import Debug.Trace

someFunc :: IO ()
someFunc = interactIO
    (InWindow "Covalency" (800,800) (0,0))
    white
    emptyWorld
    (return . renderWorld)
    ((return .) . handleEvent)
    (const (pure ()))

data World = World InputState (M.Map AtomLabel Atom) [Orbital]

data InputState = Typing String | Editing String

viewScale :: Float
viewScale = 400.0

viewLOD :: Int
viewLOD = 4

emptyWorld :: World
emptyWorld = World (Typing []) M.empty [Orbital []]

renderWorld :: World -> Picture
renderWorld (World input atoms orbs) = pictures $ renderOrbitals atoms orbs : renderInput input : map (uncurry renderAtom) (M.toList atoms)

renderAtom :: AtomLabel -> Atom -> Picture
renderAtom label (Atom x y _) = Color white $ Translate (x*viewScale) (y*viewScale) $ Pictures [Scale 0.1 0.1 $ Text label, Circle 20]

renderInput :: InputState -> Picture
renderInput (Typing  s) = Color red   $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('>':s)
renderInput (Editing s) = Color white $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('#':s)

renderOrbitals :: Atoms -> [Orbital] -> Picture
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
        (EventKey (MouseButton LeftButton) Down _ (x,y)) -> World inputState (M.insert s (emptyAtom (x/viewScale) (y/viewScale) (length atoms)) atoms) orbs
        (EventKey (Char 'f') Down _ _) -> World inputState atoms (testOrbs atoms)
        _ -> w

testOrbs :: Atoms -> [Orbital]
testOrbs = (:[]) . Orbital . map (,(),1) . M.keys
