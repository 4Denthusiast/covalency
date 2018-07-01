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
import Eigen

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
    black
    emptyWorld
    (return . renderWorld)
    ((return .) . handleEvent)
    (const (pure ()))

data World = World{
    inputState :: InputState,
    worldAtoms :: Atoms 3,
    worldOrbitals :: [Orbital],
    worldViewScale :: Float,
    worldValScale :: Rl,
    orbitalPicture :: Picture
}

data InputState = Typing String | Editing String

viewLOD :: Int
viewLOD = 4

emptyWorld :: World
emptyWorld = World (Typing []) M.empty [] 200 1 Blank

num :: (Real a, Fractional b) => a -> b
num = fromRational . toRational

renderWorld :: World -> Picture
renderWorld w = pictures $
    orbitalPicture w :
    renderInput (inputState w) :
    renderScales (worldViewScale w) :
    map (uncurry (renderAtom (worldViewScale w))) (M.toList (worldAtoms w))

renderAtom :: Float -> AtomLabel -> Atom n -> Picture
renderAtom viewScale label at = let (x:y:_) = atomPos at in Color white $ Translate (num x*viewScale) (num y*viewScale) $ Pictures [Scale 0.1 0.1 $ Text label, Circle 20]

renderInput :: InputState -> Picture
renderInput (Typing  s) = Color red   $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('>':s)
renderInput (Editing s) = Color white $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ('#':s)

renderScales :: Float -> Picture
renderScales vs = Color white $ Pictures [bar, t]
    where bar = Line [(-390,390),(-390,370),(-390,380),(-190,380),(-190,370),(-190,390)]
          t   = Translate (-320) (360) $ Scale 0.15 0.15 $ Text $ show (200/vs) ++ "au"

reRender :: World -> World
reRender w = w{orbitalPicture = renderOrbitals (worldViewScale w) (worldValScale w) (worldAtoms w) (worldOrbitals w)}

renderOrbitals :: KnownNat n => Float -> Rl -> Atoms n -> [Orbital] -> Picture
renderOrbitals _ _ _ [] = Blank
renderOrbitals viewScale valScale atoms (o:_) = Scale lodF lodF $ bitmapOfOrbital px px (num viewScale / lodF) valScale o atoms
    where px   = div 800 viewLOD
          lodF :: Num a => a
          lodF = fromIntegral viewLOD

handleEvent :: Event -> World -> World
handleEvent event w = case inputState w of
    (Typing s) -> case event of
        (EventKey (Char '\b') Down _ _) -> w{inputState = Typing (if null s then s else init s)}
        (EventKey (Char c) Down _ _) -> w{inputState = Typing (s++[c])}
        (EventKey (SpecialKey KeyEnter) Down _ _) -> w{inputState = Editing s}
        _ -> both s
    (Editing s) -> case event of
        (EventKey (SpecialKey KeyEnter) Down _ _) -> w{inputState = Typing s}
        (EventKey (Char 'f') Down _ _) -> r w{worldOrbitals = testOrbs atoms}
        (EventKey (Char '=') Down _ _) -> r w{worldOrbitals = drop 40 orbs}
        (EventKey (Char '[') Down _ _) -> traceShow (head orbs) $ r w{worldOrbitals = tail orbs ++ [head orbs]}
        (EventKey (Char ']') Down _ _) -> r w{worldOrbitals = last orbs : init orbs}
        (EventKey (SpecialKey KeyDelete) Down _ _) -> r w{worldOrbitals = [], worldAtoms = M.delete s atoms}
        _ -> both s
    where atoms = worldAtoms w
          orbs  = worldOrbitals w
          viewScale = worldViewScale w
          valScale  = worldValScale w
          r     = reRender
          d     = 3 --TODO explicitly link this to the dimension of the atoms.
          both s = case event of
              (EventKey (MouseButton LeftButton) Down _ (x,y)) -> r w{worldAtoms = M.insert s (emptyAtom [num $ x/viewScale, num $ y/viewScale, 0]) atoms}
              (EventKey (MouseButton WheelDown) Down _ _) -> r w{worldViewScale = viewScale / 2, worldValScale = valScale * (2 ** (d/2))}
              (EventKey (MouseButton WheelUp  ) Down _ _) -> r w{worldViewScale = viewScale * 2, worldValScale = valScale / (2 ** (d/2))} 
              _ -> w

testOrbs :: KnownNat n => Atoms n -> [Orbital]
testOrbs ats =
    --hartreeFockIterants ats 1 !! 20
    --map (return.("",)) $ M.keys $ atomOrbitals $ ats M.! ""
    negativeEigenvecs (nuclearHamiltonian ats)
    --M.elems $ nuclearHamiltonian ats
