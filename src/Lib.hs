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
import Data.Char
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
    worldOrbitals :: ([Orbital],[Orbital]),
    worldViewScale :: Float,
    worldValScale :: Rl,
    orbitalPicture :: Picture,
    worldIntegrals :: Integrals,
    worldPrevEEHamiltonian :: Matrix Label
}

data InputState = Typing String | Editing String

viewLOD :: Int
viewLOD = 4

emptyWorld :: World
emptyWorld = World (Typing []) M.empty ([],[]) 200 1 Blank undefined M.empty

num :: (Real a, Fractional b) => a -> b
num = fromRational . toRational

renderWorld :: World -> Picture
renderWorld w = pictures $
    orbitalPicture w :
    renderInput (inputState w) :
    renderAtomicNumber (inputState w) (worldAtoms w) :
    renderScales (worldViewScale w) :
    map (uncurry (renderAtom (worldViewScale w))) (M.toList (worldAtoms w))

renderAtom :: Float -> AtomLabel -> Atom n -> Picture
renderAtom viewScale label at = let (x:y:_) = atomPos at in Color white $ Translate (num x*viewScale) (num y*viewScale) $ Pictures [Scale 0.1 0.1 $ Text label, Circle 20]

renderInput :: InputState -> Picture
renderInput (Typing  s) = Color red   $ Translate (-300) (-390) $ Scale 0.2 0.2 $ Text ('>':s)
renderInput (Editing s) = Color white $ Translate (-300) (-390) $ Scale 0.2 0.2 $ Text ('#':s)

renderAtomicNumber :: InputState -> M.Map AtomLabel (Atom n) -> Picture
renderAtomicNumber i ats = Color white $ Translate (-390) (-390) $ Scale 0.2 0.2 $ Text ("Z="++z)
    where s = case i of {(Typing s) -> s; (Editing s) -> s}
          at = M.lookup s ats
          z = maybe "" (show . atomicNumber) at

renderScales :: Float -> Picture
renderScales vs = Color white $ Pictures [bar, t]
    where bar = Line [(-390,390),(-390,370),(-390,380),(-190,380),(-190,370),(-190,390)]
          t   = Translate (-320) (360) $ Scale 0.15 0.15 $ Text $ show (200/vs) ++ "au"

reRender :: World -> World
reRender w = w{orbitalPicture = renderOrbitals (worldViewScale w) (worldValScale w) (worldAtoms w) (fst $ worldOrbitals w)}

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
        (EventKey (Char 'r') Down _ _) -> w{worldIntegrals = calculateIntegrals atoms}
        (EventKey (Char 'R') Down _ _) -> rr w{worldOrbitals = ([],[]), worldPrevEEHamiltonian = M.empty}
        (EventKey (Char 'h') Down _ _) -> rr w{worldOrbitals = (hydrogenLikeOrbs ints,[])}
        (EventKey (Char n) Down _ _) | isDigit n -> rr $ hfStepWorld (0.5 ^ digitToInt n) w
        (EventKey (Char '[') Down _ _) -> traceShow (take 1 $ fst orbs) $ rr w{worldOrbitals = shiftRight orbs}
        (EventKey (Char ']') Down _ _) -> rr w{worldOrbitals = shiftLeft orbs}
        (EventKey (SpecialKey KeyDelete) Down _ _) -> rr w{worldOrbitals = ([],[]), worldAtoms = M.delete s atoms}
        (EventKey (Char '+') Down _ _) -> w{worldAtoms = M.adjust (\a -> changeZ (atomicNumber a + 1) a) s atoms}
        (EventKey (Char '-') Down _ _) -> w{worldAtoms = M.adjust (\a -> changeZ (max 1 $ atomicNumber a - 1) a) s atoms}
        _ -> both s
    where atoms = worldAtoms w
          orbs  = worldOrbitals w
          viewScale = worldViewScale w
          valScale  = worldValScale w
          ints  = worldIntegrals w
          rr    = reRender
          d     = 3 --TODO explicitly link this to the dimension of the atoms.
          pos (x,y) = [num $ x/viewScale, num $ y/viewScale, 0]
          both s = case event of
              (EventKey (MouseButton LeftButton) Down _ p) -> rr w{worldAtoms = M.alter (Just . maybe (newAtom 1 (pos p)) (\a -> a{atomPos = pos p})) s atoms}
              (EventKey (MouseButton WheelDown) Down _ _) -> rr w{worldViewScale = viewScale / 2, worldValScale = valScale * (2 ** (d/2))}
              (EventKey (MouseButton WheelUp  ) Down _ _) -> rr w{worldViewScale = viewScale * 2, worldValScale = valScale / (2 ** (d/2))} 
              _ -> w

shiftLeft :: ([a],[a]) -> ([a],[a])
shiftLeft ([],[]) = ([],[])
shiftLeft (xs,[])  = ([last xs], reverse (init xs))
shiftLeft (xs',x:xs) = (x:xs',xs)

shiftRight :: ([a],[a]) -> ([a],[a])
shiftRight ( [x],xs') = (reverse (x:xs'),[])
shiftRight (x:xs,xs') = (xs,x:xs')
shiftRight (  [], []) = ([],[])

hydrogenLikeOrbs :: Integrals -> [Orbital]
hydrogenLikeOrbs (_,h,_) = negativeEigenvecs h

hfStepWorld :: Rl -> World -> World
hfStepWorld s w = w{worldOrbitals = (orbs',[]), worldPrevEEHamiltonian = peeh'}
    where (peeh', orbs') = hartreeFockStep s (div n 2) (worldIntegrals w) (worldPrevEEHamiltonian w, (\(o,o') -> reverse o' ++ o) $ worldOrbitals w)
          n = sum $ atomicNumber <$> worldAtoms w
