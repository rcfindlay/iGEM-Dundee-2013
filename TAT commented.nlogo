;;;;Dundee iGEM Team 2013;;;;
;;;;;;;Dry Team;;;;;;;;;;;;
;;NasRach Collaborations;;




;creation of populations;
breed [PP1s PP1]                                                           ;; red PP1 molecules
breed [degPs degP]                                                         ;; degradation machinery - black squares
breed [mcs mc]                                                             ;; microcystin, purple turtles, assuming so small that kinetic interaction is negligible
breed [complexes complex]                                                  ;; PP1 bounded to microcystin, purple crosses
breed [tatAs tatA]                                                         ;; tatA proteins for gate use



;set up of all variables;
globals [  
  seperation                                                               ;; controls distance between gates
  
  
  ;variables required for loops;
  j                                                                        ;; keeps count of how  many gates being created
  k                                                                        ;; counts rows when colouring periplasm
  l                                                                        ;; counts columns when colouring periplasm
  g                                                                        ;; keeps count of gate width when colouring
  
  deg-number                                                               ;; a count of the number of degraded PP1/complexes
  PP1-transported                                                          ;; a count of PP1 transported from cytoplasm to periplasm
]



;attributes associated with all turtles;
turtles-own 
[
  speed mass energy                                                        ;; controls the speed, mass & energy of the turtles for kinetics
  last-collision                                                           ;; keeps note of when the previous collision has occured
  transport                                                                ;; a boolean variable keeping track of whether a PP1 is in the transport stage requiring the binding of tatAs
  bound-tatA                                                               ;; a counter to keep a count of the number of tatAs bound to a given transporting PP1
]













;creates a class called particles for calling-in PP1s, degPs, mcs & complexes;
to-report particles
  report (turtle-set PP1s degPs mcs complexes)
end



;creates a separate class for particles which exist in the periplasm;
to-report periparticles
  report (turtle-set degPs complexes)
end













;setup for the characteristics of PP1 proteins;
to setup-PP1
  set mass 37500                                                           ;; the masses used was the approximate mass of PP1 (all masses are in approximations of daltons)
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  
  ;physical characteristics;
  set color red
  set size 1.5
  
  ;setting up the initial transport conditions;
  set transport false                                                      
  set bound-tatA 0                                                         
end


;setup for the characteristics of degradation machinery;
to setup-degPs 
  set speed 10   
  set mass 37500                                                           ;; the mass of our degrading machinery has been taken as that of PP1 due to a lack of information. This means it will act in a similar manner to PP1
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  set size 2
  set color black
end


;setup for the characteristics of microcystin;
to setup-mcs
  set speed 10
  set mass 995                                                             ;; the mass of microcystin-LR is approx 995 daltons
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  set size 1
  set color 114
end


;setup for the characteristics of complexes (PP1-microcystin bound together);
to setup-complex
  set speed 10
  set mass 38495                                                           ;; the mass of the complex was simply defined as that of microcystin + PP1
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  set shape "x"
  set size 3
  set color 126
end

;setup for the characteristics of tatA;
to setup-tatA
  set speed 10
  set mass 37700                                                             ;; the mass of the tatA protein is ~37.7kDa
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  set shape "triangle 2"
  set size 1
  set color 45
  set transport false
end












;The setup procedure creates our environment and sets some conditions for the simulation;
to setup [both-sides?]
  clear-all
  
  ;sets the shapes of particles;
  set-default-shape PP1s "circle"                                          ;; PP1 is depicted as a red circle 
  set-default-shape degPs "square"                                         ;; degradation machinery is depicted as a black square
  set-default-shape mcs "turtle"                                           ;; microcystin is depicted as a red circle
  set-default-shape tatAs "triangle 2"                                     ;; tatA is depicted as a yellow triangle


  
  ;sets up the initial conditions;
  
  ;beginning loop values;
  set k ( 99 )
  set l ( - 99 )
  set g 0
  set j 0
  
  ;beginning with a clean slate in degradated number, transported number and gate separation;
  set deg-number 0
  set PP1-transported 0
  set seperation 0
 
 
 
 
 
 
  ;calls for the creation of the environment with a membrane and gate;
  walls
  gates
  set j 0
  set g 0
  pores
  

  
  ;creates an initial number of PP1 molecules with random positions in the cytoplasm;
  create-PP1s initial-PP1
  [setup-PP1
    randomize-right ]


  ;creates an initial number of tatA molecules with random positions in the inner membrane;
  create-tatAs tatA-number
  [ setup-tatA
    randomize-inmem ]
  
  
  ;creates an initial number of degradation machineries and microcystin numbers with random positions in the periplasm;  
  create-degPs degP-number
  [setup-degPs
    randomize-periplasm ]
    
  create-mcs mc-number
  [ setup-mcs
    randomize-left]
  

    
 
    
  reset-ticks
end







;position randomizing procedure for extracellular particles (microcystin);
to randomize-left
  setxy (- (abs random-xcor))
  random-ycor
  if any? patches in-radius 1 with [pcolor != 8]
    [ randomize-left ] ;;try again until we don't land outside the cell's outer wall
end




;position randomizing procedure for periplasmic proteins;
to randomize-periplasm
  setxy  random-xcor  random-ycor
  if any? patches in-radius 1 with [pcolor != 87]
    [ randomize-periplasm ] ;;try again until we don't land on or near the edges of the periplasm
end


;position randomizing procedure for tatAs in the inner membrane;
to randomize-inmem
  setxy  random-xcor  random-ycor
  if ([pxcor] of patch-here > 25 or [pxcor] of patch-here < 15)
    [ randomize-inmem ] ;;try again until we don't land within the inner membrane
end


;position randomizing procedure for cytoplasmic proteins;
to randomize-right
  setxy abs random-xcor
  random-ycor
  if any? patches in-radius 1 with [pcolor != 19.9]
    [ randomize-right ] ;;try again until we don't land on or near the edges of the cytoplasm
end
















;this procedure is what is run through (looping) during the operation of the simulation;
to go

    
  ;in dealing with proteins which return into unseen parts of the cell and microcystin which, after entering the periplasm, leaves the cell we must introduce losses at the boundaries;
  ask PP1s [
    loss-PP1                                                               ;; calls a procedure to deal with PP1 molecules which return into the cell
  ]
    
  ask mcs [
    loss-mc                                                                ;; calls a procedure to deal with microcystin molecules which exit the cell
    bind                                                                   ;; calls a procedure to deal with the binding of microcystin to PP1 and so the formation of complexes
  ]
  
  
  ;to deal with the production rate of PP1 and the introduction rate of microcystin, a set number of proteins and molecules are created;
  ;these are created on every fifth tick to ensure a steady rate of creation and are controlled by sliders in the interface;
  if ticks mod 5 = 0 [
    create-PP1s PP1-production
    [
      setxy abs 48
      random-ycor
      setup-PP1
    ]
    
    create-mcs mc-production
    [
      setxy ( - 48 )
      random-ycor
      setup-mcs
    ]
  ]
 
 
  degrade                                                                  ;; a procedure to degrade PP1s and complexes which encounter the degradation machinery  
  
  
  
  ;for remaining memebers of particles class, carry out the usual checks for collision and motion before carrying out any required bounces and then progress
  ask particles [
    if collide? [check-for-collision]
    if random-walk? [rt random-float 360]
    bounce
    if breed != PP1s [fd 1]
  ]
  
  ;for tatA proteins, they must only continue moving through the inner membrane if they have not been captured to aid in the transportation of a PP1 molecule;
  ask tatAs [
    if transport = false [
    if collide? [check-for-collision]
    if random-walk? [rt random-float 360]
    stick
    bounce
    fd 1 ] 
    
  ]
  
  ;only free PP1 proteins which are not being transported can continue moving;
  ask PP1s [
    if transport = false [fd 1]]
  
  ;if a PP1 molecule which is being transported gains a sufficiently high number of tatA proteins, it is transported across the inner membrane;
  ask PP1s [
    if bound-tatA = tatA-required [translocation]]
  
  
  
  ;the simulation continues forward one tick;
  tick
end


















;;;;Procedures to create the walls and gates of chosen width;;;;

;procedure used to create the environment with coloured areas representing the cytoplasm and membrane;
to walls
  ask patches with [pxcor > 25 and pxcor < 51][ set pcolor 19.9]            ;; colour the cytoplasm (right half) white                      
  ask patches with [ pxcor < 26 and  pxcor > 14 ][ set pcolor 75 ]          ;; colour the membrane green      
  ask patches with [ pxcor < ( - 29 ) ] [ set pcolor 8 ]                    ;; colour the region indicating the rest of the cell (right) light grey
  ask patches with [ pxcor > -19 and pxcor < 15] [ set pcolor 87 ]          ;; colour the periplasm (left side) light blue
  ask patches with [ pxcor > (- 30) and pxcor < ( - 18 ) ] [ set pcolor 3 ] ;; colour the outermembrane (left) dark grey  
end



;Deals with the creation of pores in the outer membrane;
to gates
 
  ;for an odd number of gates;
  if (gate-number mod 2 ) = 1 [                                             
    set seperation ( round ( 99 / ( (gate-number + 1) * 2 ) ) )            ;; calculate the required separation between consecutive gate centres (odd number of gates)
    if j != ((gate-number + 1) / 2) [                                      ;; run through the pairs of gates (above and below the centre) using j as our counting variable and 
                                                                           ;;only stop once all gates have been considered
      
      ask patches with [ (pxcor > (14) and pxcor < ( 26 )) and 
        (pycor = 0 + (seperation * j * 2) or 
          pycor = 0 - (seperation * j * 2) )] [                            ;; ask patches in the centre of our current gate pair (j)
      set pcolor 78                                                        ;; colour these patches a light green
      set g 0                                                              ;; set the gate-size count of the current gate pair (g)  = 0              
      gate-width                                                           ;; call the gate-width procedure to colour the width of the gate this same light green colour
          ]
      
      set j (j + 1)                                                        ;; progress to the next set of gate pairs
      gates                                                                ;; re-run gate until all gate pairs have been considered
    ]
  ]
  
  
  ;for an even number of gates;
  if (gate-number mod 2) = 0 [
    set seperation ( round ( 99 / ( gate-number * 2 ) ) )                  ;; calculate the required separation between consecutive gate centres (even number of gates)
    if j != (gate-number / 2) [                                            ;; run through the pairs of gates (above and below the centre) using j as our counting variable and
                                                                           ;; only stop once all gates have been considered
      ask patches with [ (pxcor > (14) and pxcor < ( 26 )) and 
        ( pycor = ( - seperation) + (seperation * (j + 1) * 2) or       
          pycor = seperation - (seperation * (j + 1) * 2) )] [             ;; ask patches in the centre of our current gate pair (j)
      set pcolor 78                                                        ;; colour these patches a light green
      set g 0                                                              ;; set the gate-size count of the current gate pair (g)  = 0
      gate-width                                                           ;; call the gate-width procedure to colour the width of the gate this same light green colour
          ] 
      
      set j (j + 1)                                                        ;; progress to the next set of gate pairs
      gates                                                                ;; re-run gate until all gate pairs have been considered
    ]
  ]
    
    
end



; procedure to alter the sizes of the gates to the specified widths;
to gate-width
  
  ;for odd gate-sizes;
  if ( gate-size mod 2) = 1 [                                                  
    if g !=(( gate-size + 1 ) / 2  ) [                                     ;; running through g up to half the width of the required gate 
      ask patch-at-heading-and-distance 0 g [set pcolor 78]                ;; for this gate, set the colour of the patch at an incrementing co-ordinate (g) above = light green
      ask patch-at-heading-and-distance 180 g [set pcolor 78]              ;; for this gate, set the colour of the patch at an incrementing co-ordinate (g) below = light green
      set g ( g + 1 )                                                      ;; increment g by 1 to make the gate wider if required
      gate-width                                                           ;; re-run this script until the current gate (j) is of the correct width
    ]
  ]
  
  ;for odd gate-sizes;
  if ( gate-size mod 2) = 0 [
    if g !=(( gate-size ) / 2  ) [                                         ;; running through g up to half the width of the required gate
      ask patch-at-heading-and-distance 0 1 [set pcolor 78]                ;; due to a lack of symmetry for an even gate width case, we must push the patches up by 1
      ask patch-at-heading-and-distance 0 (g + 1) [set pcolor 78]          ;; for this gate, set the colour of the patch at an incrementing co-ordinate (g) above = light green
      ask patch-at-heading-and-distance 180 g [set pcolor 78]              ;; for this gate, set the colour of the patch at an incrementing co-ordinate (g) below = light green
      set g ( g + 1 )                                                      ;; increment g by 1 to make the gate wider if required
      gate-width                                                           ;; re-run this script until the current gate (j) is of the correct width
      
    ]    
  ]
end





;Deals with the creation of pores in the outer membrane;
to pores
 
  ;for an odd number of pores;
  if (pore-number mod 2 ) = 1 [                                             
    set seperation ( round ( 99 / ( (pore-number + 1) * 2 ) ) )            ;; calculate the required separation between consecutive pore centres (odd number of pores)
    if j != ((pore-number + 1) / 2) [                                      ;; run through the pairs of pores (above and below the centre) using j as our counting variable and 
                                                                           ;; only stop once all pores have been considered
      
      ask patches with [ (pxcor > (- 30) and pxcor < ( - 18 )) and 
        (pycor = 0 + (seperation * j * 2) or 
          pycor = 0 - (seperation * j * 2) )] [                            ;; ask patches in the positions of our current pore pair (j)
      set pcolor 8                                                         ;; colour these patches a light grey
      set g 0                                                              ;; set the pore-size count of the current pore pair (g)  = 0              
      pore-width                                                           ;; call the pore-width procedure to colour the width of the gate this same light grey colour
          ]
      
      set j (j + 1)                                                        ;; progress to the next set of pore pairs
      pores                                                                ;; re-run "pores" until all the required pores have been created
    ]
  ]
  
  
  ;for an even number of pores;
  if (pore-number mod 2) = 0 [
    set seperation ( round ( 99 / ( pore-number * 2 ) ) )                  ;; calculate the required separation between consecutive pore centres (even number of pores)
    if j != (pore-number / 2) [                                            ;; run through the pairs of pores (above and below the centre) using j as our counting variable and
                                                                           ;; only stop once all gates have been considered
      ask patches with [ ((pxcor > (- 30) and pxcor < ( - 18 ))) and 
        ( pycor = ( - seperation) + (seperation * (j + 1) * 2) or       
          pycor = seperation - (seperation * (j + 1) * 2) )] [             ;; ask patches at the positions of our current pore pair (j)
      set pcolor 8                                                         ;; colour these patches a light grey colour
      set g 0                                                              ;; set the pore-size count of the current pore pair (g)  = 0
      pore-width                                                           ;; call the pore-width procedure to colour the width of the gate this same light grey colour
          ] 
      
      set j (j + 1)                                                        ;; progress to the next set of pore pairs
      pores                                                                ;; re-run "pores" until all gate pairs have been considered
    ]
  ]
    
    
end



; procedure to alter the sizes of the pores to the specified widths;
to pore-width
  
  ;for odd pore-sizes;
  if ( pore-size mod 2) = 1 [                                                  
    if g !=(( pore-size + 1 ) / 2  ) [                                     ;; running through g up to half the width of the required pore 
      ask patch-at-heading-and-distance 0 g [set pcolor 8]                 ;; for this pore, set the colour of the patch at an incrementing co-ordinate (g) above = light grey
      ask patch-at-heading-and-distance 180 g [set pcolor 8]               ;; for this pore, set the colour of the patch at an incrementing co-ordinate (g) below = light grey
      set g ( g + 1 )                                                      ;; increment g by 1 to make the gate wider if required
      pore-width                                                           ;; re-run this script until the current gate (j) is of the correct width
    ]
  ]
  
  ;for odd pore-sizes;
  if ( pore-size mod 2) = 0 [
    if g !=(( pore-size ) / 2  ) [                                         ;; running through g up to half the width of the required pore
      ask patch-at-heading-and-distance 0 1 [set pcolor 8]                 ;; due to a lack of symmetry for an even gate width case, we must push the patches up by 1
      ask patch-at-heading-and-distance 0 (g + 1) [set pcolor 8]           ;; for this pore, set the colour of the patch at an incrementing co-ordinate (g) above = light grey
      ask patch-at-heading-and-distance 180 g [set pcolor 8]               ;; for this gate, set the colour of the patch at an incrementing co-ordinate (g) below = light grey
      set g ( g + 1 )                                                      ;; increment g by 1 to make the pore wider if required
      pore-width                                                           ;; re-run this script until the current pore (j) is of the correct width
      
    ]    
  ]
end
    
   













;;;;Procedures for protein/particle death and loss;;;;


; procedure to kill PP1 proteins which disappear back into the cell;
to loss-PP1
  if [pxcor] of patch-ahead 1 = 50 [                                       ;; if the patch ahead is that which goes to the the inner cell (right wall - colour light grey))
    die                                                                    ;; if the PP1 is not linked to a tat gate, it returns to the cell and dissapears (dies)
    ]
  
end




; procedure to kill microcystin which disappear back out into the environment;
to loss-mc 
  if [pcolor] of patch-ahead 2 = 19.9 [die]                                   ;; if the patch ahead is that of the outer membrane (dark grey) microcystin escape and so die
end




; procedure to degrade proteins in the periplasm using degradation proteins;
to degrade
  
  ask degPs [ let prey one-of PP1s-here                                    ;; ask degradation machinery to target PP1s on their patch
    if prey != nobody [                                                    ;; if a PP1 exists on their patch
      ask PP1s-here [ 
        
        if random 100 <= degprob [                                         ;; depending upon the degradation probability defined in the interface
          ask prey [die]                                                   ;; have the PP1 degrade and die
        ]
      ]set deg-number ( deg-number + 1 )]                                  ;; increase the degraded tally by one
  ]
  
  
  ask degPs [ let prey one-of complexes-here                               ;; similarly ask the degradation machinery to target complexes on their patch
    if prey != nobody [                                                    ;; if there is such a complex
      ask complexes-here [ 
        if random 100 <= degprob [                                         ;; depending on the degradation probability
          ask prey [die]                                                   ;; have the complexes become degraded and so die
        ]
      ]set deg-number ( deg-number + 1 )]                                  ;; increase the degraded tally by one
  ]
  
  
end















;procedure dealing with the binding of microcystin to PP1 in the periplasm
to bind
  
  if count PP1s-here = 1 [                                                 ;; ask tatAs to check if there is a single PP1 on their current position
    if random 100 <= bindprob [                                            ;; continue depending upon the binding probability set in the interface
      ask PP1s-here [
        set breed complexes                                                ;; have the PP1 convert into a differenct agent type - a PP1-microcystin complex
        setup-complex                                                      ;; set the characteristics of this complex
        ]
      
      die]                                                                 ;; have the original microcystin disappear (die)
  ]
  
end
    





;procedure controling the bouncing of particles from walls assuming perfect reflectivity
to bounce
  
  ;for all particle types other than microcystin and tatA, always bounce off the solid parts of the inner and outer membrane
  if breed != mcs and breed != tatAs [
    if [pcolor] of patch-ahead 1 = 75 or [pcolor] of patch-ahead 1 = 3 [
      set heading (- heading) 
    ]
  ]
  
  ;specifically for microcystin, a more complex bounce command is needed for bouncing off the outer membrane;
  ;this complex bounce command ensures that when the microcystin bounces off a horizontal surface, it is correctly reflected;
  if breed = mcs [
    if [pcolor] of patch-here = 8 and [pcolor] of patch-ahead 1 = 3 [
      ifelse (xcor > (- 29) and xcor < (-17)) [
          set heading (180 - heading)
      ]
      [
        set heading (- heading)
      ]
    ]
    ;after entry of microcystin into the periplasm, it is no longer allowed to leave;
    if [pcolor] of patch-here = 87 and [pcolor] of patch-ahead 1 = 3 [
      set heading (- heading)
    ]
    ;microcystin cannot penetrate the inner membrane and must bounce here;
    if [pcolor] of patch-ahead 1 = 75 [set heading (- heading)]
  ]
  
  
  
  ;degradation machinery and PP1-microcystin complexes only reside within the periplasm and so should reflect from any surfaces;
  if breed = complexes or breed = degPs [
    if [pcolor] of patch-ahead 1 != 87 [set heading (- heading)]]
  
  ;tatA proteins must remain within the inner membrane and so reflect when entering the cytoplasm or periplasm;
  if breed = tatAs [
    if [pcolor] of patch-ahead 1 = 19.9  [set heading (- heading)]
    if [pcolor] of patch-ahead 1 = 87[set heading (- heading)] ]
  
  ;microcystin should reflect from gates also, they cannot enter a gate and be transported to the cytoplasm;
  if breed = mcs [
    if [pcolor] of patch-ahead 1 = 78  [set heading (- heading)]]
  
  
  if breed = PP1s [
    ;if PP1 molecules find themselves at an available gate for transport, they must set their transport true and await tatA proteins;
    if [pcolor] of patch-here = 78 [set transport true transporter]
    ;PP1 molecules in the periplasm cannot enter pores and exit the cell;
    if [pcolor] of patch-ahead 1 = 8 [set heading (- heading)]
    ;PP1 molecules in the periplasm can not be transported back to the cytoplasm, they will instead by reflected at the gates;
    if [pcolor] of patch-here = 87 and [pcolor] of patch-ahead 1 = 78 [set heading (- heading)]
    ]
    
  ;no particle other than tatA proteins are allowed to enter an active gate region, they are reflected;
  if breed = PP1s or breed = complexes or breed = degPs or breed = mcs [
    if [pcolor] of patch-ahead 1 = 129 [set heading (- heading)]]
    
    
end




;when a PP1 molecule encounters a free gate in the periplasm, the gate is made active and the colour of the gate altered from light green to light pink;
to transporter
  ask patches in-radius 1 with [pcolor = 78][
    set pcolor 129 ]
  repeat (gate-size * 12) [
    ask patches with [pcolor = 129] [
     ask patches in-radius 1 with [pcolor = 78] [
      set pcolor 129]]]
end








; a procedure detecting collisions between particles in the environment
to check-for-collision 
  if count other turtles-here = 1                                          ;; if two turtles occupy the same current space    
  [
    let candidate one-of other turtles-here with                           ;; choose one as a candidate
      [who < [who] of myself and myself != last-collision]                 ;; who was not the last candidate for collision
    
    if (candidate != nobody) and (speed > 0 or [speed] of candidate > 0)   ;; if there is a candidate who is not stationary
    [
      collide-with candidate                                               ;; run the collide-with procedure with the candidate
      set last-collision candidate                                         ;; having taken part in a collision, set the asking turtle's last collision
      ask candidate [ set last-collision myself ]                          ;; and the candidate turtle's last collision
    ]
  ]
  
  
end





; a procedure controlling the physics of collisions between particles which are colliding
; the equations below take into account the masses, speeds and heading of the colliding particles and carries out momentum conservation calculations
to collide-with [ other-particle ] 
  let mass2 [mass] of other-particle
  let speed2 [speed] of other-particle
  let heading2 [heading] of other-particle
  let theta (random-float 360)
  let v1t (speed * cos (theta - heading))
  let v1l (speed * sin (theta - heading))
  let v2t (speed2 * cos (theta - heading2))
  let v2l (speed2 * sin (theta - heading2))
  let vcm (((mass * v1t) + (mass2 * v2t)) / (mass + mass2) )
  set v1t (2 * vcm - v1t)
  set v2t (2 * vcm - v2t)
  set speed sqrt ((v1t ^ 2) + (v1l ^ 2))
  set energy (0.5 * mass * speed ^ 2)
  if v1l != 0 or v1t != 0
    [ set heading (theta - (atan v1l v1t)) ]
  ask other-particle [
    set speed sqrt ((v2t ^ 2) + (v2l ^ 2))
    set energy (0.5 * mass * (speed ^ 2))
    if v2l != 0 or v2t != 0
      [ set heading (theta - (atan v2l v2t)) ]
  ]
end





;the stick procedure indicates that should a tatA encounter a PP1 molecule which awaits transporting and; 
;does not have enough tatA molecules, it should stop and "stick to the PP1;
to stick
  let candidate one-of other PP1s-here
  if candidate != nobody [
    ask candidate [
      if bound-tatA < tatA-required and transport = true [
        ask tatAs-here [
          set transport true]
        set bound-tatA (bound-tatA + 1)
      ]
    ]
  ]
end


;the translocation procedure causes a PP1 molecule which has enough tatA molecules to move horizontally through the gate;
;upon reaching the periplasm, the tatA molecules used for transport are released for re-use and the PP1 molecule allowed to move freely once more;
to translocation
  if transport = true [
  set heading 270 
  fd 1
  ask tatAs in-radius 3 with [transport = true] [ set heading 270
    fd 1 
    if pxcor = 15 [ 
      set transport false
      rt random 180
      ]
    
  ]
  if pxcor < 13 [ 
    fd 1
    rt (360 - random-float 180)
    set transport false
    reset-gate
    ]
    
  ]
end


;a gate used to transport a PP1 molecule must be reset and made available to future PP1 molecules;
;the reset-gate script carries out this resetting of the gates;
to reset-gate 
    ask patches in-radius 5 with [pcolor = 129] [
      set pcolor 78]
    repeat (gate-size * 12) [
     ask patches  with [pcolor = 78] [
       ask patches in-radius 1 with [pcolor = 129] [
          set pcolor 78]]]
end
@#$#@#$#@
GRAPHICS-WINDOW
347
10
862
546
50
50
5.0
1
10
1
1
1
0
1
1
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
5
10
97
43
setup
setup false\n;; false = reds all on\n;; one side
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
101
10
170
43
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
866
63
948
108
PP1s-right
count PP1s with\n[xcor >= 0]
3
1
11

MONITOR
867
15
945
60
PP1-left
count PP1s with\n[xcor < 0]
3
1
11

PLOT
866
260
1116
442
Progression Graph
Time
Number
0.0
50.0
-100.0
100.0
true
true
"" ""
PENS
"PP1" 1.0 0 -2674135 true "" "let PP1s-left count PP1s with [xcor < 0]\nlet PP1s-right count PP1s with [xcor >= 0]\nplot (PP1s-right + PP1s-left)"
"mc" 1.0 0 -13345367 true "" "let mcs-left count mcs with [xcor < 0]\nlet mcs-right count mcs with [xcor >= 0]\nplot (mcs-right + mcs-left)"
"midline" 1.0 0 -16777216 false ";; draw the x axis in black\nauto-plot-off\nset-current-plot-pen \"midline\"\nplot 0\nplotxy 10000 0\nauto-plot-on" ""

SLIDER
5
308
177
341
gate-number
gate-number
1
17
5
1
1
NIL
HORIZONTAL

SLIDER
4
161
176
194
PP1-production
PP1-production
0
50
5
1
1
NIL
HORIZONTAL

SLIDER
4
124
176
157
degP-number
degP-number
0
42
3
1
1
NIL
HORIZONTAL

MONITOR
867
209
947
254
deg-number
deg-number
17
1
11

MONITOR
866
111
970
156
PP1-transported
count PP1s with [xcor < 0] + count complexes with [xcor < 0]
17
1
11

SLIDER
5
50
177
83
Initial-PP1
Initial-PP1
1
100
61
1
1
NIL
HORIZONTAL

SWITCH
6
458
109
491
collide?
collide?
1
1
-1000

SWITCH
6
495
138
528
random-walk?
random-walk?
1
1
-1000

SLIDER
4
87
176
120
mc-number
mc-number
0
100
56
1
1
NIL
HORIZONTAL

MONITOR
866
161
980
206
Bound PP1
count complexes
17
1
11

SLIDER
5
344
177
377
gate-size
gate-size
0
10
4
1
1
NIL
HORIZONTAL

SLIDER
4
197
176
230
mc-production
mc-production
0
100
7
1
1
NIL
HORIZONTAL

INPUTBOX
5
538
65
598
secBprob
100
1
0
Number

INPUTBOX
70
539
130
599
gateprob
100
1
0
Number

INPUTBOX
135
539
193
599
degprob
100
1
0
Number

INPUTBOX
198
539
257
599
bindprob
100
1
0
Number

TEXTBOX
196
10
346
38
1. Press 'setup'\n2. Press 'go'
11
0.0
1

TEXTBOX
194
59
344
77
Initial number PP1
11
0.0
1

TEXTBOX
195
95
345
113
Initial number microcystin
11
0.0
1

TEXTBOX
195
133
345
161
Initial number degradation\nproteins
11
0.0
1

TEXTBOX
194
167
344
185
Number PP1 produced per tick
11
0.0
1

TEXTBOX
195
202
345
230
Number microcystin produced\nper tick
11
0.0
1

TEXTBOX
195
315
345
333
Number of tatB-C complexes
11
0.0
1

TEXTBOX
195
344
345
372
Width of tatB-C complexes (pixels)
11
0.0
1

TEXTBOX
126
460
321
488
Particles collide off one another
11
0.0
1

TEXTBOX
157
498
340
526
Particles undergo brownians motion
11
0.0
1

TEXTBOX
15
605
67
651
% secB binds to PP1
11
0.0
1

TEXTBOX
76
603
133
668
% PP1-secB binds to gate
11
0.0
1

TEXTBOX
141
602
202
655
% degradation occurs
11
0.0
1

TEXTBOX
208
602
266
649
% PP1 binds to microcystin
11
0.0
1

TEXTBOX
958
16
1108
34
Number PP1 in periplasm
11
0.0
1

TEXTBOX
965
64
1115
82
Number PP1 in cytoplasm
11
0.0
1

TEXTBOX
986
111
1136
139
Number PP1 gone through tatB-C complexes
11
0.0
1

TEXTBOX
996
162
1146
190
Number complexes
11
0.0
1

TEXTBOX
961
209
1111
237
Number of degradation proteins
11
0.0
1

TEXTBOX
924
445
1074
473
Graph showing number PP1 versus number Microcystin
11
0.0
1

SLIDER
5
416
177
449
pore-size
pore-size
0
10
5
1
1
NIL
HORIZONTAL

TEXTBOX
196
417
346
445
Widths of the pores in the outer-membrane
11
0.0
1

SLIDER
5
380
177
413
pore-number
pore-number
0
15
6
1
1
NIL
HORIZONTAL

TEXTBOX
197
380
347
408
Number of pores in the outer-membrane
11
0.0
1

SLIDER
4
234
176
267
tatA-number
tatA-number
0
100
56
1
1
NIL
HORIZONTAL

SLIDER
5
270
177
303
tatA-required
tatA-required
0
30
2
1
1
NIL
HORIZONTAL

TEXTBOX
196
239
346
257
Initial number tatA proteins
11
0.0
1

TEXTBOX
198
268
348
296
Number tatA needed to export PP1 through tatB-C complex
11
0.0
1

@#$#@#$#@
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

key
false
0
Rectangle -7500403 true true 90 120 285 150
Rectangle -7500403 true true 255 135 285 195
Rectangle -7500403 true true 180 135 210 195
Circle -7500403 true true 0 60 150
Circle -16777216 true false 30 90 90

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.0.4
@#$#@#$#@
setup false
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
