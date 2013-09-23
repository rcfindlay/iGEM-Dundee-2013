;creation of populations;
breed [PP1s PP1]                                                           ;; red PP1 molecules
breed [degPs degP]                                                         ;; degradation machinery - black squares
breed [mcs mc]                                                             ;; microcystin, purple turtles, assuming so small that kinetic interaction is negligible
breed [complexes complex]                                                  ;; PP1 bounded to microcystin, purple crosses
breed [secBs secB]                                                         ;; secB proteins, green keys, which bind to microcystin allowing them to enter sec gate



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
  linked?                                                                  ;; note of whether or not a secB protein is bound to the turtle
]













;creates a class called particles for calling-in PP1s, degPs, mcs & complexes;
to-report particles
  report (turtle-set PP1s degPs mcs complexes)
end



;creates a separate class for particles which exist in the periplasm;
to-report periparticles
  report (turtle-set degPs mcs complexes)
end













;setup for the characteristics of PP1 proteins;
to setup-PP1
  set mass 37500                                                           ;; the masses used was the approximate mass of PP1 (all masses are in approximations of daltons)
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  
  ;physical characteristics;
  set color red
  set size 1.5
  
  set linked? false                                                        ;; linked is a variable used to describe whether a PP1 is bound to a secB protein - All PP1 molecules begin unbound
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


;setup for the characteristics for the secB proteins;
to setup-secB
  set speed 10
  set mass 17278                                                           ;; secBs have a referenced mass of 17278
  set energy (0.5 * mass * (speed ^ 2))
  set last-collision nobody
  set color 55
  set size 2
  set linked? false                                                        ;; linked, again, defines whether the secB protein is bound or free to become bound to a PP1 molecule
end











;The setup procedure creates our environment and sets some conditions for the simulation;
to setup [both-sides?]
  clear-all
  
  ;sets the shapes of particles;
  set-default-shape PP1s "circle"                                          ;; PP1 is depicted as a red circle 
  set-default-shape degPs "square"                                         ;; degradation machinery is depicted as a black square
  set-default-shape mcs "turtle"                                           ;; microcystin is depicted as a purple turtle
  set-default-shape secBs "key"                                            ;; secB proteins are depicted as green keys
  
  
  
  
  ;sets up the initial conditions;
  
  ;initial loop values;
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
  

  
  ;creates an initial number of PP1 molecules and secB molecules with random positions in the cytoplasm;
  create-PP1s initial-PP1
  [setup-PP1
    randomize-right ]
  
  create-secBs initial-secB
  [setup-secB
    randomize-right ]
    
    
    
    
    
  ;creates an initial number of degradation machineries and microcystin molecules with random positions in the periplasm;  
  create-degPs degP-number
  [setup-degPs
    randomize-left ]
    
  create-mcs mc-number
  [ setup-mcs
    randomize-left]
    
 
    
  reset-ticks
end













;position randomizing procedure for periplasmic proteins;
to randomize-left
  setxy (- (abs random-xcor))
  random-ycor
  if any? patches in-radius 1 with [pcolor != 87]
    [ randomize-left ] ;;try again until we don't land on or near the edges of the periplasm
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
 
  
  
  ;due to the binding nature of secBs, they require separate coding for their independant vs bound motion;
  ask secBs [
    stick                                                                  ;; calls a procedure to deal with the binding of secB proteins to PP1 molecules in the cytoplasm and the revoking of secB independant motion
    if collide? [check-for-collision]                                      ;; checks and acts if collisions have been activated in the interface (used for all particle types)
    if random-walk? [rt random-float 360]                                  ;; checks and acts if random (brownian) motion has been activated for particle motion
    bounce                                                                 ;; procedure to allow particles to bounce from walls (if required) and the membrane 
    if linked? = false [fd 1]                                              ;; this statement allows secB proteins to carry out independant motion ONLY if they are unbound (i.e. if they are not linked)
  ]
  
  
  
  
  ;for each secB which is used to transport a PP1, we kill these and re-create that same number at random gate positions - this is for ease of coding and;
  ;for simulation of the lack of knowledge on how long these proteins would be bound at specific gates;
  create-secBs (initial-secB - (count secBs)) [
    setup-secB
    setxy 1 (one-of [random-pxcor] of patches with [pcolor = 78])
  ]
  
  
  
  degrade                                                                  ;; a procedure to degrade PP1s and complexes which encounter the degradation machinery  
  
  
  
  ;for remaining memebers of particles class, carry out the usual checks for collision and motion before carrying out any required bounces and then progress
  ask particles [
    if collide? [check-for-collision]
    if random-walk? [rt random-float 360]
    bounce
    fd 1
  ]
  
  
  
  ;the simulation continues forward one tick;
  tick
end


















;;;;Procedures to create the walls and gates of chosen width;;;;

;procedure used to create the environment with coloured areas representing the cytoplasm and membrane;
to walls
  ask patches with [pxcor > 0 and pxcor < 50][ set pcolor 19.9]             ;; colour the cytoplasm (right half) white                      
  ask patches with [ pxcor = (- 50)] [ set pcolor 3 ]                       ;; colour the outermembrane (left) dark grey     
  ask patches with [ pxcor = 0] [ set pcolor 75 ]                           ;; colour the membrane (middle) green      
  ask patches with [ pxcor = 50 ] [ set pcolor 8 ]                          ;; colour the region indicating the rest of the cell (right) light grey
  ask patches with [ pxcor > -50 and pxcor < 0] [ set pcolor 87 ]           ;; colour the periplasm (left side) light blue
end


;procedure to create a specified number of gates;
to gates
 
  ;for an odd number of gates;
  if (gate-number mod 2 ) = 1 [                                             
    set seperation ( round ( 99 / ( (gate-number + 1) * 2 ) ) )            ;; calculate the required separation between consecutive gate centres (odd number of gates)
    if j != ((gate-number + 1) / 2) [                                      ;; run through the pairs of gates (above and below the centre) using j as our counting variable and 
                                                                           ;; only stop once all gates have been considered
      
      ask patches with [ pxcor = 0 and 
        (pycor = 0 + (seperation * j * 2) or 
          pycor = 0 - (seperation * j * 2) )] [                            ;; ask patches in the centre of our current gate pair (j)
      set pcolor 78                                                        ;; colour these patches a light green
      set g 0                                                              ;; set the gate-size count of the current gate pair (j)  = 0              
      gate-width                                                           ;; call the gate-width procedure to colour the width of the gate this same light green colour
          ]
      
      set j (j + 1)                                                        ;; progress to the next set of gate pairs
     gates                                                                 ;; re-run gate until all gate pairs have been considered
    ]
  ]
  
  
  ;for an even number of gates;
  if (gate-number mod 2) = 0 [
    set seperation ( round ( 99 / ( gate-number * 2 ) ) )                  ;; calculate the required separation between consecutive gate centres (even number of gates)
    if j != (gate-number / 2) [                                            ;; run through the pairs of gates (above and below the centre) using j as our counting variable and
                                                                           ;; only stop once all gates have been considered
      ask patches with [ pxcor = 0 and 
        ( pycor = ( - seperation) + (seperation * (j + 1) * 2) or       
          pycor = seperation - (seperation * (j + 1) * 2) )] [             ;; ask patches in the centre of our current gate pair (j)
      set pcolor 78                                                        ;; colour these patches a light green
      set g 0                                                              ;; set the gate-size count of the current gate pair (j)  = 0
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
    
   













;;;;Procedures for protein/particle death and loss;;;;


; procedure to kill PP1 proteins which disappear back into the cell;
to loss-PP1
  if [pcolor] of patch-ahead 1 = 8 [                                       ;; if the patch ahead is that which goes to the the inner cell (right wall - colour light grey))
    if linked? = false [die                                                ;; if the PP1 is not linked to a secB, it returns to the cell and dissapears (dies)
    ]
  ]
end




; procedure to kill microcystin which disappear back out into the environment;
to loss-mc 
  if [pcolor] of patch-ahead 2 = 3 [die]                                   ;; if the patch ahead is that of the outer membrane (dark grey) microcystin escape and so die
end




; procedure to degrade proteins in the periplasm using degradation proteins;
to degrade
  
  ask degPs [ let prey one-of PP1s-here                                    ;; ask degradation machinery to target PP1s on their patch
    if prey != nobody [                                                    ;; if a PP1 exists on their patch
      ask PP1s-here [ 
        
        if random 100 <= degprob [                                         ;; depending upon the degradation probability defined in the interface
          ask prey [die]                                                   ;; have the PP1 degrade and die
        ]
      ]set deg-number ( deg-number + 1 )]                                  ;; increase the degraded counter by one
  ]
  
  
  ask degPs [ let prey one-of complexes-here                               ;; similarly ask the degradation machinery to target complexes on their patch
    if prey != nobody [                                                    ;; if there is such a complex
      ask complexes-here [ 
        if random 100 <= degprob [                                         ;; depending on the degradation probability
          ask prey [die]                                                   ;; have the complexes become degraded and so die
        ]
      ]set deg-number ( deg-number + 1 )]                                  ;; increase the degraded counter by one
  ]
  
  
end















;procedure dealing with the binding of microcystin to PP1 in the periplasm
to bind
  
  if count PP1s-here = 1 [                                                 ;; ask secBs to check if there is a single PP1 on their current position
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
  
  ;turtles should always bounce off the solid parts of the inner and outer membrane
  if [pcolor] of patch-ahead 1 = 75 or [pcolor] of patch-ahead 1 = 3 [
    set heading (- heading) 
  ]
  
  
  ;turtles should always bounce off the gates if they are currently in the periplasm - no return
  if [pcolor] of patch-ahead 1 = 78 and [pcolor] of patch-here = 87 [
    set heading (- heading)
  ]
  
  
  ;turtles of periparticles class only belong in the periplasm and so if they encounter any patches not of the periplasm colour, they should reflect
  ask periparticles [
    if [pcolor] of patch-ahead 1 != 87 [set heading (- heading)]
  ]
  
  
  ;for secBs specifically, they should reflect from the periplasm (staying in the cytoplasm) and we simulate them also reflecting when approaching the inner cell
  ;we can use this approximation as on average, an equal number will migrate to the inner section of the cell as migrate into the outer part
  ask secBs [
    if [pcolor] of patch-ahead 1 = 87 or [pcolor] of patch-ahead 1 = 8  [set heading (- heading)]
  ]
  
  
  ;for PP1s, they should only be considered to continue into the periplasm if they are bound to secB
  ask PP1s [
    if [pcolor] of patch-ahead 1 = 78 and [pcolor] of patch-here = 19.9 [
      if linked? = false [                                                 ;; if the PP1 is not bound to a secB protein         
        set heading (- heading)]                                           ;; the PP1 must reflect from the gate
      if linked? = true [                                                  ;; however, if the PP1 is bound to a secB it should be allowed to transfer
        if gateprob <= (random 100) [                                      ;; But depending upon the gate transfer probability
          set heading (- heading)                                          ;; it may still be reflected
        ]
      ]
    ]
    
    
    
    if [pcolor] of patch-ahead 1 = 87 and [pcolor] of patch-here = 78 [    ;; if the PP1 does continue on - i.e. it is allowed by the probability
      
      set linked? false                                                    ;; unlink the PP1 from it's secB
      unstick                                                              ;; unstick the two proteins
                                                                           ;; the PP1 direction is unchanged so it is transported
    ]
  ]
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






; a procedure for the attachment of secBs to PP1 molecules in the cytoplasm
to stick
  if (count other PP1s-here = 1 ) and (linked? = false) [                  ;; ask secBs if there are PP1s in their current positions ensure that the secBs are not linked     
    if random 100 <= secBprob [                                            ;; if this is true and probability assigned in the interface agrees, continue
      let candidate one-of other PP1s-here                                 ;; assign the PP1 that is here as a candidate
      if (candidate != nobody)                                             ;; as long as the candidate exists, continue
      [
        ask candidate [
          if linked? = false [                                             ;; ensure that the candidate is also not linked and if so continue
            create-link-to myself[                                         ;; create a link with the candidate
              tie                                                          ;; create a tie with the candidate
              hide-link                                                    ;; hide the created (visual) link
            ] ask one-of other secBs-here [set linked? true]               ;; set the secBs linked? property to true (now that it has linked)
            set linked? true                                               ;; set the PP1s linked? property to true (again as it is linked)
          ]
        ]
        
      ]
    ]
  ]
end




; a procedure for the unattachment of secBs on transport through the gates
to unstick
  ask out-link-neighbors [                                                 ;; ask the linked neighbour
    die                                                                    ;; die - this also removes the link, leaving the PP1 free to progress
  ]
  
end
@#$#@#$#@
GRAPHICS-WINDOW
346
10
861
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
306
1116
488
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
3
269
175
302
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
5
196
177
229
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
23
1
1
NIL
HORIZONTAL

MONITOR
866
257
946
302
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
56
1
1
NIL
HORIZONTAL

SWITCH
1
341
104
374
collide?
collide?
0
1
-1000

SWITCH
1
378
133
411
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
54
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
0
304
177
337
gate-size
gate-size
0
10
7
1
1
NIL
HORIZONTAL

SLIDER
5
232
177
265
mc-production
mc-production
0
100
7
1
1
NIL
HORIZONTAL

SLIDER
4
160
176
193
initial-secB
initial-secB
0
100
20
1
1
NIL
HORIZONTAL

MONITOR
865
210
922
255
secBs
count secBs
17
1
11

INPUTBOX
1
416
61
476
secBprob
100
1
0
Number

INPUTBOX
71
416
131
476
gateprob
100
1
0
Number

INPUTBOX
141
417
199
477
degprob
100
1
0
Number

INPUTBOX
206
417
265
477
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
195
169
345
187
Initial number secB proteins
11
0.0
1

TEXTBOX
195
202
345
220
Number PP1 produced per tick
11
0.0
1

TEXTBOX
196
237
346
265
Number microcystin produced\nper tick
11
0.0
1

TEXTBOX
194
276
344
294
Number of sec gates
11
0.0
1

TEXTBOX
194
305
344
323
Width of sec gates (pixels)
11
0.0
1

TEXTBOX
121
343
316
371
Particles collide off one another
11
0.0
1

TEXTBOX
152
381
335
409
Particles undergo brownians motion
11
0.0
1

TEXTBOX
3
478
55
524
% secB binds to PP1
11
0.0
1

TEXTBOX
71
480
128
545
% PP1-secB binds to gate
11
0.0
1

TEXTBOX
142
486
203
539
% degradation occurs
11
0.0
1

TEXTBOX
211
485
269
532
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
Number PP1 gone through gates
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
943
212
1093
230
Number of secB proteins
11
0.0
1

TEXTBOX
960
257
1110
285
Number of degradation proteins
11
0.0
1

TEXTBOX
924
491
1074
519
Graph showing number PP1 versus number Microcystin
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model follows on the Simple Kinetics 2 model.  In Simple Kinetics 2, we saw how changes to variables such as temperature, volume, and concentration affected the rate at which a chemical reaction reached an equilibrium state.  Here we model the same phenomenon based upon a physical separation.

This model investigates some laboratory methods that chemists use to purify chemicals.  Most of these methods are based upon physical properties of molecular separation. The same principles that affect chemical equilibrium affect physical equilibrium as well. By changing the variables of temperature, volume, and concentration, we can affect not only the speed at which a system reaches equilibrium, but also the nature of the distribution ratio. In this model, we watch how these factors affect the physical distribution of red molecules that are considered "dissolved" in a blue solvent.

## HOW TO USE IT

Setup the model by pressing either the SETUP-RANDOM or the SETUP-SIDE buttons at the top of the Interface tab.  SETUP-RANDOM distributes all the molecules randomly around the world. SETUP-SIDE distributes the blue molecules evenly, while placing the red molecules on the right side of the world.

Press GO to watch the molecules move about the world as they achieve equilibrium. The plot tracks the relative concentrations of each color molecule on each side of the central divider. If the red line dips below 0, there are more red molecules on the left side of the divider than on the right. If it rises above 0, there are more red molecules on the right side of the divider than on the left. The blue line plots the same relationship for blue molecules.

You can add more red molecules to the right side of the world by pressing ADD RED.

Similarly, you can shrink or expand the right side of the box with the buttons SHRINK RIGHT and EXPAND RIGHT, respectively.

Finally, to change the size of the connection window, move the WINDOW slider to your desired size and then press the CHANGE WINDOW button.

## THINGS TO NOTICE

Pay attention to the plot and compare it what you see in the world. Is there an equal number of blue and red molecules on each side of the divider according to the plot and according to what you see in the view?

## THINGS TO TRY

Run the model with several different states for each variable. Do you observe similar equilibrium effects to those seen in Simple Kinetics 2?  Are there significant differences?

Does the temperature affect the system in the same way it affected the chemical reaction in Simple Kinetics 2?  Why or why not?

How does changing the concentration affect the rate at which the molecules achieve equilibrium? Does this make sense?

## EXTENDING THE MODEL

The system we have established here always comes to an approximately identical equilibrium state, no matter how you change the variables. In the lab, this is not useful to chemists, who want to separate one type of molecule from another. Can you extend the model to separate all of the red molecules from the blue molecules?

Try adding another color of molecule to the system and randomly distributing all the molecules around the world. Can you devise a way to separate the new molecules from the red molecules?

Add a slider that allows you to alter the temperature of the system. Think about what effect cooling and heating the system would have on the molecules. Be sure to include a command in the procedures window that will execute your proposed effect.

## RELATED MODELS

Simple Kinetics 1, Simple Kinetics 2

## RELATED MODELS

Simple Kinetics 1
Simple Kinetics 2

## CREDITS AND REFERENCES

Thanks to Mike Stieff for his work on this model.


## HOW TO CITE

If you mention this model in a publication, we ask that you include these citations for the model itself and for the NetLogo software:

* Stieff, M. and Wilensky, U. (2001).  NetLogo Simple Kinetics 3 model.  http://ccl.northwestern.edu/netlogo/models/SimpleKinetics3.  Center for Connected Learning and Computer-Based Modeling, Northwestern Institute on Complex Systems, Northwestern University, Evanston, IL.
* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern Institute on Complex Systems, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2001 Uri Wilensky.

![CC BY-NC-SA 3.0](http://i.creativecommons.org/l/by-nc-sa/3.0/88x31.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

This model was created as part of the projects: PARTICIPATORY SIMULATIONS: NETWORK-BASED DESIGN FOR SYSTEMS LEARNING IN CLASSROOMS and/or INTEGRATED SIMULATION AND MODELING ENVIRONMENT. The project gratefully acknowledges the support of the National Science Foundation (REPP & ROLE programs) -- grant numbers REC #9814682 and REC-0126227.
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
