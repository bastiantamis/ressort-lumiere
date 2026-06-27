open Graphics

(*Manipulation de la matrice de tableaux*)
let make_espace a b c x =
  Array.init a (fun i -> Array.make_matrix b c x) 

let copy_espace arr =
  Array.init (Array.length arr) (fun x ->
    Array.init (Array.length arr.(0)) (fun y ->
      Array.init (Array.length arr.(0).(0)) (fun ch -> arr.(x).(y).(ch))))

(*Variables*)
let w, h, channels, cell = 201, 201, 3, 3

let dt, amp, min_damp = 1., 1000., 1.
let coef = Array.init channels (fun i -> 1.4 +. float_of_int i /. 30.)
let z = make_espace w h channels 0.
let v = make_espace w h channels 0.
let k = Array.make_matrix w h 1.

let sqr x = x * x

let gaussienne i j r = 
    exp (-. 8. *. float_of_int (sqr i + sqr j) /. float_of_int (sqr r))

let pulsion_blanche_en ox oy r =
  for i = -r to r do 
      for j = -r to r do 
        for ch = 0 to channels - 1 do
          let (x, y) = (i + ox, j + oy) in
          if 0 <= x && x < w && 0 <= y && y < h then
            (*let fi, fj, fr = float_of_int i, float_of_int j, float_of_int r in
            let f = amp *. gaussienne i j r *. sin (fi *. 2. /. 3.14) *. cos (fj /. fr /. 3.14) in*)
            let f = gaussienne i j r *. amp in
            z.(x).(y).(ch) <- z.(x).(y).(ch) +. f (*if f > 0. then  z.(x).(y).(ch) +. f else 0.*);
            v.(x).(y).(ch) <- v.(x).(y).(ch) +. f
        done 
      done 
    done

(*Mise en place de l'état initial*)
let () =
  pulsion_blanche_en (w / 4) (h / 2) 8;
  
  let r = w / 4 in
  let ox, oy = (w / 2, h / 2) in
  (*Disque en verre*)
  for i = -r to r do
    for j = -r to r do
      let (x, y) = (i + ox, j + oy) in
      if sqr i + sqr j < sqr r then k.(x).(y) <- 1.5
    done
  done

  (*Prisme en verre
  for i = -r to r do
    for j = -r to r do
      let (x, y) = (i + ox, j + oy) in
      if abs i - r / 2 < j / 2 then k.(x).(y) <- 1.5 
    done
  done*)

(*Boucle principale*)
let prev_max_amp = ref amp

let avgs = Array.make_matrix w h 0.

let loop () =
  let z0 = copy_espace z in
  let max_amp = ref 0. in 
  
  let safe_get arr x y ch = 
    try arr.(x).(y).(ch) with | Invalid_argument _ -> 0. 
  in

  (*Pour chaque cellule*)
  for x = 0 to w - 1 do
    for y = 0 to h - 1 do
      
      let moyenne_voisinage ch =
        safe_get z0 x y ch -. 0.25 *.
          (safe_get z0 (x - 1) y ch 
          +. safe_get z0 (x + 1) y ch 
          +. safe_get z0 x (y + 1) ch 
          +. safe_get z0 x (y - 1) ch)
      in

      let absf f = if f > 0. then f else -.f in
      
      (*Pour chaque chaîne de couleur*)
      let avg = ref 0. in
      
      for ch = 0 to channels -1 do
        if x mod (w - 1) = 0 || y mod (h - 1) = 0 then z.(x).(y).(ch) <- z.(x).(y).(ch) *. -0.5
        else
          (
          let time_step_add arr arr0 ch f = 
            (*if absf f >= min_damp then*) arr.(x).(y).(ch) <- safe_get arr0 x y ch +. dt *. f 
            in

          (*Formule de la force du rappel d'un ressort: a = -k(l - l0)*)
          let dv = -. moyenne_voisinage ch *. coef.(ch) /. k.(x).(y) in
          time_step_add v v ch dv;
          
          let damp = safe_get v x y ch in
          time_step_add z z0 ch damp;
          avg := !avg +. absf z.(x).(y).(ch) /. float_of_int channels
          )
      done;

      avgs.(x).(y) <- !avg;
      if absf !avg > !max_amp then max_amp := absf !avg;
      
      (*Affichage*)
      let col ch = 
        min 255 (int_of_float ((k.(x).(y) -. 1.) *. 50. +. absf (safe_get z x y ch) /. !prev_max_amp *. 255.))
      in

      set_color (rgb (col 0) (col 1) (col 2));
      fill_rect (x * cell) (y * cell) cell cell;
    done
  done;
  prev_max_amp := !max_amp

(*Boucle principale*)
let () =
  open_graph ""; resize_window (w * cell) (h * cell); auto_synchronize false;
  
  set_color black;
  for x = 0 to w - 1 do for y = 0 to h - 1 do
    fill_rect (x * cell) (y * cell) cell cell;
  done done;

  while not (key_pressed () && read_key () = ' ') do
    loop ();
    synchronize ()
  done;
  close_graph ()